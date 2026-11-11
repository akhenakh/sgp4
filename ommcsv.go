package sgp4

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// requiredOMMCSVFields are the columns a GP CSV record needs to describe an
// element set. Optional fields (BSTAR, the derivatives, the counters) may be
// absent and default to zero; a column that is present but malformed is always
// an error.
var requiredOMMCSVFields = []string{
	"EPOCH", "MEAN_MOTION", "ECCENTRICITY", "INCLINATION",
	"RA_OF_ASC_NODE", "ARG_OF_PERICENTER", "MEAN_ANOMALY", "NORAD_CAT_ID",
}

// ParseOMMsCSV parses an OMM/General-Perturbations CSV stream into OMM objects.
//
// The stream must start with a header row naming the OMM fields (as produced by
// CelesTrak's gp.php?...&FORMAT=csv): the column order is not significant, and
// unknown columns are ignored. Because the catalog number is a numeric field,
// this format is not limited to the five digits of the legacy TLE layout.
//
// A missing required column or a value that does not parse is returned as an
// error naming the row and field, instead of being silently read as zero.
func ParseOMMsCSV(r io.Reader) ([]OMM, error) {
	cr := csv.NewReader(r)
	cr.FieldsPerRecord = -1
	cr.LazyQuotes = true

	header, err := cr.Read()
	if err != nil {
		return nil, fmt.Errorf("sgp4: reading OMM CSV header: %w", err)
	}
	col := make(map[string]int, len(header))
	for i, h := range header {
		col[strings.ToUpper(strings.TrimSpace(h))] = i
	}
	for _, name := range requiredOMMCSVFields {
		if _, ok := col[name]; !ok {
			return nil, fmt.Errorf("sgp4: OMM CSV header is missing the %s column", name)
		}
	}

	// field returns the trimmed cell, or false when the column is absent or the
	// cell is empty (or the JSON-style null CelesTrak sometimes emits).
	field := func(rec []string, name string) (string, bool) {
		i, ok := col[name]
		if !ok || i >= len(rec) {
			return "", false
		}
		v := strings.TrimSpace(rec[i])
		if v == "" || strings.EqualFold(v, "null") {
			return "", false
		}
		return v, true
	}

	var omms []OMM
	for line := 2; ; line++ {
		rec, err := cr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("sgp4: reading OMM CSV row %d: %w", line, err)
		}

		floatField := func(name string, required bool) (float64, error) {
			s, ok := field(rec, name)
			if !ok {
				if required {
					return 0, fmt.Errorf("sgp4: OMM CSV row %d: missing %s", line, name)
				}
				return 0, nil
			}
			v, err := strconv.ParseFloat(s, 64)
			if err != nil {
				return 0, fmt.Errorf("sgp4: OMM CSV row %d: invalid %s %q", line, name, s)
			}
			return v, nil
		}
		intField := func(name string, required bool) (int, error) {
			s, ok := field(rec, name)
			if !ok {
				if required {
					return 0, fmt.Errorf("sgp4: OMM CSV row %d: missing %s", line, name)
				}
				return 0, nil
			}
			v, err := strconv.Atoi(s)
			if err != nil {
				return 0, fmt.Errorf("sgp4: OMM CSV row %d: invalid %s %q", line, name, s)
			}
			return v, nil
		}

		epoch, ok := field(rec, "EPOCH")
		if !ok {
			return nil, fmt.Errorf("sgp4: OMM CSV row %d: missing EPOCH", line)
		}

		var omm OMM
		numFloats := []struct {
			name     string
			required bool
			dst      *float64
		}{
			{"MEAN_MOTION", true, &omm.MeanMotion},
			{"ECCENTRICITY", true, &omm.Eccentricity},
			{"INCLINATION", true, &omm.Inclination},
			{"RA_OF_ASC_NODE", true, &omm.RAOfAscNode},
			{"ARG_OF_PERICENTER", true, &omm.ArgOfPericenter},
			{"MEAN_ANOMALY", true, &omm.MeanAnomaly},
			{"BSTAR", false, &omm.BStar},
			{"MEAN_MOTION_DOT", false, &omm.MeanMotionDot},
			{"MEAN_MOTION_DDOT", false, &omm.MeanMotionDDot},
		}
		for _, f := range numFloats {
			v, err := floatField(f.name, f.required)
			if err != nil {
				return nil, err
			}
			*f.dst = v
		}
		numInts := []struct {
			name     string
			required bool
			dst      *int
		}{
			{"NORAD_CAT_ID", true, &omm.NoradCatID},
			{"EPHEMERIS_TYPE", false, &omm.EphemerisType},
			{"ELEMENT_SET_NO", false, &omm.ElementSetNo},
			{"REV_AT_EPOCH", false, &omm.RevAtEpoch},
		}
		for _, f := range numInts {
			v, err := intField(f.name, f.required)
			if err != nil {
				return nil, err
			}
			*f.dst = v
		}

		omm.ObjectName, _ = field(rec, "OBJECT_NAME")
		omm.ObjectID, _ = field(rec, "OBJECT_ID")
		omm.ClassificationType, _ = field(rec, "CLASSIFICATION_TYPE")
		omm.EpochStr = epoch
		omms = append(omms, omm)
	}
	return omms, nil
}

// ParseOMMsReader parses an OMM JSON array, a single OMM JSON object, or an
// OMM/GP CSV stream, choosing by the first non-space byte.
func ParseOMMsReader(r io.Reader) ([]OMM, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	trimmed := strings.TrimSpace(strings.TrimPrefix(string(data), "\ufeff"))
	if trimmed == "" {
		return nil, fmt.Errorf("sgp4: empty OMM input")
	}
	switch trimmed[0] {
	case '[':
		return ParseOMMs([]byte(trimmed))
	case '{':
		var omm OMM
		if err := json.Unmarshal([]byte(trimmed), &omm); err != nil {
			return nil, fmt.Errorf("sgp4: parsing OMM JSON object: %w", err)
		}
		return []OMM{omm}, nil
	default:
		return ParseOMMsCSV(strings.NewReader(trimmed))
	}
}
