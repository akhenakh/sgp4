package sgp4

import (
	"encoding/csv"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// ParseOMMsCSV parses an OMM/General-Perturbations CSV stream into OMM objects.
//
// The stream must start with a header row naming the OMM fields (as produced by
// CelesTrak's gp.php?...&FORMAT=csv): the column order is not significant, and
// unknown columns are ignored. Because the catalog number is a numeric field,
// this format is not limited to the five digits of the legacy TLE layout.
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
	get := func(rec []string, name string) string {
		if i, ok := col[name]; ok && i < len(rec) {
			return strings.TrimSpace(rec[i])
		}
		return ""
	}

	var omms []OMM
	for {
		rec, err := cr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("sgp4: reading OMM CSV: %w", err)
		}
		omms = append(omms, OMM{
			ObjectName:         get(rec, "OBJECT_NAME"),
			ObjectID:           get(rec, "OBJECT_ID"),
			EpochStr:           get(rec, "EPOCH"),
			MeanMotion:         csvFloat(get(rec, "MEAN_MOTION")),
			Eccentricity:       csvFloat(get(rec, "ECCENTRICITY")),
			Inclination:        csvFloat(get(rec, "INCLINATION")),
			RAOfAscNode:        csvFloat(get(rec, "RA_OF_ASC_NODE")),
			ArgOfPericenter:    csvFloat(get(rec, "ARG_OF_PERICENTER")),
			MeanAnomaly:        csvFloat(get(rec, "MEAN_ANOMALY")),
			EphemerisType:      csvInt(get(rec, "EPHEMERIS_TYPE")),
			ClassificationType: get(rec, "CLASSIFICATION_TYPE"),
			NoradCatID:         csvInt(get(rec, "NORAD_CAT_ID")),
			ElementSetNo:       csvInt(get(rec, "ELEMENT_SET_NO")),
			RevAtEpoch:         csvInt(get(rec, "REV_AT_EPOCH")),
			BStar:              csvFloat(get(rec, "BSTAR")),
			MeanMotionDot:      csvFloat(get(rec, "MEAN_MOTION_DOT")),
			MeanMotionDDot:     csvFloat(get(rec, "MEAN_MOTION_DDOT")),
		})
	}
	return omms, nil
}

// ParseOMMsReader parses either an OMM JSON array or an OMM/GP CSV stream,
// choosing by the first non-space byte.
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
	case '[', '{':
		return ParseOMMs([]byte(trimmed))
	default:
		return ParseOMMsCSV(strings.NewReader(trimmed))
	}
}

func csvFloat(s string) float64 {
	if s == "" || strings.EqualFold(s, "null") {
		return 0
	}
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

func csvInt(s string) int {
	if s == "" || strings.EqualFold(s, "null") {
		return 0
	}
	v, _ := strconv.Atoi(s)
	return v
}
