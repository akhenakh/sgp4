package sgp4

import (
	"encoding/json"
	"fmt"
	"strings"
	"time"
)

// OMM represents a single Orbit Mean-elements Message object from a JSON representation.
// Fields are based on the CCSDS OMM standard and common JSON outputs (e.g., from space-track.org).
type OMM struct {
	ObjectName         string  `json:"OBJECT_NAME"`
	ObjectID           string  `json:"OBJECT_ID"`   // e.g., "1998-067A"
	EpochStr           string  `json:"EPOCH"`       // ISO 8601 e.g., "2025-05-26T13:06:57.824640"
	MeanMotion         float64 `json:"MEAN_MOTION"` // rev/day
	Eccentricity       float64 `json:"ECCENTRICITY"`
	Inclination        float64 `json:"INCLINATION"`         // degrees
	RAOfAscNode        float64 `json:"RA_OF_ASC_NODE"`      // degrees (Right Ascension of Ascending Node)
	ArgOfPericenter    float64 `json:"ARG_OF_PERICENTER"`   // degrees
	MeanAnomaly        float64 `json:"MEAN_ANOMALY"`        // degrees
	EphemerisType      int     `json:"EPHEMERIS_TYPE"`      // Typically 0 for TLE-derived SGP4 elements
	ClassificationType string  `json:"CLASSIFICATION_TYPE"` // e.g., "U" for unclassified
	NoradCatID         int     `json:"NORAD_CAT_ID"`
	ElementSetNo       int     `json:"ELEMENT_SET_NO"`   // Element Set Number for TLE
	RevAtEpoch         int     `json:"REV_AT_EPOCH"`     // Revolution number at epoch for TLE
	BStar              float64 `json:"BSTAR"`            // B* drag term, 1/EarthRadii
	MeanMotionDot      float64 `json:"MEAN_MOTION_DOT"`  // First time derivative of mean motion (n-dot/2 for TLE), rev/day^2
	MeanMotionDDot     float64 `json:"MEAN_MOTION_DDOT"` // Second time derivative of mean motion (n-ddot/6 for TLE), rev/day^3

	// Optional fields that might appear in more complete CCSDS OMM JSON representations
	CenterName        string `json:"CENTER_NAME,omitempty"`
	RefFrame          string `json:"REF_FRAME,omitempty"`
	TimeSystem        string `json:"TIME_SYSTEM,omitempty"`
	MeanElementTheory string `json:"MEAN_ELEMENT_THEORY,omitempty"`
}

// ParseOMMs parses a JSON byte slice containing an array of OMM objects.
func ParseOMMs(jsonData []byte) ([]OMM, error) {
	var omms []OMM
	err := json.Unmarshal(jsonData, &omms)
	if err != nil {
		return nil, fmt.Errorf("error unmarshalling OMM JSON: %w", err)
	}
	return omms, nil
}

// ommEpochToTleEpoch converts an OMM epoch string (ISO 8601) to TLE epoch components.
// It returns the full year, the day of the year (1.0 to 366.xxxx), the parsed time.Time (UTC), and any error.
func ommEpochToTleEpoch(epochStr string) (epochYear int, epochDayFloat float64, epochTimeUTC time.Time, err error) {
	var t time.Time
	var parseErr error

	effectiveEpochStr := epochStr
	hasOriginalZ := strings.HasSuffix(epochStr, "Z")
	hasOriginalOffset := false
	if !hasOriginalZ {
		lastPlus := strings.LastIndex(epochStr, "+")
		lastMinus := strings.LastIndex(epochStr, "-")
		if lastPlus > 7 && strings.Contains(epochStr[lastPlus:], ":") {
			hasOriginalOffset = true
		} else if lastMinus > 7 && strings.Contains(epochStr[lastMinus:], ":") {
			hasOriginalOffset = true
		}
		if !hasOriginalOffset {
			effectiveEpochStr = epochStr + "Z" // Assume UTC if no explicit Z or offset
		}
	}

	layouts := []string{
		time.RFC3339Nano,
		time.RFC3339,
		"2006-01-02T15:04:05.999999999Z", // For fixed Z
		"2006-01-02T15:04:05Z",           // For fixed Z
	}

	for _, layout := range layouts {
		t, parseErr = time.Parse(layout, effectiveEpochStr)
		if parseErr == nil {
			break
		}
	}

	// If parsing with explicit/appended Z failed, and original string had neither Z nor offset,
	// try parsing the original string "as-is", interpreting it as UTC.
	if parseErr != nil && !hasOriginalZ && !hasOriginalOffset {
		plainLayouts := []string{
			"2006-01-02T15:04:05.999999999",
			"2006-01-02T15:04:05.999999",
			"2006-01-02T15:04:05.999",
			"2006-01-02T15:04:05",
		}
		for _, layout := range plainLayouts {
			tempTime, plainParseErr := time.ParseInLocation(layout, epochStr, time.UTC)
			if plainParseErr == nil {
				t = tempTime
				parseErr = nil
				break
			}
		}
	}

	if parseErr != nil {
		return 0, 0, time.Time{}, fmt.Errorf("error parsing OMM epoch string '%s' (tried as '%s'): %w", epochStr, effectiveEpochStr, parseErr)
	}

	epochTimeUTC = t.In(time.UTC)
	epochYear = epochTimeUTC.Year()

	dayOfYearIntegerPart := float64(epochTimeUTC.YearDay())
	nanosInStandardDay := float64(24 * 60 * 60 * 1e9)
	startOfCurrentDay := time.Date(epochYear, epochTimeUTC.Month(), epochTimeUTC.Day(), 0, 0, 0, 0, time.UTC)
	durationIntoCurrentDayNs := epochTimeUTC.Sub(startOfCurrentDay).Nanoseconds()
	fractionOfDay := float64(durationIntoCurrentDayNs) / nanosInStandardDay
	epochDayFloat = dayOfYearIntegerPart + fractionOfDay

	return epochYear, epochDayFloat, epochTimeUTC, nil
}

// ommObjectIDToTleInternational converts an OMM OBJECT_ID to TLE International Designator.
// OMM format: "YYYY-NNNP{PP}" (e.g., "1998-067A")
// TLE format: "YYNNNP{PP}" (e.g., "98067A")
func ommObjectIDToTleInternational(objectID string) (string, error) {
	parts := strings.Split(objectID, "-")
	if len(parts) != 2 {
		return "", fmt.Errorf("invalid OBJECT_ID format: expected 'YYYY-NNNPPP', got '%s'", objectID)
	}
	yearStr := parts[0]
	launchNumPiece := parts[1]

	if len(yearStr) < 2 { // Should be 4 (e.g., "1998") or at least 2 for "YY"
		return "", fmt.Errorf("invalid year part in OBJECT_ID: '%s'", yearStr)
	}

	// TLE international designator uses last two digits of launch year
	tleLaunchYear := yearStr[len(yearStr)-2:]

	// Basic validation for launchNumPiece (e.g., NNN must be digits, P{PP} must be letters)
	if len(launchNumPiece) < 4 { // Minimum length e.g. "001A"
		return "", fmt.Errorf("invalid launch number/piece part in OBJECT_ID: '%s', too short", launchNumPiece)
	}
	// Example: For 001A, NNN is "001", P is "A".
	// Could add more detailed validation (e.g. NNN part is numeric, piece is alpha) if strictness is needed.

	return tleLaunchYear + launchNumPiece, nil
}

// ToTLE converts an OMM object to a TLE object.
// Note: TLE checksums are not part of OMM data and will be set to 0 in the TLE struct.
// They would typically be calculated when formatting the TLE struct into text lines.
func (o *OMM) ToTLE() (*TLE, error) {
	tle := &TLE{}

	tle.Name = o.ObjectName
	tle.SatelliteNumber = o.NoradCatID

	if len(o.ClassificationType) > 0 {
		tle.Classification = rune(o.ClassificationType[0])
	} else {
		tle.Classification = 'U' // Default if not specified
	}

	var err error
	tle.International, err = ommObjectIDToTleInternational(o.ObjectID)
	if err != nil {
		return nil, fmt.Errorf("failed to convert ObjectID to TLE International: %w", err)
	}

	epochFullYear, epochDayFrac, _, err := ommEpochToTleEpoch(o.EpochStr)
	if err != nil {
		return nil, fmt.Errorf("failed to parse OMM epoch: %w", err)
	}
	tle.EpochYear = epochFullYear // TLE struct stores the full year (e.g., 2024)
	tle.EpochDay = epochDayFrac

	// These fields in OMM are generally equivalent to the TLE fields
	// MeanMotionDot in TLE is n-dot/2; OMM MEAN_MOTION_DOT is assumed to be this value.
	// MeanMotionDDot in TLE is n-ddot/6; OMM MEAN_MOTION_DDOT is assumed to be this value.
	tle.MeanMotionDot = o.MeanMotionDot
	tle.MeanMotionDot2 = o.MeanMotionDDot
	tle.Bstar = o.BStar
	tle.ElementNumber = o.ElementSetNo

	// Checksums are not available from OMM.
	// They are calculated based on the TLE line string format.
	tle.CheckSum1 = 0
	tle.CheckSum2 = 0

	// Line 2 fields
	tle.Inclination = o.Inclination
	tle.RightAscension = o.RAOfAscNode
	tle.Eccentricity = o.Eccentricity
	tle.ArgOfPerigee = o.ArgOfPericenter
	tle.MeanAnomaly = o.MeanAnomaly
	tle.MeanMotion = o.MeanMotion
	tle.RevolutionNumber = o.RevAtEpoch

	// Basic validation of converted TLE fields
	if tle.Eccentricity >= 1.0 || tle.Eccentricity < 0.0 {
		return nil, fmt.Errorf("eccentricity from OMM (%.10f) is out of TLE bounds [0,1)", tle.Eccentricity)
	}
	if tle.Inclination < 0.0 || tle.Inclination > 180.0 {
		return nil, fmt.Errorf("inclination from OMM (%.4f) is out of TLE bounds [0,180]", tle.Inclination)
	}
	// Other TLE field constraints (e.g., range of RAAN, ArgP, Mean Anomaly) are generally 0-360 deg.
	// Mean Motion should be positive.

	return tle, nil
}
