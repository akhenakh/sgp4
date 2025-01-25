package spg4

import (
	"fmt"
	"math"
	"strconv"
	"strings"
	"time"
)

// TLE represents a Two-Line Element set used for satellite tracking
type TLE struct {
	// Line 0 (optional name)
	Name string

	// Line 1 fields
	SatelliteNumber int
	Classification  rune
	International   string // International Designator
	EpochYear       int
	EpochDay        float64
	MeanMotionDot   float64
	MeanMotionDot2  float64
	Bstar           float64
	ElementNumber   int
	CheckSum1       int

	// Line 2 fields
	Inclination      float64
	RightAscension   float64
	Eccentricity     float64
	ArgOfPerigee     float64
	MeanAnomaly      float64
	MeanMotion       float64
	RevolutionNumber int
	CheckSum2        int
}

// EpochTime returns the time.Time representation of the TLE epoch
func (tle *TLE) EpochTime() time.Time {
	year := tle.EpochYear

	wholeDays := int(tle.EpochDay)
	fractionalDay := tle.EpochDay - float64(wholeDays)

	// Convert to hours, minutes, seconds
	hoursInDay := fractionalDay * 24
	hours := int(hoursInDay)
	minutesInHour := (hoursInDay - float64(hours)) * 60
	minutes := int(minutesInHour)
	secondsInMinute := (minutesInHour - float64(minutes)) * 60
	seconds := int(secondsInMinute)
	nanos := int((secondsInMinute - float64(seconds)) * 1e9)

	return time.Date(year, 1, 1, hours, minutes, seconds, nanos, time.UTC).
		Add(time.Duration(wholeDays*24) * time.Hour)
}

// ParseTLE parses a two-line element set string and returns a TLE struct.
// It accepts either a two-line or three-line format (with satellite name).
func ParseTLE(input string) (*TLE, error) {
	lines := strings.Split(strings.TrimSpace(input), "\n")
	if len(lines) < 2 || len(lines) > 3 {
		return nil, fmt.Errorf("invalid TLE: must contain 2 or 3 lines")
	}

	tle := &TLE{}

	// Handle optional name line
	startLine := 0
	if len(lines) == 3 {
		tle.Name = strings.TrimSpace(lines[0])
		startLine = 1
	}

	// Validate line lengths
	line1 := lines[startLine]
	line2 := lines[startLine+1]
	if len(line1) != 69 || len(line2) != 69 {
		return nil, fmt.Errorf("invalid TLE: lines must be 69 characters")
	}

	// Parse Line 1
	var err error
	if err = tle.parseLine1(line1); err != nil {
		return nil, fmt.Errorf("error parsing line 1: %w", err)
	}

	// Parse Line 2
	if err = tle.parseLine2(line2); err != nil {
		return nil, fmt.Errorf("error parsing line 2: %w", err)
	}

	return tle, nil
}

func (tle *TLE) parseLine1(line string) error {
	if line[0] != '1' {
		return fmt.Errorf("line 1 must begin with '1'")
	}

	var err error
	tle.SatelliteNumber, err = strconv.Atoi(strings.TrimSpace(line[2:7]))
	if err != nil {
		return fmt.Errorf("invalid satellite number: %w", err)
	}

	tle.Classification = rune(line[7])
	tle.International = strings.TrimSpace(line[9:17])

	year, err := strconv.Atoi(line[18:20])
	if err != nil {
		return fmt.Errorf("invalid epoch year: %w", err)
	}
	tle.EpochYear = 2000 + year
	if year > 50 {
		tle.EpochYear = 1900 + year
	}

	tle.EpochDay, err = strconv.ParseFloat(line[20:32], 64)
	if err != nil {
		return fmt.Errorf("invalid epoch day: %w", err)
	}

	// Trim spaces before parsing mean motion dot
	tle.MeanMotionDot, err = strconv.ParseFloat(strings.TrimSpace(line[33:43]), 64)
	if err != nil {
		return fmt.Errorf("invalid mean motion dot: %w", err)
	}

	// Parse Mean Motion Dot 2 (decimal point assumed)
	mantissa := strings.TrimSpace(line[44:50])
	if mantissa == "00000" {
		tle.MeanMotionDot2 = 0
	} else {
		exp, err := strconv.ParseFloat(mantissa, 64)
		if err != nil {
			return fmt.Errorf("invalid mean motion dot 2: %w", err)
		}
		tle.MeanMotionDot2 = exp * 1e-5
	}

	// Parse B* (decimal point assumed)
	bstarStr := strings.TrimSpace(line[53:59])
	if bstarStr == "" {
		return fmt.Errorf("invalid B* value: empty field")
	}
	bstar, err := strconv.ParseFloat(bstarStr, 64)
	if err != nil {
		return fmt.Errorf("invalid B* value: %w", err)
	}

	bstarExp := strings.TrimSpace(line[59:61])
	if bstarExp == "" {
		return fmt.Errorf("invalid B* exponent: empty field")
	}
	exp, err := strconv.Atoi(bstarExp)
	if err != nil {
		return fmt.Errorf("invalid B* exponent: %w", err)
	}
	tle.Bstar = (bstar * 1e-5) * math.Pow(10, float64(exp))

	tle.ElementNumber, err = strconv.Atoi(strings.TrimSpace(line[64:68]))
	if err != nil {
		return fmt.Errorf("invalid element number: %w", err)
	}

	checksum, err := strconv.Atoi(line[68:69])
	if err != nil {
		return fmt.Errorf("invalid checksum: %w", err)
	}
	tle.CheckSum1 = checksum

	return nil
}

func (tle *TLE) parseLine2(line string) error {
	if line[0] != '2' {
		return fmt.Errorf("line 2 must begin with '2'")
	}

	var err error
	satNum, err := strconv.Atoi(strings.TrimSpace(line[2:7]))
	if err != nil {
		return fmt.Errorf("invalid satellite number in line 2: %w", err)
	}
	if satNum != tle.SatelliteNumber {
		return fmt.Errorf("satellite numbers do not match between lines")
	}

	tle.Inclination, err = strconv.ParseFloat(strings.TrimSpace(line[8:16]), 64)
	if err != nil {
		return fmt.Errorf("invalid inclination: %w", err)
	}

	tle.RightAscension, err = strconv.ParseFloat(strings.TrimSpace(line[17:25]), 64)
	if err != nil {
		return fmt.Errorf("invalid right ascension: %w", err)
	}

	// Eccentricity (decimal point assumed)
	ecc, err := strconv.ParseFloat(strings.TrimSpace(line[26:33]), 64)
	if err != nil {
		return fmt.Errorf("invalid eccentricity: %w", err)
	}
	tle.Eccentricity = ecc * 1e-7

	tle.ArgOfPerigee, err = strconv.ParseFloat(strings.TrimSpace(line[34:42]), 64)
	if err != nil {
		return fmt.Errorf("invalid argument of perigee: %w", err)
	}

	tle.MeanAnomaly, err = strconv.ParseFloat(strings.TrimSpace(line[43:51]), 64)
	if err != nil {
		return fmt.Errorf("invalid mean anomaly: %w", err)
	}

	tle.MeanMotion, err = strconv.ParseFloat(strings.TrimSpace(line[52:63]), 64)
	if err != nil {
		return fmt.Errorf("invalid mean motion: %w", err)
	}

	tle.RevolutionNumber, err = strconv.Atoi(strings.TrimSpace(line[63:68]))
	if err != nil {
		return fmt.Errorf("invalid revolution number: %w", err)
	}

	checksum, err := strconv.Atoi(line[68:69])
	if err != nil {
		return fmt.Errorf("invalid checksum: %w", err)
	}
	tle.CheckSum2 = checksum

	return nil
}
