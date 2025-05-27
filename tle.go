package sgp4

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
	days := int(tle.EpochDay) // Integer part of the day of the year
	fractionalDay := tle.EpochDay - float64(days)

	// Base date is Jan 1st of the epoch year
	baseDate := time.Date(year, 1, 1, 0, 0, 0, 0, time.UTC)
	// Add (days - 1) because day 1 means 0 full days passed from Jan 1st 00:00
	epochBaseDay := baseDate.AddDate(0, 0, days-1)

	// Calculate nanoseconds from the fractional day
	// 86400 seconds in a day. 1e9 nanoseconds in a second.
	// totalNanosInDayFloat := fractionalDay * 86400.0 * 1e9
	// Using math.Round helps to get the closest integer nanosecond value,
	// mitigating floating point inaccuracies that might cause off-by-one errors.
	totalNanosInDay := int64(math.Round(fractionalDay * 86400.0 * 1e9))

	return epochBaseDay.Add(time.Duration(totalNanosInDay))
}

// ParseTLE parses a two-line element set string and returns a TLE struct.
// It accepts either a two-line or three-line format (with satellite name).
func ParseTLE(input string) (*TLE, error) {
	lines := strings.Split(strings.TrimSpace(input), "\n")
	for i, line := range lines { // Trim spaces from each line, useful if input has trailing spaces per line
		lines[i] = strings.TrimSpace(line)
	}

	if len(lines) < 2 || len(lines) > 3 {
		return nil, fmt.Errorf("invalid TLE: must contain 2 or 3 lines")
	}

	tle := &TLE{}

	// Handle optional name line
	startLine := 0
	if len(lines) == 3 {
		tle.Name = lines[0] // Already trimmed
		startLine = 1
	}

	// Validate line lengths
	line1 := lines[startLine]
	line2 := lines[startLine+1]
	if len(line1) != 69 {
		return nil, fmt.Errorf("invalid TLE: line 1 must be 69 characters, got %d", len(line1))
	}
	if len(line2) != 69 {
		return nil, fmt.Errorf("invalid TLE: line 2 must be 69 characters, got %d", len(line2))
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

	calculatedCheckSum1, err := calculateChecksum(line1)
	if err != nil {
		return nil, fmt.Errorf("checksum calculation failed for line 1: %w", err)
	}
	if calculatedCheckSum1 != tle.CheckSum1 {
		return nil, fmt.Errorf("checksum mismatch in line 1: expected %d (from TLE), got %d (calculated)", tle.CheckSum1, calculatedCheckSum1)
	}

	calculatedCheckSum2, err := calculateChecksum(line2)
	if err != nil {
		return nil, fmt.Errorf("checksum calculation failed for line 2: %w", err)
	}
	if calculatedCheckSum2 != tle.CheckSum2 {
		return nil, fmt.Errorf("checksum mismatch in line 2: expected %d (from TLE), got %d (calculated)", tle.CheckSum2, calculatedCheckSum2)
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

	yearVal, err := strconv.Atoi(strings.TrimSpace(line[18:20]))
	if err != nil {
		return fmt.Errorf("invalid epoch year: %w", err)
	}
	if yearVal < 57 { // YY < 57 represents 20YY (e.g., 24 -> 2024)
		tle.EpochYear = 2000 + yearVal
	} else { // YY >= 57 represents 19YY (e.g., 98 -> 1998)
		tle.EpochYear = 1900 + yearVal
	}

	tle.EpochDay, err = strconv.ParseFloat(strings.TrimSpace(line[20:32]), 64)
	if err != nil {
		return fmt.Errorf("invalid epoch day: %w", err)
	}

	// Mean Motion Dot (First derivative of Mean Motion / 2)
	// Field 34-43: ".XXXXX" or " S.XXXXXFEE" (S=sign, X=digit, F=fraction, E=exponent)
	// Example: " .00007749"
	mmdStr := line[33:43]
	if strings.HasPrefix(mmdStr, ".") { // Implicit leading zero before decimal
		mmdStr = "0" + mmdStr
	} else if strings.HasPrefix(mmdStr, " ") && (mmdStr[1] == '.' || mmdStr[1] == '-') { // " .xxx" or " -.xxx"
		// If " .xxx" then " 0.xxx", if " -.xxx" then "-0.xxx"
		if mmdStr[1] == '.' {
			mmdStr = strings.Replace(mmdStr, " .", " 0.", 1)
		} else if mmdStr[1] == '-' && mmdStr[2] == '.' {
			mmdStr = strings.Replace(mmdStr, "-.", "-0.", 1)
		}
	}
	tle.MeanMotionDot, err = strconv.ParseFloat(strings.TrimSpace(mmdStr), 64)
	if err != nil {
		return fmt.Errorf("invalid mean motion dot ('%s'): %w", line[33:43], err)
	}

	// Mean Motion Dot Dot (Second derivative of Mean Motion / 6) (decimal point assumed)
	// Field 45-52: " SXXXXX±E" (S=sign, X=digit, E=exponent digit) e.g., " 00000+0"
	mmd2Str := line[44:50]
	mmd2ExpStr := line[50:52]

	mmd2Val, err := strconv.ParseFloat(strings.TrimSpace(mmd2Str), 64)
	if err != nil {
		return fmt.Errorf("invalid mean motion dot 2 mantissa ('%s'): %w", mmd2Str, err)
	}
	mmd2Exp, err := strconv.ParseInt(strings.TrimSpace(mmd2ExpStr), 10, 64)
	if err != nil {
		return fmt.Errorf("invalid mean motion dot 2 exponent ('%s'): %w", mmd2ExpStr, err)
	}
	tle.MeanMotionDot2 = mmd2Val * 1e-5 * math.Pow(10, float64(mmd2Exp))

	// B* Drag Term (decimal point assumed)
	// Field 54-61: " SXXXXX±E" (S=sign, X=digit, E=exponent digit) e.g., " 14567-3"
	bstarMantissaStr := line[53:59]
	bstarExponentStr := line[59:61]

	bstarMantissa, err := strconv.ParseFloat(strings.TrimSpace(bstarMantissaStr), 64)
	if err != nil {
		return fmt.Errorf("invalid B* mantissa ('%s'): %w", bstarMantissaStr, err)
	}
	bstarExponent, err := strconv.ParseInt(strings.TrimSpace(bstarExponentStr), 10, 64)
	if err != nil {
		return fmt.Errorf("invalid B* exponent ('%s'): %w", bstarExponentStr, err)
	}
	tle.Bstar = bstarMantissa * 1e-5 * math.Pow(10, float64(bstarExponent))

	// Element Set Type (usually 0) and Element Number
	// Field 63 is Ephemeris Type, field 65-68 is Element Number.
	// We can ignore ephemeris type for now or parse if needed.
	tle.ElementNumber, err = strconv.Atoi(strings.TrimSpace(line[64:68]))
	if err != nil {
		return fmt.Errorf("invalid element number: %w", err)
	}

	checksum, err := strconv.Atoi(strings.TrimSpace(line[68:69]))
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
		return fmt.Errorf("satellite numbers do not match between lines (%d vs %d)", tle.SatelliteNumber, satNum)
	}

	tle.Inclination, err = strconv.ParseFloat(strings.TrimSpace(line[8:16]), 64)
	if err != nil {
		return fmt.Errorf("invalid inclination: %w", err)
	}

	tle.RightAscension, err = strconv.ParseFloat(strings.TrimSpace(line[17:25]), 64)
	if err != nil {
		return fmt.Errorf("invalid right ascension: %w", err)
	}

	// Eccentricity (decimal point assumed: XXXXXXX -> 0.XXXXXXX)
	eccStr := line[26:33]
	ecc, err := strconv.ParseFloat("0."+strings.TrimSpace(eccStr), 64)
	if err != nil {
		return fmt.Errorf("invalid eccentricity ('%s'): %w", eccStr, err)
	}
	tle.Eccentricity = ecc

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

	checksum, err := strconv.Atoi(strings.TrimSpace(line[68:69]))
	if err != nil {
		return fmt.Errorf("invalid checksum: %w", err)
	}
	tle.CheckSum2 = checksum

	return nil
}

// RecoveredSemiMajorAxis returns the semi-major axis (in Earth radii)
// calculated from the mean motion using Kepler's Third Law
func (tle *TLE) RecoveredSemiMajorAxis() float64 {
	// Mean motion in radians per minute
	meanMotionRadPerMin := tle.MeanMotion * twoPi / minutesPerDay

	// Calculate semi-major axis using the relation between mean motion and semi-major axis
	// from Kepler's Third Law: n² * a³ = μ
	// where:
	// n is mean motion in rad/min
	// a is semi-major axis in Earth radii
	// μ is GM * (1/xke)^2 where GM is the gravitational parameter
	// xke = sqrt(GM) / (EarthRadius^(3/2)) * minutesPerTU (approximately)
	// a = ( (GM / n^2)^(1/3) ) / EarthRadius
	// a_ER = ( (xke / n_rad_min)^2 )^(1/3)
	return math.Pow(xke/meanMotionRadPerMin, 2.0/3.0)
}

// IsGeostationary checks if the satellite's orbital elements from the TLE
// suggest it is likely a geostationary satellite.
// This is based on mean elements and does not guarantee perfect station-keeping.
func (tle *TLE) IsGeostationary() bool {
	// Thresholds for geostationary characteristics
	const meanMotionGeoMin = 0.99 // revs/day (slightly less than 1 sidereal rotation)
	const meanMotionGeoMax = 1.01 // revs/day (slightly more than 1 sidereal rotation)
	// Ideal GEO mean motion is approx 1.0027 rev/day relative to fixed stars.
	// If TLE mean motion is mean solar day relative, then it's closer to 1.0.
	// Let's use a slightly wider band around 1.0027.
	// Period = 1436.068 minutes -> Mean motion = 1440 / 1436.068 = 1.0027379 revs/day
	const idealGeoMeanMotion = minutesPerDay / (minutesPerDay / 1.0027379093509) // approx 1.0027379 revs/day
	const meanMotionTolerance = 0.05                                             // Allow +/- this much from ideal mean motion

	const maxInclinationDeg = 5.0 // degrees (some GEO sats have slight inclination)
	const maxEccentricity = 0.05  // (some GEO sats have slight eccentricity)

	// 1. Check Mean Motion (period)
	// tle.MeanMotion is in revolutions per day.
	if tle.MeanMotion < (idealGeoMeanMotion-meanMotionTolerance) || tle.MeanMotion > (idealGeoMeanMotion+meanMotionTolerance) {
		return false
	}

	// 2. Check Inclination
	if tle.Inclination > maxInclinationDeg { // Inclination from TLE is already in degrees
		return false
	}

	// 3. Check Eccentricity
	if tle.Eccentricity > maxEccentricity { // Eccentricity from TLE (0.XXXXXXX format)
		return false
	}

	// If all checks pass, it's likely geostationary or geosynchronous.
	// True geostationary implies inclination is very close to zero.
	// Geosynchronous can have inclination but still a 24-hour (sidereal) period.
	// The term "geostationary" usually implies near-zero inclination.
	return true
}

// SGP4Constants holds values used in SGP4 orbital calculations
type SGP4Constants struct {
	Sinio  float64 // Sine of inclination
	Cosio  float64 // Cosine of inclination
	X3thm1 float64 // 3*cos²(i) - 1
	X1mth2 float64 // 1 - cos²(i)
	X7thm1 float64 // 7*cos²(i) - 1
	Xlcof  float64 // Long period periodic coefficient
	Aycof  float64 // Another long period periodic coefficient
}

func NewSGP4Constants(xinc float64) *SGP4Constants {
	s := &SGP4Constants{}
	s.Sinio, s.Cosio, s.X3thm1, s.X1mth2, s.X7thm1, s.Xlcof, s.Aycof = RecomputeConstants(xinc)
	return s
}

// RecomputeConstants calculates values needed for the SGP4 propagator
// `sinio` - Sine of inclination
// `cosio` - Cosine of inclination
// `x3thm1` - Three times theta squared minus one
// `x1mth2` - One minus theta squared
// `x7thm1` - Seven times theta squared minus one
// `xlcof` - Long period periodic coefficient
// `aycof` - Another long period periodic coefficient
func RecomputeConstants(xinc float64) (sinio, cosio, x3thm1, x1mth2, x7thm1, xlcof, aycof float64) {
	sinio = math.Sin(xinc)
	cosio = math.Cos(xinc)

	theta2 := cosio * cosio

	// Calculate periodic terms
	x3thm1 = 3.0*theta2 - 1.0 // Three times theta squared minus one
	x1mth2 = 1.0 - theta2     // One minus theta squared
	x7thm1 = 7.0*theta2 - 1.0 // Seven times theta squared minus one

	// Calculate xlcof with protection against division by zero
	if math.Abs(cosio+1.0) > 1.5e-12 {
		xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0*cosio) / (1.0 + cosio)
	} else {
		xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0*cosio) / 1.5e-12
	}

	// Calculate aycof
	aycof = 0.25 * a3ovk2 * sinio

	return
}

// calculateChecksum calculates the modulo-10 checksum for a TLE line.
// It sums all numerical digits, with '-' counting as 1. Other characters are ignored.
// The checksum is calculated for the first 68 characters of the line.
func calculateChecksum(line string) (int, error) {
	if len(line) != 69 {
		return 0, fmt.Errorf("line must be 69 characters long, got %d", len(line))
	}
	sum := 0
	for i := range 68 { // Iterate over characters 0-67
		char := line[i]
		if char >= '0' && char <= '9' {
			sum += int(char - '0')
		} else if char == '-' {
			sum += 1
		}
		// All other characters (letters, space, '.', '+') are ignored.
	}
	return sum % 10, nil
}
