package sgp4

import (
	"fmt"
	"math"
	"strconv"
	"strings"
	"time"
	"unicode/utf8"
)

// The catalog-number field of a TLE holds five characters and so cannot carry a
// six-digit number as digits. From 100000 the leading digit is replaced by a
// letter -- the Alpha-5 scheme Space-Track uses -- with A=10, B=11, and so on to
// Z=33. I and O are skipped so they cannot be mistaken for 1 and 0, which puts
// the ceiling at 339999 (Z9999). Objects above that have no TLE representation
// and must be carried as OMM/GP data instead.
const (
	alpha5Letters     = "ABCDEFGHJKLMNPQRSTUVWXYZ"
	alpha5FirstLetter = 10
	// The last prefix letter is Z, so the ceiling is Z9999: 33*10000 + 9999.
	alpha5LastLetter = alpha5FirstLetter + len(alpha5Letters) - 1
	alpha5Max        = alpha5LastLetter*10000 + 9999
)

// ParseAlpha5 decodes a five-character TLE catalog-number field. Below 100000
// the field is plain digits with leading zeros ("00964" is 964); from 100000 to
// 339999 the leading digit is a letter, so "A0000" is 100000 and "Z9999" is
// 339999. Lowercase letters, the letters I and O, and non-digit characters after
// the first are rejected, matching Space-Track's definition.
func ParseAlpha5(field string) (int, error) {
	if len(field) != 5 {
		return 0, fmt.Errorf("sgp4: invalid catalog field %q: want 5 characters, got %d", field, len(field))
	}
	if allDigits(field) {
		// Five ASCII digits always fit an int.
		return strconv.Atoi(field)
	}

	prefix := strings.IndexByte(alpha5Letters, field[0])
	if prefix < 0 {
		return 0, fmt.Errorf("sgp4: invalid catalog field %q: want digits or an uppercase Alpha-5 prefix (A-Z without I and O)", field)
	}
	if !allDigits(field[1:]) {
		return 0, fmt.Errorf("sgp4: invalid catalog field %q: the four characters after the Alpha-5 prefix must be digits", field)
	}
	tail, _ := strconv.Atoi(field[1:])
	return (alpha5FirstLetter+prefix)*10000 + tail, nil
}

// FormatAlpha5 encodes a catalog number as a five-character TLE field: five
// digits below 100000, Alpha-5 from 100000 to 339999. Larger or negative numbers
// have no TLE representation and return an error.
func FormatAlpha5(n int) (string, error) {
	if n < 0 || n > alpha5Max {
		return "", fmt.Errorf("sgp4: catalog number %d cannot be written as a TLE catalog field: the field holds 0 to %d", n, alpha5Max)
	}
	if n < alpha5FirstLetter*10000 {
		return fmt.Sprintf("%05d", n), nil
	}
	return fmt.Sprintf("%c%04d", alpha5Letters[n/10000-alpha5FirstLetter], n%10000), nil
}

// parseCatalogNumber decodes the catalog field of a TLE line. Older element sets
// sometimes space-pad short numbers instead of zero-padding them, so surrounding
// spaces are ignored; anything else goes through ParseAlpha5 unchanged.
func parseCatalogNumber(field string) (int, error) {
	field = strings.TrimSpace(field)
	if field == "" {
		return 0, fmt.Errorf("sgp4: empty catalog number field")
	}
	if len(field) <= 5 && allDigits(field) {
		return strconv.Atoi(field)
	}
	return ParseAlpha5(field)
}

func allDigits(s string) bool {
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return false
		}
	}
	return len(s) > 0
}

// TLELines renders the element set as TLE lines 1 and 2, computing the
// checksums. A value the fixed-column layout cannot hold -- a catalog number
// above 339999, an epoch outside 1957-2056, an angle out of range -- returns an
// error rather than being silently truncated.
//
// Values are rounded to the precision of their field, so a TLE parsed by this
// package renders back byte for byte when its fields are canonical.
func (tle *TLE) TLELines() ([2]string, error) {
	var lines [2]string

	sat, err := FormatAlpha5(tle.SatelliteNumber)
	if err != nil {
		return lines, err
	}

	classification := tle.Classification
	if classification == 0 {
		classification = 'U'
	}
	if classification < 0x20 || classification > 0x7e {
		return lines, fmt.Errorf("sgp4: classification %q is not a printable character", classification)
	}

	intl := strings.TrimSpace(tle.International)
	if len(intl) > 8 {
		return lines, fmt.Errorf("sgp4: international designator %q is longer than the 8 characters a TLE holds", tle.International)
	}

	if tle.EpochYear < 1957 || tle.EpochYear > 2056 {
		return lines, fmt.Errorf("sgp4: epoch year %d cannot be written in a TLE: the two-digit year covers 1957-2056", tle.EpochYear)
	}
	// The day-of-year field is 1-based and the year may have 365 or 366 days.
	lastDay := float64(time.Date(tle.EpochYear, time.December, 31, 0, 0, 0, 0, time.UTC).YearDay())
	if tle.EpochDay < 1 || tle.EpochDay >= lastDay+1 {
		return lines, fmt.Errorf("sgp4: epoch day %.8f is outside the %d-day year %d", tle.EpochDay, int(lastDay), tle.EpochYear)
	}

	ndot, err := tleNdotField(tle.MeanMotionDot)
	if err != nil {
		return lines, err
	}
	ddot, err := tleExpField(tle.MeanMotionDot2, "MEAN_MOTION_DDOT")
	if err != nil {
		return lines, err
	}
	bstar, err := tleExpField(tle.Bstar, "BSTAR")
	if err != nil {
		return lines, err
	}

	if tle.EphemerisType < 0 || tle.EphemerisType > 9 {
		return lines, fmt.Errorf("sgp4: ephemeris type %d does not fit the single-column field", tle.EphemerisType)
	}
	if tle.ElementNumber < 0 || tle.ElementNumber > 9999 {
		return lines, fmt.Errorf("sgp4: element number %d is outside the TLE range 0-9999", tle.ElementNumber)
	}

	if tle.Inclination < 0 || tle.Inclination > 180 {
		return lines, fmt.Errorf("sgp4: inclination %.4f is outside the TLE range [0, 180]", tle.Inclination)
	}
	if tle.RightAscension < 0 || tle.RightAscension > 360 {
		return lines, fmt.Errorf("sgp4: right ascension %.4f is outside the TLE range [0, 360]", tle.RightAscension)
	}
	if tle.ArgOfPerigee < 0 || tle.ArgOfPerigee > 360 {
		return lines, fmt.Errorf("sgp4: argument of perigee %.4f is outside the TLE range [0, 360]", tle.ArgOfPerigee)
	}
	if tle.MeanAnomaly < 0 || tle.MeanAnomaly > 360 {
		return lines, fmt.Errorf("sgp4: mean anomaly %.4f is outside the TLE range [0, 360]", tle.MeanAnomaly)
	}
	ecc, err := tleEccField(tle.Eccentricity)
	if err != nil {
		return lines, err
	}
	if tle.MeanMotion <= 0 || tle.MeanMotion >= 100 {
		return lines, fmt.Errorf("sgp4: mean motion %.8f is outside the TLE range (0, 100)", tle.MeanMotion)
	}
	if tle.RevolutionNumber < 0 || tle.RevolutionNumber > 99999 {
		return lines, fmt.Errorf("sgp4: revolution number %d is outside the TLE range 0-99999", tle.RevolutionNumber)
	}

	l1 := fmt.Sprintf("1 %s%c %-8s %02d%012.8f %s %s %s %d %4d",
		sat, classification, intl, tle.EpochYear%100, tle.EpochDay,
		ndot, ddot, bstar, tle.EphemerisType, tle.ElementNumber)
	l2 := fmt.Sprintf("2 %s %8.4f %8.4f %s %8.4f %8.4f %11.8f%5d",
		sat, tle.Inclination, tle.RightAscension, ecc,
		tle.ArgOfPerigee, tle.MeanAnomaly, tle.MeanMotion, tle.RevolutionNumber)

	ck1, err := calculateChecksum(l1 + "0")
	if err != nil {
		return lines, fmt.Errorf("sgp4: computing line 1 checksum: %w", err)
	}
	ck2, err := calculateChecksum(l2 + "0")
	if err != nil {
		return lines, fmt.Errorf("sgp4: computing line 2 checksum: %w", err)
	}

	lines[0] = l1 + strconv.Itoa(ck1)
	lines[1] = l2 + strconv.Itoa(ck2)
	return lines, nil
}

// TLEString renders the element set as three-line TLE text: the object name
// (padded or cut to the name line's 24 characters) followed by lines 1 and 2.
// When the name is empty only the two element lines are returned.
func (tle *TLE) TLEString() (string, error) {
	lines, err := tle.TLELines()
	if err != nil {
		return "", err
	}
	if strings.TrimSpace(tle.Name) == "" {
		return lines[0] + "\n" + lines[1], nil
	}
	name := tle.Name
	if len(name) > 24 {
		name = name[:24]
		// Do not cut a multi-byte character in half.
		for !utf8.ValidString(name) {
			name = name[:len(name)-1]
		}
	}
	return name + strings.Repeat(" ", 24-len(name)) + "\n" + lines[0] + "\n" + lines[1], nil
}

// tleNdotField renders MEAN_MOTION_DOT as the ten-column implied-decimal field:
// a sign (space or '-'), a '.', and eight digits.
func tleNdotField(v float64) (string, error) {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return "", fmt.Errorf("sgp4: MEAN_MOTION_DOT %v cannot be written in a TLE", v)
	}
	if math.Abs(v) >= 1 {
		return "", fmt.Errorf("sgp4: MEAN_MOTION_DOT %g is too large for the TLE field", v)
	}
	units := math.Round(math.Abs(v) * 1e8)
	if units >= 1e8 {
		return "", fmt.Errorf("sgp4: MEAN_MOTION_DOT %g rounds to 1.0, which the TLE field cannot hold", v)
	}
	sign := byte(' ')
	if v < 0 {
		sign = '-'
	}
	return fmt.Sprintf("%c.%08d", sign, int64(units)), nil
}

// tleExpField renders an implied-decimal field with an exponent (BSTAR,
// MEAN_MOTION_DDOT) as a sign, five mantissa digits, an exponent sign and one
// exponent digit, for example " 14567-3".
func tleExpField(v float64, name string) (string, error) {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return "", fmt.Errorf("sgp4: %s %v cannot be written in a TLE", name, v)
	}
	sign := byte(' ')
	abs := v
	if v < 0 {
		sign = '-'
		abs = -v
	}
	if abs == 0 {
		return " 00000+0", nil
	}

	// Normalise to 0.ddddd * 10^exponent; the field holds -9 to 9.
	exponent := int(math.Floor(math.Log10(abs))) + 1
	mantissa := math.Round(abs / math.Pow(10, float64(exponent)) * 1e5)
	if mantissa >= 1e5 {
		mantissa = math.Round(mantissa / 10)
		exponent++
	}
	if mantissa < 1e4 || exponent < -9 || exponent > 9 {
		return "", fmt.Errorf("sgp4: %s %g is outside the range the TLE exponent field holds", name, v)
	}

	expSign := byte('+')
	if exponent < 0 {
		expSign = '-'
		exponent = -exponent
	}
	return fmt.Sprintf("%c%05d%c%d", sign, int64(mantissa), expSign, exponent), nil
}

// tleEccField renders eccentricity as the seven implied-decimal digits.
func tleEccField(e float64) (string, error) {
	if math.IsNaN(e) || e < 0 || e >= 1 {
		return "", fmt.Errorf("sgp4: eccentricity %v is outside the TLE range [0, 1)", e)
	}
	digits := math.Round(e * 1e7)
	if digits >= 1e7 {
		return "", fmt.Errorf("sgp4: eccentricity %g rounds to 1.0000000, which the TLE field cannot hold", e)
	}
	return fmt.Sprintf("%07d", int64(digits)), nil
}
