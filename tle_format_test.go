package sgp4

import (
	"strconv"
	"strings"
	"testing"
)

// The first 68 characters of a known-good ISS element set; the checksum is
// dropped so tests can rewrite the catalog field and recompute it.
const (
	line1Prefix = "1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  999"
	line2Prefix = "2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.4958795751053"
)

// tleWithCatalog replaces the catalog field of both lines and recomputes the
// checksums, so the result is a valid TLE carrying the given field.
func tleWithCatalog(line1, line2, field string) (string, string) {
	l1 := line1[:2] + field + line1[7:]
	l2 := line2[:2] + field + line2[7:]
	c1, _ := calculateChecksum(l1 + "0")
	c2, _ := calculateChecksum(l2 + "0")
	return l1 + strconv.Itoa(c1), l2 + strconv.Itoa(c2)
}

func TestParseAlpha5(t *testing.T) {
	valid := []struct {
		field string
		want  int
	}{
		{"00000", 0},
		{"00005", 5},
		{"00964", 964},
		{"25544", 25544},
		{"99999", 99999},
		{"A0000", 100000},
		{"A9999", 109999},
		{"B0000", 110000},
		{"E8493", 148493},
		{"J2931", 182931},
		{"P4018", 234018},
		{"W1928", 301928},
		{"Z0000", 330000},
		{"Z9999", 339999},
	}
	for _, tt := range valid {
		got, err := ParseAlpha5(tt.field)
		if err != nil {
			t.Errorf("ParseAlpha5(%q): %v", tt.field, err)
			continue
		}
		if got != tt.want {
			t.Errorf("ParseAlpha5(%q) = %d, want %d", tt.field, got, tt.want)
		}
	}

	invalid := []string{
		"", "A000", "A00000", " 964 ", "a0000", "I0000", "O0000",
		"8E493", "E84A3", "E84 3", "A-001", "+1234", "A.000",
	}
	for _, field := range invalid {
		if got, err := ParseAlpha5(field); err == nil {
			t.Errorf("ParseAlpha5(%q) = %d, want an error", field, got)
		}
	}
}

func TestFormatAlpha5(t *testing.T) {
	// The documented examples, including the top of every letter band.
	examples := map[int]string{
		0:      "00000",
		964:    "00964",
		99999:  "99999",
		100000: "A0000",
		148493: "E8493",
		182931: "J2931",
		234018: "P4018",
		301928: "W1928",
		339999: "Z9999",
	}
	for n, want := range examples {
		got, err := FormatAlpha5(n)
		if err != nil || got != want {
			t.Errorf("FormatAlpha5(%d) = %q, %v; want %q", n, got, err, want)
		}
	}

	for n := 0; n <= alpha5Max; n++ {
		field, err := FormatAlpha5(n)
		if err != nil {
			t.Fatalf("FormatAlpha5(%d): %v", n, err)
		}
		if len(field) != 5 {
			t.Fatalf("FormatAlpha5(%d) = %q, want 5 characters", n, field)
		}
		got, err := ParseAlpha5(field)
		if err != nil {
			t.Fatalf("ParseAlpha5(%q): %v", field, err)
		}
		if got != n {
			t.Fatalf("ParseAlpha5(FormatAlpha5(%d)) = %d", n, got)
		}
	}

	for _, n := range []int{-1, 340000, 799501621} {
		if field, err := FormatAlpha5(n); err == nil {
			t.Errorf("FormatAlpha5(%d) = %q, want an error", n, field)
		}
	}
}

func TestParseTLEAlpha5(t *testing.T) {
	tests := []struct {
		field   string
		want    int
		wantErr bool
	}{
		{"A0000", 100000, false},
		{"E8493", 148493, false},
		{"Z9999", 339999, false},
		{"I0000", 0, true},
		{"O0000", 0, true},
		{"a0000", 0, true},
		{"8E493", 0, true},
		{"E84A3", 0, true},
	}
	for _, tt := range tests {
		l1, l2 := tleWithCatalog(line1Prefix, line2Prefix, tt.field)
		tle, err := ParseTLE(l1 + "\n" + l2)
		if tt.wantErr {
			if err == nil {
				t.Errorf("ParseTLE with catalog field %q: want an error, got number %d", tt.field, tle.SatelliteNumber)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseTLE with catalog field %q: %v", tt.field, err)
			continue
		}
		if tle.SatelliteNumber != tt.want {
			t.Errorf("catalog field %q parsed as %d, want %d", tt.field, tle.SatelliteNumber, tt.want)
		}
	}
}

func TestParseTLEAlpha5Mismatch(t *testing.T) {
	// Line 2 keeps digits while line 1 is Alpha-5: same number, different
	// encoding, and the TLE format requires the fields to match.
	l1, _ := tleWithCatalog(line1Prefix, line2Prefix, "A0000")
	_, l2 := tleWithCatalog(line1Prefix, line2Prefix, "A0001")
	if _, err := ParseTLE(l1 + "\n" + l2); err == nil {
		t.Fatal("ParseTLE accepted catalog numbers that differ between the lines")
	}
}

func TestAlpha5Propagation(t *testing.T) {
	numeric, err := ParseTLE(line1Prefix + "4\n" + line2Prefix + "3")
	if err != nil {
		t.Fatalf("parsing numeric TLE: %v", err)
	}
	l1, l2 := tleWithCatalog(line1Prefix, line2Prefix, "E8493")
	alpha, err := ParseTLE(l1 + "\n" + l2)
	if err != nil {
		t.Fatalf("parsing Alpha-5 TLE: %v", err)
	}

	for _, tsince := range []float64{0, 45, 1440} {
		want, err := numeric.FindPosition(tsince)
		if err != nil {
			t.Fatalf("numeric FindPosition(%v): %v", tsince, err)
		}
		got, err := alpha.FindPosition(tsince)
		if err != nil {
			t.Fatalf("Alpha-5 FindPosition(%v): %v", tsince, err)
		}
		if got.Position != want.Position || got.Velocity != want.Velocity {
			t.Errorf("propagation depends on the catalog field at t=%v: got %v/%v, want %v/%v",
				tsince, got.Position, got.Velocity, want.Position, want.Velocity)
		}
	}
}

func TestParseTLEAlpha5Ticket(t *testing.T) {
	// The element set from the issue: SARAMAGO carries catalog number 100000.
	tle, err := ParseTLE("SARAMAGO\n" +
		"1 A0000U 26067CY  26195.90649229  .00004770  00000+0  22159-3 0  9994\n" +
		"2 A0000  97.4593 154.0970 0005590 270.5113  89.5482 15.20467281 15911")
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	if tle.SatelliteNumber != 100000 {
		t.Errorf("SatelliteNumber = %d, want 100000", tle.SatelliteNumber)
	}
	if _, err := tle.FindPosition(0); err != nil {
		t.Errorf("FindPosition: %v", err)
	}
	if tle.Name != "SARAMAGO" {
		t.Errorf("Name = %q, want SARAMAGO", tle.Name)
	}
}

func TestTLELinesRoundTrip(t *testing.T) {
	tles := []string{
		// ISS: negative mean-motion derivative and a B* with an exponent.
		line1Prefix + "4\n" + line2Prefix + "3",
		// Empty BALLAST: a GEO element set, and a negative B* below.
		"1 02866U 67066E   26313.83598808 -.00000063  00000+0  00000+0 0  9999\n" +
			"2 02866   2.9580  92.2224 0054005 216.5249 167.3784  1.09425390132636",
		"1 16908U 86061A   26314.28578718 -.00000112  00000+0 -10232-3 0  9993\n" +
			"2 16908  50.0104  58.1417 0011348 191.5678 201.4379 12.44517677498580",
		// The issue's Alpha-5 set.
		"1 A0000U 26067CY  26195.90649229  .00004770  00000+0  22159-3 0  9994\n" +
			"2 A0000  97.4593 154.0970 0005590 270.5113  89.5482 15.20467281 15911",
	}

	for _, input := range tles {
		parsed, err := ParseTLE(input)
		if err != nil {
			t.Errorf("parsing %q: %v", input, err)
			continue
		}
		want := strings.Split(input, "\n")
		lines, err := parsed.TLELines()
		if err != nil {
			t.Errorf("TLELines for %q: %v", want[0], err)
			continue
		}
		if lines[0] != want[0] || lines[1] != want[1] {
			t.Errorf("round trip of %q:\n got %q\n     %q\nwant %q\n     %q",
				want[0], lines[0], lines[1], want[0], want[1])
		}
	}
}

func TestTLELinesAlpha5(t *testing.T) {
	// Rendering a set whose number is carried by a letter other than A: the
	// writer must accept the whole Alpha-5 range, not just the first band.
	for _, field := range []string{"E8493", "J2931", "P4018", "W1928", "Z9999"} {
		l1, l2 := tleWithCatalog(line1Prefix, line2Prefix, field)
		tle, err := ParseTLE(l1 + "\n" + l2)
		if err != nil {
			t.Errorf("parsing catalog field %q: %v", field, err)
			continue
		}
		lines, err := tle.TLELines()
		if err != nil {
			t.Errorf("rendering catalog field %q (number %d): %v", field, tle.SatelliteNumber, err)
			continue
		}
		if lines[0] != l1 || lines[1] != l2 {
			t.Errorf("round trip of catalog field %q:\n got %q\n     %q\nwant %q\n     %q",
				field, lines[0], lines[1], l1, l2)
		}
	}
}

func TestTLEString(t *testing.T) {
	tle, err := ParseTLE("1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994\n" +
		line2Prefix + "3")
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	tle.Name = "ISS (ZARYA)"
	out, err := tle.TLEString()
	if err != nil {
		t.Fatalf("TLEString: %v", err)
	}
	lines := strings.Split(out, "\n")
	if len(lines) != 3 {
		t.Fatalf("TLEString returned %d lines, want 3", len(lines))
	}
	if len(lines[0]) != 24 || !strings.HasPrefix(lines[0], "ISS (ZARYA)") {
		t.Errorf("name line = %q, want 24 columns starting with the name", lines[0])
	}
	if _, err := ParseTLE(out); err != nil {
		t.Errorf("TLEString output does not parse back: %v", err)
	}
}

func TestTLELinesRefusals(t *testing.T) {
	base, err := ParseTLE(line1Prefix + "4\n" + line2Prefix + "3")
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}

	tests := []struct {
		name   string
		mutate func(*TLE)
	}{
		{"catalog number above the Alpha-5 ceiling", func(t *TLE) { t.SatelliteNumber = 340000 }},
		{"nine-digit catalog number", func(t *TLE) { t.SatelliteNumber = 799501621 }},
		{"negative catalog number", func(t *TLE) { t.SatelliteNumber = -1 }},
		{"epoch year after 2056", func(t *TLE) { t.EpochYear = 2057 }},
		{"epoch year before 1957", func(t *TLE) { t.EpochYear = 1956 }},
		{"epoch day 367", func(t *TLE) { t.EpochDay = 367 }},
		{"epoch day 366 in a non-leap year", func(t *TLE) { t.EpochYear = 2025; t.EpochDay = 366 }},
		{"eccentricity of 1", func(t *TLE) { t.Eccentricity = 1 }},
		{"inclination above 180", func(t *TLE) { t.Inclination = 181 }},
		{"zero mean motion", func(t *TLE) { t.MeanMotion = 0 }},
		{"mean motion of 100", func(t *TLE) { t.MeanMotion = 100 }},
		{"element number too large", func(t *TLE) { t.ElementNumber = 10000 }},
		{"revolution number too large", func(t *TLE) { t.RevolutionNumber = 100000 }},
		{"mean motion derivative too large", func(t *TLE) { t.MeanMotionDot = 1 }},
		{"B* beyond the exponent field", func(t *TLE) { t.Bstar = 1e10 }},
		{"international designator too long", func(t *TLE) { t.International = "98067ABCD" }},
		{"ephemeris type not a single digit", func(t *TLE) { t.EphemerisType = 10 }},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			bad := *base
			tt.mutate(&bad)
			if lines, err := bad.TLELines(); err == nil {
				t.Errorf("TLELines accepted %s: %q / %q", tt.name, lines[0], lines[1])
			}
		})
	}
}
