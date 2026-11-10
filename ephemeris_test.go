package sgp4

import (
	"errors"
	"strconv"
	"testing"
)

func TestUnsupportedEphemerisType(t *testing.T) {
	// The issue's second example: column 63 holds 4, which marks an SGP4-XP set.
	line1 := "1 20453U 90008B   26263.56914166  .00350177  49935-4  75989-3 4  9990"
	line2 := "2 20453  35.5934 307.3932 0022512 310.3690  49.5094 15.96788691956793"

	tle, err := ParseTLE("DELTA 2 R/B(1)\n" + line1 + "\n" + line2)
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	if tle.EphemerisType != 4 {
		t.Fatalf("EphemerisType = %d, want 4", tle.EphemerisType)
	}
	if _, err := tle.FindPosition(0); err == nil {
		t.Fatal("FindPosition accepted an SGP4-XP element set")
	} else {
		var typeErr *UnsupportedEphemerisTypeError
		if !errors.As(err, &typeErr) {
			t.Fatalf("FindPosition error = %v, want UnsupportedEphemerisTypeError", err)
		}
		if typeErr.Type != 4 {
			t.Errorf("error reports type %d, want 4", typeErr.Type)
		}
	}

	// The same set declared as SGP4 (type 0) propagates, so the refusal is about
	// the declared theory and not the element values.
	line1 = line1[:62] + "0" + line1[63:68]
	ck, _ := calculateChecksum(line1 + "0")
	line1 += strconv.Itoa(ck)
	tle, err = ParseTLE(line1 + "\n" + line2)
	if err != nil {
		t.Fatalf("ParseTLE of type-0 set: %v", err)
	}
	if _, err := tle.FindPosition(0); err != nil {
		t.Errorf("FindPosition of type-0 set: %v", err)
	}
}

func TestEphemerisTypeOMMRoundTrip(t *testing.T) {
	tle, err := ParseTLE("1 20453U 90008B   26263.56914166  .00350177  49935-4  75989-3 4  9990\n" +
		"2 20453  35.5934 307.3932 0022512 310.3690  49.5094 15.96788691956793")
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	omm := tle.ToOMM()
	if omm.EphemerisType != 4 {
		t.Errorf("ToOMM EphemerisType = %d, want 4", omm.EphemerisType)
	}
	back, err := omm.ToTLE()
	if err != nil {
		t.Fatalf("ToTLE: %v", err)
	}
	if back.EphemerisType != 4 {
		t.Errorf("ToTLE EphemerisType = %d, want 4", back.EphemerisType)
	}
	if _, err := back.Initialize(); err == nil {
		t.Error("Initialize accepted an SGP4-XP set converted from OMM")
	}
}
