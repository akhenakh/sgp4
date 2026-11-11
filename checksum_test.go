package sgp4

import "testing"

func TestChecksumKnownLines(t *testing.T) {
	// Real element sets: their declared checksum must match the computed one.
	lines := []string{
		"1 25544U 98067A   20101.90690972 -.00000449  00000-0  00000+0 0  9993",
		"2 25544  51.6446 321.4198 0003848 108.5166  84.0719 15.48680394221581",
		"1 00900U 64063C   26314.18770255  .00000675  00000+0  67559-3 0  9991",
		"2 00900  90.2236  75.1584 0026167 206.8891 227.1268 13.76759504 91510",
		// Alpha-5 catalog field (letters count zero).
		"1 A0000U 26067CY  26195.90649229  .00004770  00000+0  22159-3 0  9994",
		"2 A0000  97.4593 154.0970 0005590 270.5113  89.5482 15.20467281 15911",
	}
	for _, line := range lines {
		computed, err := calculateChecksum(line)
		if err != nil {
			t.Fatalf("calculateChecksum(%q): %v", line, err)
		}
		declared := int(line[68] - '0')
		if computed != declared {
			t.Errorf("line %q: computed checksum %d, declared %d", line, computed, declared)
		}
	}
}

func TestChecksumRules(t *testing.T) {
	// Digits count as their value, '-' counts 1, and every other character
	// (letters, '+', '.', space) counts 0.
	line := "1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994"
	computed, err := calculateChecksum(line)
	if err != nil {
		t.Fatalf("calculateChecksum: %v", err)
	}

	// Recompute independently with the spec rules.
	sum := 0
	for i := 0; i < 68; i++ {
		switch c := line[i]; {
		case c >= '0' && c <= '9':
			sum += int(c - '0')
		case c == '-':
			sum++
		}
	}
	if want := sum % 10; computed != want {
		t.Errorf("computed checksum %d, want %d", computed, want)
	}

	// A letter in the catalog field must not contribute.
	alpha, err := calculateChecksum("1 A0000U 26067CY  26195.90649229  .00004770  00000+0  22159-3 0  9994")
	if err != nil {
		t.Fatalf("calculateChecksum: %v", err)
	}
	numeric, err := calculateChecksum("1 00000U 26067CY  26195.90649229  .00004770  00000+0  22159-3 0  9994")
	if err != nil {
		t.Fatalf("calculateChecksum: %v", err)
	}
	if alpha != numeric {
		t.Errorf("letter changed the checksum: %d vs %d", alpha, numeric)
	}
}

func TestChecksumRejectsWrongChecksum(t *testing.T) {
	line := "1 25544U 98067A   20101.90690972 -.00000449  00000-0  00000+0 0  9993"
	bad := line[:68] + string('0'+(line[68]-'0'+1)%10)
	if _, err := ParseTLE(bad + "\n" + "2 25544  51.6446 321.4198 0003848 108.5166  84.0719 15.48680394221581"); err == nil {
		t.Fatal("ParseTLE accepted a line with a wrong checksum")
	}
}

func TestChecksumRejectsWrongLength(t *testing.T) {
	if _, err := calculateChecksum("1 25544U 98067A"); err == nil {
		t.Fatal("calculateChecksum accepted a short line")
	}
}
