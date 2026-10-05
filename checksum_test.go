package sgp4

import (
	"testing"
)

func TestChecksum(t *testing.T) {
	// ISS TLE from sgp4_test.go

	lines := []string{
		"1 25544U 98067A   20101.90690972 -.00000449  00000-0  00000+0 0  9993",
		"2 25544  51.6446 321.4198 0003848 108.5166  84.0719 15.48680394221581",
	}

	for i, line := range lines {
		computedChecksum, err := calculateChecksum(line)
		if err != nil {
			t.Errorf("Error calculating checksum for line %d: %v", i+1, err)
			continue
		}

		actualChecksum := int(line[68] - '0')
		t.Logf("Line %d: Computed checksum = %d, Actual checksum = %d",
			i+1, computedChecksum, actualChecksum)

		// Print character by character examination of first line for debugging
		if i == 0 {
			t.Log("Character by character examination of line 1:")
			sum := 0
			for j := 0; j < len(line)-1; j++ {
				char := line[j]
				charValue := 0

				switch {
				case char >= '0' && char <= '9':
					charValue = int(char - '0')
					sum += charValue
				case char == '-' || char == '+':
					charValue = 1
					sum += charValue
				case char >= 'A' && char <= 'Z':
					charValue = int(char - 'A' + 1)
					sum += charValue
				case char == ' ' || char == '.':
					// Ignore spaces and decimal points
					charValue = 0
				}

				t.Logf("Pos %2d: '%c' -> Value: %2d, Running sum: %3d",
					j, char, charValue, sum)
			}
			t.Logf("Final sum: %d, Checksum (mod 10): %d", sum, sum%10)
		}
	}
}

// A simplified version of the main checksum calculator for comparison
func manualChecksum(line string) int {
	sum := 0
	for i := 0; i < len(line)-1; i++ {
		char := line[i]
		switch {
		case char >= '0' && char <= '9':
			sum += int(char - '0')
		case char == '-' || char == '+':
			sum += 1
		case char >= 'A' && char <= 'Z':
			sum += int(char - 'A' + 1)
			// Spaces and decimal points are ignored
		}
	}
	return sum % 10
}
