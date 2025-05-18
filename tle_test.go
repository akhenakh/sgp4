package sgp4

import (
	"math"
	"strconv"
	"strings"
	"testing"
	"time"
)

func TestParseTLE(t *testing.T) {
	// This TLE is from the user's original test, corresponding to the Sattrack output provided.
	// Line 1 Checksum: 4 (Calculated: 4)
	// Line 2 Checksum: 3 (Calculated: 2 - Note: The TLE string's checksum for L2 is non-standard)
	issTLE := `1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994
2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.49587957510533`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		// If the checksum for Line 2 is indeed '3' in the TLE and our standard calculation gives '2',
		// this ParseTLE will fail here if the fixed calculateChecksum is used.
		// This is the behavior indicated by "TestParseTLE is failing because of the checksum not passing".
		// The test setup itself will cause a failure if tle.CheckSum2 (3) != calculated (2).
		t.Fatalf("Failed to parse ISS TLE: %v", err)
	}

	// Expected values from Sattrack output for the TLE above
	// Epoch: 2025-05-18 08:53:29.535936 UTC
	expectedEpoch := time.Date(2025, 5, 18, 8, 53, 29, 535936000, time.UTC)

	tests := []struct {
		name    string
		got     interface{}
		want    interface{}
		epsilon float64 // Epsilon for float comparison, or nanoseconds for time comparison
		compare func(got, want interface{}, epsilon float64) bool
	}{
		// Line 1 fields
		{"Satellite Number", tle.SatelliteNumber, 25544, 0, compareExact},
		{"Classification", string(tle.Classification), "U", 0, compareExact},
		{"International ID", tle.International, "98067A", 0, compareExact},
		{"Epoch Year", tle.EpochYear, 2025, 0, compareExact},
		{"Epoch Day", tle.EpochDay, 138.37048074, 1e-9, compareFloat}, // Increased precision for EpochDay
		{"Mean Motion Dot", tle.MeanMotionDot, 0.00007749, 1e-9, compareFloat},
		{"Mean Motion Dot2", tle.MeanMotionDot2, 0.0, 1e-9, compareFloat},
		{"B* Drag Term", tle.Bstar, 0.00014567, 1e-9, compareFloat}, // B* from " 14567-3"
		{"Element Number", tle.ElementNumber, 999, 0, compareExact},
		{"Checksum1", tle.CheckSum1, 4, 0, compareExact},                 // From TLE string
		{"EpochTime", tle.EpochTime(), expectedEpoch, 1000, compareTime}, // Epsilon 1000ns = 1µs for time

		// Line 2 fields
		{"Inclination", tle.Inclination, 51.6369, 1e-4, compareFloat},
		{"Right Ascension", tle.RightAscension, 94.7823, 1e-4, compareFloat},
		{"Eccentricity", tle.Eccentricity, 0.0002558, 1e-7, compareFloat}, // From "0002558"
		{"Argument of Perigee", tle.ArgOfPerigee, 120.7586, 1e-4, compareFloat},
		{"Mean Anomaly", tle.MeanAnomaly, 15.7840, 1e-4, compareFloat},
		{"Mean Motion", tle.MeanMotion, 15.49587957, 1e-8, compareFloat},
		{"Revolution Number", tle.RevolutionNumber, 51053, 0, compareExact},
		{"Checksum2", tle.CheckSum2, 3, 0, compareExact}, // From TLE string (note: standard calc yields 2)
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if !tt.compare(tt.got, tt.want, tt.epsilon) {
				t.Errorf("%s = %v, want %v (epsilon: %g)", tt.name, tt.got, tt.want, tt.epsilon)
			}
		})
	}
}

func compareExact(got, want interface{}, _ float64) bool {
	return got == want
}

func compareFloat(got, want interface{}, epsilon float64) bool {
	g, ok1 := got.(float64)
	w, ok2 := want.(float64)
	if !ok1 || !ok2 {
		// Allow int to float comparison for convenience, e.g. Bstar can be 0
		gInt, okGInt := got.(int)
		wInt, okWInt := want.(int)
		if okGInt && okWInt {
			return math.Abs(float64(gInt)-float64(wInt)) < epsilon
		}
		if okGInt {
			g = float64(gInt)
			ok1 = true
		}
		if okWInt {
			w = float64(wInt)
			ok2 = true
		}
		if !ok1 || !ok2 {
			return false
		}
	}
	return math.Abs(g-w) < epsilon
}

func compareTime(got, want interface{}, epsilonNano float64) bool {
	g, ok1 := got.(time.Time)
	w, ok2 := want.(time.Time)
	if !ok1 || !ok2 {
		return false
	}
	diffNano := math.Abs(float64(g.UnixNano() - w.UnixNano()))
	return diffNano < epsilonNano
}

func TestInvalidTLE(t *testing.T) {
	tests := []struct {
		name    string
		tle     string
		wantErr bool
	}{
		{
			name:    "Empty input",
			tle:     "",
			wantErr: true,
		},
		{
			name:    "Single line",
			tle:     "1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996",
			wantErr: true,
		},
		{
			name: "Invalid line length",
			tle: `1 25544U 98067A   25025.00048859
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`,
			wantErr: true,
		},
		{
			name: "Invalid line numbers", // Line 1 starts with '3'
			tle: `3 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`,
			wantErr: true,
		},
		{
			name: "Line 2 satellite number mismatch",
			tle: `1 25544U 98067A   24001.00000000  .00000000  00000-0  00000-0 0  001Checksum1
2 00002  51.6000 000.0000 0000000 000.0000 000.0000 15.0000000000000Checksum2`,
			wantErr: true, // Satellite numbers 25544 vs 00002
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Replace Checksum1 and Checksum2 placeholders with actual valid checksums if needed for these specific error tests
			// For now, these TLEs are malformed enough that they should fail before or at checksum.
			// Let's make the checksums "valid" for the dummy data to ensure other errors are caught.
			// This is a bit of a workaround for generic error tests.
			// For the purpose of these tests, we only care that *an* error is returned.

			testTLE := tt.tle
			if strings.Contains(testTLE, "Checksum1") {
				l1 := strings.Split(testTLE, "\n")[0]
				l1 = l1[:68]
				c1, _ := calculateChecksum(l1 + "0")
				testTLE = strings.Replace(testTLE, "Checksum1", strconv.Itoa(c1), 1)
			}
			if strings.Contains(testTLE, "Checksum2") {
				l2 := strings.Split(testTLE, "\n")[1]
				l2 = l2[:68]
				c2, _ := calculateChecksum(l2 + "0")
				testTLE = strings.Replace(testTLE, "Checksum2", strconv.Itoa(c2), 1)
			}

			_, err := ParseTLE(testTLE)
			if (err != nil) != tt.wantErr {
				t.Errorf("ParseTLE() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}
