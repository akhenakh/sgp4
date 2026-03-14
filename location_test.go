package sgp4

import (
	"math"
	"testing"
	"time"
)

// TestAcTanQuadrants verifies AcTan returns correct quadrant angles.
// AcTan returns angles in range [−π/2, 3π/2] which caller may normalize to [0, 2π].
func TestAcTanQuadrants(t *testing.T) {
	tests := []struct {
		name        string
		sinx        float64
		cosx        float64
		expectedDeg float64
		tol         float64
	}{
		// Quadrant I: sin >0, cos >0
		{"Q1", 1.0, 1.0, 45.0, 0.01},
		// Quadrant II: sin >0, cos <0
		{"Q2", 1.0, -1.0, 135.0, 0.01},
		// Quadrant III: sin <0, cos <0
		{"Q3", -1.0, -1.0, 225.0, 0.01},
		// Quadrant IV: sin <0, cos >0 (returns negative angle)
		{"Q4", -1.0, 1.0, -45.0, 0.01},
		// Cardinal directions
		{"North", 0.0, 1.0, 0.0, 0.01},
		{"East", 1.0, 0.0, 90.0, 0.01},
		{"South", 0.0, -1.0, 180.0, 0.01},
		{"West", -1.0, 0.0, 270.0, 0.01},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			angle := AcTan(tt.sinx, tt.cosx)
			angleDeg := angle * 180.0 / math.Pi
			if math.Abs(angleDeg-tt.expectedDeg) > tt.tol {
				t.Errorf("AcTan(%.2f, %.2f) = %.2f°, want %.2f°",
					tt.sinx, tt.cosx, angleDeg, tt.expectedDeg)
			}
		})
	}
}

// TestAcTanAzimuthCalculation tests the specific azimuth calculation formula.
// This directly tests that AcTan(topE, -topS) produces correct azimuth values.
// The bug was: AcTan(-topE, topS) should be AcTan(topE, -topS).
func TestAcTanAzimuthCalculation(t *testing.T) {
	// For SEZ coordinates:
	// topS = South component (positive toward South)
	// topE = East component (positive toward East)
	// topZ = Zenith component (positive upward)
	//
	// Azimuth = AcTan(topE, -topS) normalized to [0, 360)
	//
	// Test cases: (topE, topS) -> expected azimuth direction
	// - (0, >0): satellite South -> azimuth 180°
	// - (0, <0): satellite North -> azimuth 0°
	// - (>0, 0): satellite East -> azimuth 90°
	// - (<0, 0): satellite West -> azimuth 270°
	// - (>0, >0): satellite SE -> azimuth 90°-180°
	// - (<0, >0): satellite SW -> azimuth 180°-270°
	// - (>0, <0): satellite NE -> azimuth 0°-90°
	// - (<0, <0): satellite NW -> azimuth 270°-360°

	tests := []struct {
		name       string
		topE       float64
		topS       float64
		expectedAz float64
		tol        float64
	}{
		// Satellite directly South (topS >0, topE =0)
		{"South", 0.0, 1.0, 180.0, 0.1},
		// Satellite directly North (topS <0, topE =0)
		{"North", 0.0, -1.0, 0.0, 0.1},
		// Satellite directly East (topE >0, topS =0)
		{"East", 1.0, 0.0, 90.0, 0.1},
		// Satellite directly West (topE <0, topS =0)
		{"West", -1.0, 0.0, 270.0, 0.1},
		// Satellite SE quadrant (topE >0, topS >0)
		{"SE", 1.0, 1.0, 135.0, 0.1},
		// Satellite SW quadrant (topE <0, topS >0)
		{"SW", -1.0, 1.0, 225.0, 0.1},
		// Satellite NE quadrant (topE >0, topS <0)
		{"NE", 1.0, -1.0, 45.0, 0.1},
		// Satellite NW quadrant (topE <0, topS <0)
		{"NW", -1.0, -1.0, 315.0, 0.1},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Azimuth = AcTan(topE, -topS) normalized to [0, 360)
			azRad := AcTan(tt.topE, -tt.topS)
			azDeg := azRad * 180.0 / math.Pi
			if azDeg < 0 {
				azDeg += 360.0
			}

			if math.Abs(azDeg-tt.expectedAz) > tt.tol {
				t.Errorf("Azimuth(topE=%.2f, topS=%.2f) = %.2f°, want %.2f°",
					tt.topE, tt.topS, azDeg, tt.expectedAz)
			}
		})
	}
}

// TestAzimuthAgainstCppReference validates azimuth values against C++ implementation.
// The C++ code uses: az = atan(-top_e / top_s); if (top_s > 0) az += PI; if (az < 0) az += 2*PI;
// The Go code should produce identical results with AcTan(topE, -topS).
func TestAzimuthAgainstCppReference(t *testing.T) {
	// Test cases verify equivalence between:
	// C++: az = atan(-top_e/top_s); if (top_s > 0) az += PI; if (az < 0) az += 2*PI
	// Go: az = AcTan(topE, -topS); if (az < 0) az += 2*PI
	tests := []struct {
		name  string
		topS  float64
		topE  float64
		cppAz float64 // Expected azimuth from C++ formula
	}{
		// Satellite directly South: topS >0, topE = 0
		// C++: az = atan(0) = 0; topS >0 so az += PI; result: 180°
		{"South", 1.0, 0.0, 180.0},
		// Satellite directly North: topS <0, topE = 0
		// C++: az = atan(0) = 0; topS <=0 so no PI added; result: 0°
		{"North", -1.0, 0.0, 0.0},
		// Satellite directly East: topS =0, topE >0
		// AcTan(topE, -topS) = AcTan(1, 0) = PI/2 = 90°
		{"East", 0.0, 1.0, 90.0},
		// Satellite directly West: topS =0, topE <0
		// AcTan(-1, 0) = 3PI/2 = 270°
		{"West", 0.0, -1.0, 270.0},
		// Satellite SE: topS >0, topE >0
		// C++: az = atan(-topE/topS) = atan(-1); topS >0 so +PI; result: 135°
		// Go: AcTan(topE, -topS) = AcTan(1, -1) = PI + atan(-1) = 135°
		{"SE", 1.0, 1.0, 135.0},
		// Satellite SW: topS >0, topE <0
		// C++: az = atan(1) = PI/4; topS >0 so +PI = 5PI/4 = 225°
		// Go: AcTan(topE, -topS) = AcTan(-1, -1) = PI + atan(1) = 225°
		{"SW", 1.0, -1.0, 225.0},
		// Satellite NE: topS <0, topE >0
		// C++: az = atan(-topE/topS) = atan(-1/-1) = atan(1) = PI/4 = 45°
		// Go: AcTan(topE, -topS) = AcTan(1, 1) = atan(1) = 45°
		{"NE", -1.0, 1.0, 45.0},
		// Satellite NW: topS <0, topE <0
		// C++: az = atan(-topE/topS) = atan(1/-1) = atan(-1) = -PI/4; az < 0 so +2PI = 315°
		// Go: AcTan(topE, -topS) = AcTan(-1, 1) = atan(-1) = -45°; normalize = 315°
		{"NW", -1.0, -1.0, 315.0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Go formula
			azRad := AcTan(tt.topE, -tt.topS)
			azDeg := azRad * 180.0 / math.Pi
			if azDeg < 0 {
				azDeg += 360.0
			}

			if math.Abs(azDeg-tt.cppAz) > 0.1 {
				t.Errorf("Azimuth mismatch: got %.2f°, want %.2f° (topS=%.2f, topE=%.2f)",
					azDeg, tt.cppAz, tt.topS, tt.topE)
			}
		})
	}
}

// TestGetLookAngleAzimuthReference validates azimuth values against known geometry.
// This test uses reference values computed with the correct formula AcTan(topE, -topS).
// If the formula is changed to AcTan(-topE, topS), the azimuth values will be wrong
// and this test will fail.
//
// The azimuth formula must map SEZ coordinates as follows:
//   - Satellite South of observer (topS > 0) → azimuth in [90°, 270°]
//   - Satellite North of observer (topS < 0) → azimuth in [0°, 90°] ∪ [270°, 360°]
//
// The buggy formula AcTan(-topE, topS) would incorrectly map:
//   - Satellite South (topS > 0) → azimuth in [0°, 180°] ← WRONG
//   - Satellite North (topS < 0) → azimuth in [180°, 360°] ← WRONG
func TestGetLookAngleAzimuthReference(t *testing.T) {
	// ISS TLE from TestGeneratePasses_Reference (validated against C++)
	issTLE := `1 25544U 98067A   25247.10182809  .00011777  00000-0  21333-3 0  9997
2 25544  51.6327 275.9345 0004179 299.5263  60.5309 15.50088696528307`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	// Observer in Montreal (from TestGeneratePasses_Reference)
	observer := &Location{
		Latitude:  45.51,
		Longitude: -73.59,
		Altitude:  60.0,
	}

	// Reference azimuth values computed with the correct formula AcTan(topE, -topS).
	// These values were validated to be consistent with C++ libsgp4 implementation.
	// During this pass, the ISS moves from SW to NE, meaning:
	//   - At AOS: satellite is SW of observer, azimuth ~210°
	//   - At mid-pass: satellite is SE of observer, azimuth ~150°
	//   - At LOS: satellite is NE of observer, azimuth ~70°
	testCases := []struct {
		name        string
		timeUTC     string
		expectedAz  float64
		expectedEl  float64
		azTolerance float64
		elTolerance float64
	}{
		{
			name:        "Pass_AOS_Time",
			timeUTC:     "2025-09-11T00:09:30Z",
			expectedAz:  211.75, // Satellite SW of observer
			expectedEl:  1.95,
			azTolerance: 1.0,
			elTolerance: 0.5,
		},
		{
			name:        "Pass_Mid_Time",
			timeUTC:     "2025-09-11T00:14:00Z",
			expectedAz:  146.63, // Satellite SE of observer
			expectedEl:  29.72,
			azTolerance: 1.0,
			elTolerance: 0.5,
		},
		{
			name:        "Pass_LOS_Time",
			timeUTC:     "2025-09-11T00:19:00Z",
			expectedAz:  67.93, // Satellite NE of observer
			expectedEl:  1.47,
			azTolerance: 1.0,
			elTolerance: 0.5,
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			testTime, err := time.Parse(time.RFC3339, tc.timeUTC)
			if err != nil {
				t.Fatalf("Failed to parse time: %v", err)
			}

			epoch := tle.EpochTime()
			tsince := testTime.Sub(epoch).Minutes()

			eci, err := tle.FindPosition(tsince)
			if err != nil {
				t.Fatalf("FindPosition failed: %v", err)
			}

			sv := StateVector{
				X:  eci.Position.X,
				Y:  eci.Position.Y,
				Z:  eci.Position.Z,
				VX: eci.Velocity.X,
				VY: eci.Velocity.Y,
				VZ: eci.Velocity.Z,
			}

			obs, err := sv.GetLookAngle(observer, eci.DateTime)
			if err != nil {
				t.Fatalf("GetLookAngle failed: %v", err)
			}

			// Check azimuth matches expected value
			azDiff := math.Abs(obs.LookAngles.Azimuth - tc.expectedAz)
			if azDiff > tc.azTolerance {
				t.Errorf("Azimuth mismatch: got %.2f°, want %.2f° (diff %.2f°)",
					obs.LookAngles.Azimuth, tc.expectedAz, azDiff)
			}

			// Check elevation matches expected value
			elDiff := math.Abs(obs.LookAngles.Elevation - tc.expectedEl)
			if elDiff > tc.elTolerance {
				t.Errorf("Elevation mismatch: got %.2f°, want %.2f° (diff %.2f°)",
					obs.LookAngles.Elevation, tc.expectedEl, elDiff)
			}
		})
	}
}

// TestGetLookAngleValidatesRange ensures azimuth and elevation are in valid ranges.
func TestGetLookAngleValidatesRange(t *testing.T) {
	// ISS TLE from TestGeneratePasses_Reference
	issTLE := `1 25544U 98067A   25247.10182809  .00011777  00000-0  21333-3 0  9997
2 25544  51.6327 275.9345 0004179 299.5263  60.5309 15.50088696528307`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	observer := &Location{
		Latitude:  45.51,
		Longitude: -73.59,
		Altitude:  60.0,
	}

	// Test at multiple times during a known pass
	testTimes := []string{
		"2025-09-11T00:09:30Z", // Near AOS
		"2025-09-11T00:14:00Z", // Mid pass
		"2025-09-11T00:19:00Z", // Near LOS
	}

	for _, timeStr := range testTimes {
		t.Run(timeStr, func(t *testing.T) {
			testTime, _ := time.Parse(time.RFC3339, timeStr)
			epoch := tle.EpochTime()
			tsince := testTime.Sub(epoch).Minutes()

			eci, err := tle.FindPosition(tsince)
			if err != nil {
				t.Fatalf("FindPosition failed: %v", err)
			}

			sv := StateVector{
				X:  eci.Position.X,
				Y:  eci.Position.Y,
				Z:  eci.Position.Z,
				VX: eci.Velocity.X,
				VY: eci.Velocity.Y,
				VZ: eci.Velocity.Z,
			}

			obs, err := sv.GetLookAngle(observer, eci.DateTime)
			if err != nil {
				t.Fatalf("GetLookAngle failed: %v", err)
			}

			// Azimuth must be in [0, 360)
			if obs.LookAngles.Azimuth < 0 || obs.LookAngles.Azimuth >= 360 {
				t.Errorf("Azimuth %.2f° out of valid range [0, 360)", obs.LookAngles.Azimuth)
			}

			// Elevation must be in [-90, 90]
			if obs.LookAngles.Elevation < -90 || obs.LookAngles.Elevation > 90 {
				t.Errorf("Elevation %.2f° out of valid range [-90, 90]", obs.LookAngles.Elevation)
			}

			// Range must be positive
			if obs.LookAngles.Range <= 0 {
				t.Errorf("Range %.2f km must be positive", obs.LookAngles.Range)
			}

			t.Logf("Time %s: Az=%.2f°, El=%.2f°, Range=%.1f km",
				timeStr, obs.LookAngles.Azimuth, obs.LookAngles.Elevation, obs.LookAngles.Range)
		})
	}
}

// TestGetLookAngleNilLocation tests error handling for nil location.
func TestGetLookAngleNilLocation(t *testing.T) {
	sv := StateVector{X: 1000, Y: 1000, Z: 1000}
	_, err := sv.GetLookAngle(nil, time.Now())
	if err != ErrLocationNil {
		t.Errorf("Expected ErrLocationNil, got %v", err)
	}
}

// TestGetLookAngleInvalidLatitude tests error handling for invalid latitude.
func TestGetLookAngleInvalidLatitude(t *testing.T) {
	sv := StateVector{X: 1000, Y: 1000, Z: 1000}

	// Latitude too high
	loc := &Location{Latitude: 91, Longitude: 0, Altitude: 0}
	_, err := sv.GetLookAngle(loc, time.Now())
	if err != ErrInvalidLocationLatitude {
		t.Errorf("Expected ErrInvalidLocationLatitude for lat=91, got %v", err)
	}

	// Latitude too low
	loc.Latitude = -91
	_, err = sv.GetLookAngle(loc, time.Now())
	if err != ErrInvalidLocationLatitude {
		t.Errorf("Expected ErrInvalidLocationLatitude for lat=-91, got %v", err)
	}
}
