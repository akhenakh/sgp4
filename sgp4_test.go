package sgp4

import (
	"errors"
	"math"
	"testing"
	"time"
)

func TestGetPositionWithLatLngAlt(t *testing.T) {
	// ISS TLE
	issTLE := `1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994
2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.49587957510533`

	/*
	   Norad Number:         25544
	   Int. Designator:      98067A
	   Epoch:                2025-05-18 08:53:29.535936 UTC
	   Orbit Number:         51053
	   Mean Motion Dt2:        0.00007749
	   Mean Motion Ddt6:       0.00000000
	   Eccentricity:           0.00025580
	   BStar:                  0.00014567
	   Inclination:           51.63690000
	   Right Ascending Node:  94.78230000
	   Argument Perigee:     120.75860000
	   Mean Anomaly:          15.78400000
	   Mean Motion:           15.49587957

	   2025-05-18 08:53:29.535936 UTC Az:  287.678, El:  -13.965, Rng:   4348.664, Rng Rt:  -4.246 Lat:   32.740, Lon: -125.293, Alt:    418.256
	*/
	// Parse TLE
	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	// Get latitude, longitude, and altitude
	eci, err := tle.FindPosition(0.0)
	if err != nil {
		t.Fatalf("Failed to get position: %v", err)
	}

	lat, lng, alt := eci.ToGeodetic()

	// Expected values (verified)
	expectedLat := 32.740
	expectedLng := -125.293
	expectedAlt := 418.256

	// Tolerance for comparison
	const tolerance = 0.005

	// Compare latitude
	if math.Abs(lat-expectedLat) > tolerance {
		t.Errorf("Latitude mismatch: got %.3f, want %.3f (±%.3f)", lat, expectedLat, tolerance)
	}

	// Compare longitude
	if math.Abs(lng-expectedLng) > tolerance {
		t.Errorf("Longitude mismatch: got %.3f, want %.3f (±%.3f)", lng, expectedLng, tolerance)
	}

	// Compare altitude
	if math.Abs(alt-expectedAlt) > tolerance {
		t.Errorf("Altitude mismatch: got %.3f, want %.3f (±%.3f)", alt, expectedAlt, tolerance)
	}
}

// Helper function to convert Geodetic (Lat, Lon, Alt in degrees/km) to ECI
// for a given DateTime. This is essentially the inverse of ToGeodetic for testing.
func geodeticToEciForTest(latDeg, lonDeg, altKm float64, dt time.Time, R_earth, flattening float64) Vector {
	latRad := latDeg * deg2rad
	lonRad := lonDeg * deg2rad

	e2 := flattening * (2.0 - flattening)
	sinLat := math.Sin(latRad)
	cosLat := math.Cos(latRad)

	var N float64
	if math.Abs(1.0-e2*sinLat*sinLat) < 1e-14 {
		N = R_earth / math.Sqrt(1e-14)
	} else {
		N = R_earth / math.Sqrt(1.0-e2*sinLat*sinLat)
	}

	// ECEF coordinates
	ecefX := (N + altKm) * cosLat * math.Cos(lonRad)
	ecefY := (N + altKm) * cosLat * math.Sin(lonRad)
	ecefZ := (N*(1.0-e2) + altKm) * sinLat

	// Rotate ECEF to ECI
	tempEci := Eci{DateTime: dt}
	gmst := tempEci.GreenwichSiderealTime() // radians

	eciX := ecefX*math.Cos(gmst) - ecefY*math.Sin(gmst)
	eciY := ecefX*math.Sin(gmst) + ecefY*math.Cos(gmst)
	eciZ := ecefZ

	return Vector{X: eciX, Y: eciY, Z: eciZ}
}

func TestToGeodetic(t *testing.T) {
	// Use a fixed DateTime for all test cases to make GMST predictable. J2000.0 is a good choice.
	j2000 := time.Date(2000, 1, 1, 12, 0, 0, 0, time.UTC)

	// Constants used by ToGeodetic (these must match what ToGeodetic actually uses)
	const R_earth_for_test = reSGP4
	const flattening_for_test = fSGP4

	tests := []struct {
		name        string
		geodeticIn  struct{ lat, lon, alt float64 } // Input geodetic to derive ECI
		expectedLat float64                         // Expected output Latitude
		expectedLng float64                         // Expected output Longitude
		expectedAlt float64                         // Expected output Altitude
		tolerance   float64
	}{
		{
			name:        "Equatorial position, 0 longitude",
			geodeticIn:  struct{ lat, lon, alt float64 }{0.0, 0.0, 0.0},
			expectedLat: 0.0,
			expectedLng: 0.0,
			expectedAlt: 0.0,
			tolerance:   0.0001, // Should be very precise
		},
		{
			name:        "Equatorial position, 90E longitude",
			geodeticIn:  struct{ lat, lon, alt float64 }{0.0, 90.0, 0.0},
			expectedLat: 0.0,
			expectedLng: 90.0,
			expectedAlt: 0.0,
			tolerance:   0.0001,
		},
		{
			name:        "North Polar position",
			geodeticIn:  struct{ lat, lon, alt float64 }{90.0, 0.0, 700.0}, // Lon can be anything at pole
			expectedLat: 90.0,
			expectedLng: 0.0, // atan2(0,0) in ECI often yields 0, then GMST subtraction. Result can be arbitrary.
			expectedAlt: 700.0,
			tolerance:   0.0001,
		},
		{
			name:        "South Polar position",
			geodeticIn:  struct{ lat, lon, alt float64 }{-90.0, 180.0, 200.0}, // Lon can be anything at pole
			expectedLat: -90.0,
			expectedLng: 180.0, // Similar to North Pole, output lon can be arbitrary.
			expectedAlt: 200.0,
			tolerance:   0.0001,
		},
		{
			name:        "Mid-latitude position",
			geodeticIn:  struct{ lat, lon, alt float64 }{34.35, 46.30, 100.0},
			expectedLat: 34.35,
			expectedLng: 46.30,
			expectedAlt: 100.0,
			tolerance:   0.0001,
		},
		{
			name:        "Mid-latitude position, negative longitude",
			geodeticIn:  struct{ lat, lon, alt float64 }{-22.5, -75.25, 50.5},
			expectedLat: -22.5,
			expectedLng: -75.25,
			expectedAlt: 50.5,
			tolerance:   0.0001,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Derive ECI input from the geodeticIn for the fixed J2000 time
			inputEciPos := geodeticToEciForTest(tt.geodeticIn.lat, tt.geodeticIn.lon, tt.geodeticIn.alt, j2000, R_earth_for_test, flattening_for_test)
			eciArg := Eci{DateTime: j2000, Position: inputEciPos}

			lat, lng, alt := eciArg.ToGeodetic()

			if math.Abs(lat-tt.expectedLat) > tt.tolerance {
				t.Errorf("Latitude = %.6f, want %.6f (Δ%.6f, ±%.6f)", lat, tt.expectedLat, lat-tt.expectedLat, tt.tolerance)
			}

			// For longitude, special handling for poles and wrap-around
			if math.Abs(tt.expectedLat) < 89.999 { // If not at a pole
				deltaLng := math.Abs(lng - tt.expectedLng)
				if deltaLng > 180.0 { // handles wrap-around (e.g. -179 vs 179)
					deltaLng = 360.0 - deltaLng
				}
				if deltaLng > tt.tolerance {
					t.Errorf("Longitude = %.6f, want %.6f (Δ%.6f wrapped, ±%.6f)", lng, tt.expectedLng, deltaLng, tt.tolerance)
				}
			} else {
				// At poles, longitude is not well-defined, so we don't strictly check it.
				// Or we could check if it's within a very wide range or a conventional value if any.
				// For now, skipping strict check.
				t.Logf("Skipping strict longitude check at pole for %s (Lat: %.6f)", tt.name, lat)
			}

			if math.Abs(alt-tt.expectedAlt) > tt.tolerance {
				t.Errorf("Altitude = %.6f, want %.6f (Δ%.6f, ±%.6f)", alt, tt.expectedAlt, alt-tt.expectedAlt, tt.tolerance)
			}
		})
	}
}

func TestSatelliteDecay_WithKnownDecayingTLE(t *testing.T) {
	// TLE for STARLINK-1838 (NORAD ID 47129)
	// Epoch: 2025 day 128.113... -> 2025-05-08 ~02:42 UTC
	// Reported decay on 2025-05-08 (same day as TLE epoch).
	// This TLE has extremely high B* (4.7395E-4) and Mean Motion Dot (0.1945 rev/day^2),
	// indicating very rapid decay expected shortly after the TLE epoch, likely within hours.
	// SGP4 often hits eccentricity limits before rk_spp < 1.0 for such cases.
	// Checksums: Line 1: 6, Line 2: 5 (Verified to be correct)
	decayingTLEStr := `STARLINK-1838
1 47129U 20088H   25128.11303883  .19452873  12553-4  47395-3 0  9996
2 47129  53.0081 164.1277 0004123 289.2946 169.2144 16.45991982247015`

	tle, err := ParseTLE(decayingTLEStr)
	if err != nil {
		t.Fatalf("Failed to parse TLE for decay test: %v. TLE:\n%s", err, decayingTLEStr)
	}

	// Propagation window: from epoch up to ~1.5 days after epoch.
	// Decay/model breakdown is expected very rapidly, within hours.
	timeStepMinutes := 5.0           // 5 minute step
	maxTsinceMinutes := 1.5 * 1440.0 // 1.5 days = 2160 minutes

	foundOrbitFailureIndication := false
	var actualDecayErr *SatelliteDecayedError
	var actualModelLimitsErr *SGP4ModelLimitsError

	for tsince := 0.0; tsince <= maxTsinceMinutes; tsince += timeStepMinutes {
		_, errProp := tle.FindPosition(tsince)
		if errProp != nil {
			if errors.As(errProp, &actualDecayErr) {
				t.Logf("SatelliteDecayedError (rk_spp < 1.0) caught for %s at tsince = %.2f minutes (epoch %s)",
					tle.Name, tsince, tle.EpochTime().Format(time.RFC3339))
				t.Logf("Error details: %v", actualDecayErr)
				foundOrbitFailureIndication = true

				if actualDecayErr.Tsince != tsince {
					t.Errorf("SatelliteDecayedError.Tsince mismatch: got %.2f, want %.2f", actualDecayErr.Tsince, tsince)
				}
				if actualDecayErr.Radius >= 1.0 {
					t.Errorf("SatelliteDecayedError.Radius unexpected: got %.4f, should be < 1.0", actualDecayErr.Radius)
				}
				break
			} else if errors.As(errProp, &actualModelLimitsErr) {
				t.Logf("SGP4ModelLimitsError caught for %s at tsince = %.2f minutes (epoch %s)",
					tle.Name, tsince, tle.EpochTime().Format(time.RFC3339))
				t.Logf("Error details: %v", actualModelLimitsErr)
				foundOrbitFailureIndication = true

				if actualModelLimitsErr.Tsince != tsince {
					t.Errorf("SGP4ModelLimitsError.Tsince mismatch: got %.2f, want %.2f", actualModelLimitsErr.Tsince, tsince)
				}
				// Specific checks based on reason can be added here if needed
				// e.g., if actualModelLimitsErr.Reason == ReasonEccentricityTooLow { ... }
				break
			} else {
				// Any other unexpected error
				t.Fatalf("Unexpected error during propagation for %s at tsince %.2f (epoch %s): %v", tle.Name, tsince, tle.EpochTime().Format(time.RFC3339), errProp)
			}
		}
		if tsince >= maxTsinceMinutes && !foundOrbitFailureIndication {
			break
		}
	}

	if !foundOrbitFailureIndication {
		t.Errorf("No orbit failure indication (SatelliteDecayedError or SGP4ModelLimitsError) was triggered for TLE %s (epoch %s) within %.2f minutes (%.1f days). The model might not predict failure for this TLE in the expected way within this timeframe.", tle.Name, tle.EpochTime().Format(time.RFC3339), maxTsinceMinutes, maxTsinceMinutes/1440.0)
	}
}

// TestNonDecayAtEpoch ensures that a standard, stable TLE at its epoch (tsince=0)
// does not incorrectly trigger a decay error.
func TestNonDecayAtEpoch(t *testing.T) {
	issTLE := `ISS (ZARYA)
1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994
2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.49587957510533`

	stdTle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse standard ISS TLE: %v", err)
	}
	_, err = stdTle.FindPosition(0.0) // At epoch
	if err != nil {
		var decayErrCheck *SatelliteDecayedError
		if errors.As(err, &decayErrCheck) {
			t.Errorf("Standard ISS TLE at tsince=0.0 unexpectedly reported decay: %v", err)
		} else {
			// Other errors at tsince=0 are still problems but not decay specific.
			t.Errorf("Standard ISS TLE at tsince=0.0 reported unexpected error: %v", err)
		}
	}
}
