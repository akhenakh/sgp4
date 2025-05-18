package sgp4

import (
	"math"
	"testing"
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

func TestToGeodetic(t *testing.T) {
	tests := []struct {
		name                                  string
		eci                                   Eci
		expectedLat, expectedLng, expectedAlt float64
		tolerance                             float64
	}{
		{
			name: "Equatorial position",
			eci: Eci{
				Position: Vector{X: re, Y: 0, Z: 0},
			},
			expectedLat: 0.0,
			expectedLng: 0.0,
			expectedAlt: 0.0,
			tolerance:   0.1,
		},
		{
			name: "Polar position",
			eci: Eci{
				Position: Vector{X: 0, Y: 0, Z: re + 700},
			},
			expectedLat: 90.0,
			expectedLng: 0.0,
			expectedAlt: 700.0,
			tolerance:   0.1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			lat, lng, alt := tt.eci.ToGeodetic()

			if math.Abs(lat-tt.expectedLat) > tt.tolerance {
				t.Errorf("Latitude = %f, want %f (±%f)", lat, tt.expectedLat, tt.tolerance)
			}
			if math.Abs(lng-tt.expectedLng) > tt.tolerance {
				t.Errorf("Longitude = %f, want %f (±%f)", lng, tt.expectedLng, tt.tolerance)
			}
			if math.Abs(alt-tt.expectedAlt) > tt.tolerance {
				t.Errorf("Altitude = %f, want %f (±%f)", alt, tt.expectedAlt, tt.tolerance)
			}
		})
	}
}
