package spg4

import (
	"math"
	"testing"
)

func TestGetLookAngle(t *testing.T) {
	// Test case using known ISS pass
	sv := &StateVector{
		X:  2328.97048951,
		Y:  -5995.22076416,
		Z:  1719.97067261,
		VX: 2.91207230,
		VY: -0.98341546,
		VZ: -7.09081703,
	}

	// Test location: McDonald Observatory, Texas
	loc := &Location{
		Latitude:  30.6715,
		Longitude: -104.0227,
		Altitude:  2070, // meters
	}

	topo := sv.GetLookAngle(loc)

	// Expected values (calculated using validated external software)
	// Allow for small differences due to different Earth models
	tests := []struct {
		name    string
		got     float64
		want    float64
		epsilon float64
	}{
		{"Azimuth", topo.Azimuth, 98.5, 1.0},
		{"Elevation", topo.Elevation, 23.4, 1.0},
		{"Range", topo.Range, 1245.7, 5.0},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if math.Abs(tt.got-tt.want) > tt.epsilon {
				t.Errorf("%s = %.1f, want %.1f (±%.1f)", tt.name, tt.got, tt.want, tt.epsilon)
			}
		})
	}
}

func TestLocationEdgeCases(t *testing.T) {
	// Test satellite directly overhead
	sv := &StateVector{
		X:  re, // Satellite directly above equator at prime meridian
		Y:  0,
		Z:  0,
		VX: 0,
		VY: 7.5, // Typical LEO velocity
		VZ: 0,
	}

	loc := &Location{
		Latitude:  0,
		Longitude: 0,
		Altitude:  0,
	}

	topo := sv.GetLookAngle(loc)

	if math.Abs(topo.Elevation-90) > 1e-6 {
		t.Errorf("Expected elevation near 90° for overhead pass, got %.2f", topo.Elevation)
	}

	// Test satellite on horizon
	sv = &StateVector{
		X:  0,
		Y:  re + 500, // 500km above surface
		Z:  0,
		VX: -7.5,
		VY: 0,
		VZ: 0,
	}

	loc = &Location{
		Latitude:  0,
		Longitude: 0,
		Altitude:  0,
	}

	topo = sv.GetLookAngle(loc)

	if math.Abs(topo.Elevation) > 1e-6 {
		t.Errorf("Expected elevation near 0° for horizon pass, got %.2f", topo.Elevation)
	}
	if math.Abs(topo.Azimuth-90) > 1e-6 {
		t.Errorf("Expected azimuth of 90° for satellite due East, got %.2f", topo.Azimuth)
	}
}
