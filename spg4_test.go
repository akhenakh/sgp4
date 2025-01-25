package spg4

import (
	"math"
	"testing"
)

func TestSGP4Propagator(t *testing.T) {
	// ISS TLE from the previous test
	issTLE := `1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	// Test position at epoch (t = 0)
	pos, err := tle.GetPosition(0.0)
	if err != nil {
		t.Fatalf("Failed to get position: %v", err)
	}

	// Verify basic sanity checks
	// ISS orbits at roughly 400km altitude
	radius := math.Sqrt(pos.X*pos.X + pos.Y*pos.Y + pos.Z*pos.Z)
	if radius < 6700 || radius > 6800 { // Expected range considering ISS altitude
		t.Errorf("Calculated radius %f km is outside expected range for ISS", radius)
	}

	// Verify velocity magnitude (~7.7 km/s for LEO)
	velocity := math.Sqrt(pos.VX*pos.VX + pos.VY*pos.VY + pos.VZ*pos.VZ)
	if velocity < 7.5 || velocity > 7.9 {
		t.Errorf("Calculated velocity %f km/s is outside expected range for ISS", velocity)
	}

	// Test propagation over time
	timesteps := []float64{0.0, 10.0, 30.0, 60.0} // minutes since epoch
	for _, dt := range timesteps {
		pos, err := tle.GetPosition(dt)
		if err != nil {
			t.Errorf("Failed to get position at t=%f: %v", dt, err)
			continue
		}

		// Basic sanity checks
		radius := math.Sqrt(pos.X*pos.X + pos.Y*pos.Y + pos.Z*pos.Z)
		if radius < 6700 || radius > 6800 {
			t.Errorf("At t=%f minutes: radius %f km is outside expected range", dt, radius)
		}

		velocity := math.Sqrt(pos.VX*pos.VX + pos.VY*pos.VY + pos.VZ*pos.VZ)
		if velocity < 7.5 || velocity > 7.9 {
			t.Errorf("At t=%f minutes: velocity %f km/s is outside expected range", dt, velocity)
		}
	}
}
