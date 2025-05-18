package sgp4

import "math"

// Mathematical and physical constants
const (
	twoPi         = 2 * math.Pi
	deg2rad       = math.Pi / 180.0
	rad2deg       = 180.0 / math.Pi
	xkmper        = 6378.135       // Earth's radius in km
	ae            = 1.0            // Distance units/earth radii
	xj2           = 0.001082616    // J2 harmonic
	xj3           = -0.00000253881 // J3 harmonic
	xj4           = -0.00000165597 // J4 harmonic
	torad         = math.Pi / 180.0
	minutesPerDay = 1440.0
	we            = 7.2921150e-5 // Earth's angular velocity (rad/sec)
	a3ovk2        = -0.5 * xj3 / xj2

	// WGS-84 Earth model constants
	re = 6378.137            // Earth's equatorial radius in km
	f  = 1.0 / 298.257223563 // Earth's flattening factor
)

// Computed values (non-constants)
var (
	xke    = 60.0 / math.Sqrt(xkmper*xkmper*xkmper/398600.8) // sqrt(GM/R³)
	ck2    = 0.5 * xj2 * ae * ae
	ck4    = -0.375 * xj4 * ae * ae * ae * ae
	qoms2t = math.Pow((120-78)/xkmper, 4) // (km/earth radii)^4
	s      = ae * (1.0 + 78.0/xkmper)
)

// Orbital elements intermediate values
type OrbitalElements struct {
	a     float64 // Semi-major axis (km)
	ecc   float64 // Eccentricity
	incl  float64 // Inclination (rad)
	omega float64 // Argument of perigee (rad)
	raan  float64 // Right ascension of ascending node (rad)
	m     float64 // Mean anomaly (rad)
	n     float64 // Mean motion (rad/min)
	bstar float64 // Drag term
}

// State vector components
type StateVector struct {
	X, Y, Z    float64 // Position components (km)
	VX, VY, VZ float64 // Velocity components (km/s)
}
