package sgp4

import "math"

// Mathematical and physical constants
const (
	twoPi         = 2 * math.Pi
	deg2rad       = math.Pi / 180.0
	rad2deg       = 180.0 / math.Pi
	xkmper        = 6378.135       // Earth's radius in km (SGP4 context, kXKMPER)
	ae            = 1.0            // Distance units/earth radii
	xj2           = 0.001082616    // J2 harmonic (kXJ2)
	xj3           = -0.00000253881 // J3 harmonic (kXJ3)
	xj4           = -0.00000165597 // J4 harmonic (kXJ4)
	torad         = math.Pi / 180.0
	minutesPerDay = 1440.0
	we            = 7.2921150e-5 // Earth's angular velocity (rad/sec) - Not directly kOMEGA_E
	// kOMEGA_E from libsgp4 is earth rotation per sidereal day (1.00273790934)
	// we seems to be omega_earth in rad/sec. SGP4 uses its own set of time/rotation constants.

	// Constants for ECI to Geodetic conversion, aligned with libsgp4 for comparison
	reSGP4 = 6378.135     // Earth's equatorial radius in km (matches kXKMPER for consistency)
	fSGP4  = 1.0 / 298.26 // Earth's flattening factor (matches kF for consistency)

	// WGS-84 Earth model constants (can be kept for other uses if needed, but ToGeodetic will use SGP4 ones)
	reWGS84 = 6378.137            // Earth's equatorial radius in km
	fWGS84  = 1.0 / 298.257223563 // Earth's flattening factor

	// WGS-84 Earth model constants
	re = 6378.137            // Earth's equatorial radius in km
	f  = 1.0 / 298.257223563 // Earth's flattening factor
)

// Computed values (non-constants) - these are based on SGP4's specific constants like xkmper
var (
	xke = 60.0 / math.Sqrt(xkmper*xkmper*xkmper/398600.8) // sqrt(GM_sgp4_units/R_sgp4_units^3)
	// libsgp4 kMU = 398600.8, kXKMPER = 6378.135
	ck2    = 0.5 * xj2 * ae * ae                // kCK2
	ck4    = -0.375 * xj4 * ae * ae * ae * ae   // kCK4
	qoms2t = math.Pow((120.0-78.0)/xkmper, 4.0) // kQOMS2T, based on kQ0=120, kS0=78
	s      = ae * (1.0 + 78.0/xkmper)           // kS, based on kS0=78
	a3ovk2 = -xj3 / ck2                         // kA3OVK2 is -kXJ3 / kCK2 * kAE^3, but since ae=1, this matches if xj3 is scaled by ae^3.
	// libsgp4 kA3OVK2 = -kXJ3 / kCK2 * kAE * kAE * kAE. Your xj3 is already a non-dim harmonic.
	// So a3ovk2 = -xj3 / (0.5 * xj2 * ae*ae) * ae * ae * ae = -xj3 / (0.5 * xj2) * ae
	// If xj3 in constants.go is truly dimensionless kXJ3, then a3ovk2 = -xj3 / (0.5*xj2). My original code had this right.
	// Check libsgp4: kA3OVK2 = -kXJ3 / kCK2 * kAE * kAE * kAE
	// kCK2 = 0.5 * kXJ2 * kAE * kAE
	// So kA3OVK2 = -kXJ3 / (0.5 * kXJ2 * kAE * kAE) * kAE^3 = -kXJ3 / (0.5*kXJ2) * kAE
	// Since ae=1.0, your a3ovk2 = -xj3 / (0.5*xj2) is correct.
)

// Orbital elements intermediate values (output of Initialize)
type OrbitalElements struct {
	a     float64 // Semi-major axis (aodp from SGP4 init, in Earth Radii)
	ecc   float64 // Eccentricity (e0 from TLE)
	incl  float64 // Inclination (i0 from TLE, in rad)
	omega float64 // Argument of perigee (omega_o from TLE, in rad)
	raan  float64 // Right ascension of ascending node (Omega_o from TLE, in rad)
	m     float64 // Mean anomaly (M_o from TLE, in rad)
	n     float64 // Mean motion (no_kozai from SGP4 init, rad/min)
	bstar float64 // Drag term (B* from TLE)

	// SGP4 internal elements after initialization
	// Common Constants
	sinio, cosio float64
	eta          float64
	t2cof        float64
	x1mth2       float64
	x3thm1       float64
	x7thm1       float64 // Not directly in SGP4 common consts, but used in short periodics
	aycof        float64
	xlcof        float64
	xnodcf       float64
	c1, c4       float64
	omgdot       float64 // secular rate of omega (radians/min)
	xnodot       float64 // secular rate of xnode (radians/min)
	xmdot        float64 // secular rate of M (radians/min), this is n_o_kozai

	// NearSpace Constants (for non-deep space)
	c5                  float64
	omgcof              float64
	xmcof               float64
	delmo               float64
	sinmo               float64
	d2, d3, d4          float64
	t3cof, t4cof, t5cof float64

	// Other flags
	isSimpleModel bool // if perigee < 220 km
}

// State vector components
type StateVector struct {
	X, Y, Z    float64 // Position components (km)
	VX, VY, VZ float64 // Velocity components (km/s)
}
