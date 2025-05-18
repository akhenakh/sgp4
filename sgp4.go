package sgp4

import (
	"fmt"
	"math"
	"time"
)

type Eci struct {
	DateTime time.Time
	Position Vector
	Velocity Vector
}

type Vector struct {
	X, Y, Z float64
}

// ToGeodetic converts ECI coordinates to geodetic coordinates (lat, lon, alt)
func (eci *Eci) ToGeodetic() (lat, lon, alt float64) {
	// Constants from package constants.go (re, f)
	e2 := f * (2.0 - f) // Square of eccentricity

	// Get GMST for the epoch
	gmst := eci.GreenwichSiderealTime()

	x := eci.Position.X
	y := eci.Position.Y
	z := eci.Position.Z

	// Calculate longitude
	lon = math.Atan2(y, x) - gmst
	lon = wrapLongitude(lon) // Ensure longitude is in [-Pi, Pi]

	// Calculate distance from Earth's center in the XY plane
	r := math.Sqrt(x*x + y*y)

	// Iterative calculation of geodetic latitude
	lat = math.Atan2(z, r) // Initial guess using spherical approximation (geocentric latitude)

	const maxIter = 10
	const tol = 1e-10 // Tolerance for latitude convergence in radians

	var oldLat float64
	var c_iter float64 // Term related to radius of curvature N/a

	for i := 0; i < maxIter; i++ {
		oldLat = lat
		sinLat := math.Sin(lat)
		// c_iter = N/re = 1 / sqrt(1 - e^2 sin^2(phi))
		c_iter = 1.0 / math.Sqrt(1.0-e2*sinLat*sinLat)
		lat = math.Atan2(z+re*c_iter*e2*sinLat, r)

		if math.Abs(lat-oldLat) < tol {
			break
		}
	}

	// Calculate height above ellipsoid
	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	// N_val is N = re * c_iter (Radius of curvature in prime vertical)
	N_val := re * (1.0 / math.Sqrt(1.0-e2*sinLat*sinLat)) // Recompute N_val with final lat

	if math.Abs(cosLat) < 1e-10 { // At poles or very close
		alt = math.Abs(z) - re*math.Sqrt(1.0-e2) // alt = |Z| - b
	} else {
		alt = r/cosLat - N_val
	}

	// Convert to degrees
	lat = lat * rad2deg
	lon = lon * rad2deg

	return lat, lon, alt
}

func (tle *TLE) FindPosition(tsince float64) (Eci, error) {
	// Step 1: Get the SGP4-initialized orbital elements
	elems, err := tle.Initialize()
	if err != nil {
		return Eci{}, fmt.Errorf("failed to initialize TLE for FindPosition: %w", err)
	}

	// Use elements from 'elems' struct:
	ao := elems.a         // ao (aodp for near-earth): semi-major axis corrected for J2 (ER)
	ecco := elems.ecc     // eo: eccentricity from TLE (SGP4 uses this as initial eo)
	inclo := elems.incl   // io: inclination in radians
	omegao := elems.omega // omega_o: arg of perigee in radians
	nodeo := elems.raan   // RAAN_o: RAAN in radians
	mo := elems.m         // M_o: mean anomaly in radians
	no := elems.n         // no_kozai: mean motion in rad/min, corrected for J2

	cosio := math.Cos(inclo)
	sinio := math.Sin(inclo)

	// Secular J2 Perturbations (Simplified Model - does not include J3, drag, etc.)
	// These rates are per minute since 'no' is in rad/min.
	// ck2 = 0.5 * xj2 * ae * ae (from constants.go)
	// p = ao * (1 - ecco*ecco)
	// p^2 = ao^2 * (1 - ecco^2)^2
	// Note: SGP4 defines betao = sqrt(1-ecco*ecco), so p = ao * betao^2. p^2 = ao^2 * betao^4
	// betao_sq := 1.0 - ecco*ecco // betao^2
	// p_sq := ao * ao * betao_sq * betao_sq // (ao * betao^2)^2 = ao^2 * betao^4
	// SGP4 secular rates are often simplified using specific SGP4 terms.
	// Let's use simplified Vallado-style J2 rates for Omega_dot and omega_dot.
	// x3thm1 is (3cos^2i - 1). x1mth2 is (1-cos^2i) = sin^2i
	// SGP4 uses tsquare = cosio*cosio, then x3thm1 = 3*tsquare-1, x1mth2 = 1-tsquare
	// For Omega_dot: -1.5 * no * ck2 * cosio / (ao^2 * betao^4) -- from SGP4-like derivations
	// For omega_dot:  0.75 * no * ck2 * (4 - 5*sinio^2) / (ao^2 * betao^4) -- or (3cosio^2-1) terms

	// SGP4 specific secular terms for node and argp (simplified, omitting J3 influence):
	// From Vallado SGP4 section or Spacetrack Report #3
	// Note: These are simplified. Full SGP4 uses `aodp`, `nodp`, `omgdp` which include more terms.
	// Using `ao` and `ecco` as substitutes for `aodp` and `eodp`.
	// `no` is `nodp` (Kozai mean motion)
	// epochDateTime := tle.EpochTime()
	// tempEciForGmst := Eci{DateTime: epochDateTime}
	// gsto := tempEciForGmst.GreenwichSiderealTime() // GMST at epoch

	// Secular rates based on common SGP4 structure for J2 effects:
	// These are delta_omega, delta_node per period, scaled by mean motion.
	// For rates per minute:
	// C2 = ck2
	// C4 = ck4 (not used in this simplified version)
	// A_30 = -xj3 * ae^3 * sinio (used for J3 terms, not here)

	// Standard J2 secular rates:
	// term_common = no * ck2 / (ao*ao * betao_sq*betao_sq) // no*ck2 / (p0^2) where p0=ao*betao^2
	term_common_nodal_argp := no * ck2 / (elems.a * elems.a * (1 - elems.ecc*elems.ecc) * (1 - elems.ecc*elems.ecc))

	node_dot_rad_per_min := -1.5 * term_common_nodal_argp * cosio
	argp_dot_rad_per_min := 0.75 * term_common_nodal_argp * (4.0 - 5.0*sinio*sinio)

	// Propagate RAAN and Argument of Perigee
	nodek := nodeo + node_dot_rad_per_min*tsince
	omegak := omegao + argp_dot_rad_per_min*tsince

	// Propagate Mean Anomaly
	// M_k = M_o + n_o_kozai * t + (J2 secular effect on M)
	// SGP4's `no` (elems.n here) is Kozai mean motion, already including a primary J2 effect.
	// For a simplified model, M_k = M_o + no_kozai * tsince is a common starting point.
	// More complete SGP4 would add further J2 and J3 secular terms to M's derivative.
	// The original SGP4 code also had a term for `xmdot_val` to adjust `no`.
	// For now:
	xmp := mo + no*tsince // Propagated mean anomaly M_k, using no_kozai

	// Solve Kepler's equation for Eccentric Anomaly E_k
	Ek := xmp // Initial guess for E_k is M_k
	for i := 0; i < 10; i++ {
		prev_Ek := Ek
		Ek = xmp + ecco*math.Sin(Ek)
		if math.Abs(Ek-prev_Ek) < 1e-12 {
			break
		}
	}

	// Calculate True Anomaly (fk) and radius (rk)
	sinEk := math.Sin(Ek)
	cosEk := math.Cos(Ek)

	// True anomaly fk
	fk := math.Atan2(math.Sqrt(1.0-ecco*ecco)*sinEk, cosEk-ecco)

	// Radius rk (osculating, in Earth Radii)
	rk := ao * (1.0 - ecco*cosEk)

	// Argument of Latitude (uk) = omegak + fk
	uk := omegak + fk

	// Position calculations (No short-period perturbations applied here for simplicity yet)
	// If short-period terms were added, they'd perturb rk, uk, nodek, inclo
	// For now, use secularly propagated elements for ECI transformation.

	// Orientation vectors for ECI transformation
	cos_uk := math.Cos(uk)
	sin_uk := math.Sin(uk)
	cos_nodek := math.Cos(nodek)
	sin_nodek := math.Sin(nodek)
	cos_inclo := math.Cos(inclo) // Using initial inclination; SGP4 would perturb this
	sin_inclo := math.Sin(inclo)

	// Position vector components in ECI frame (standard transformation from orbital elements)
	// Px = cos(Omega)cos(u) - sin(Omega)sin(u)cos(i)
	// Py = sin(Omega)cos(u) + cos(Omega)sin(u)cos(i)
	// Pz = sin(u)sin(i)
	// Where u is argument of latitude (omegak + fk)
	pos_ux := cos_nodek*cos_uk - sin_nodek*sin_uk*cos_inclo
	pos_uy := sin_nodek*cos_uk + cos_nodek*sin_uk*cos_inclo
	pos_uz := sin_inclo * sin_uk

	position := Vector{
		X: rk * pos_ux * xkmper, // Convert Earth radii to km
		Y: rk * pos_uy * xkmper,
		Z: rk * pos_uz * xkmper,
	}

	// Velocity calculation (ER/min then converted to km/s)
	// Using osculating elements at time t_k (ao, ecco, no, Ek, rk)

	// Radial velocity (rdot_er_min) in ER/min
	// rdot = n * a * e * sin(E) / (1 - e * cos(E))
	rdot_er_min := (no * ao * ecco * sinEk) / (1.0 - ecco*cosEk)

	// Transverse velocity (rfdot_er_min) in ER/min (this is r * v_theta_dot)
	// v_theta_dot = n * a * sqrt(1-e^2) / r = n * sqrt(1-e^2) / (1 - e*cosE)
	// rfdot = r * v_theta_dot = no * ao * sqrt(1-ecco*ecco) / (1.0 - ecco*cosEk)
	rfdot_er_min := (no * ao * math.Sqrt(1.0-ecco*ecco)) / (1.0 - ecco*cosEk)

	// ECI velocity components (transforming {rdot, rfdot} from orbital frame to ECI)
	// Qx = -cos(Omega)sin(u) - sin(Omega)cos(u)cos(i)
	// Qy = -sin(Omega)sin(u) + cos(Omega)cos(u)cos(i)
	// Qz =  cos(u)sin(i)
	vel_qx := -cos_nodek*sin_uk - sin_nodek*cos_uk*cos_inclo
	vel_qy := -sin_nodek*sin_uk + cos_nodek*cos_uk*cos_inclo
	vel_qz := cos_inclo * cos_uk // Note: sin_inclo * cos_uk in some conventions for Qz of RSW frame

	v_factor := xkmper / 60.0 // To convert ER/min to km/s

	velocity := Vector{
		X: (rdot_er_min*pos_ux + rfdot_er_min*vel_qx) * v_factor,
		Y: (rdot_er_min*pos_uy + rfdot_er_min*vel_qy) * v_factor,
		Z: (rdot_er_min*pos_uz + rfdot_er_min*vel_qz) * v_factor,
	}

	return Eci{
		DateTime: tle.EpochTime().Add(time.Duration(tsince) * time.Minute),
		Position: position,
		Velocity: velocity,
	}, nil
}

// GreenwichSiderealTime calculates the Greenwich Mean Sidereal Time
func (eci *Eci) GreenwichSiderealTime() float64 {
	jd := julianDateTime(eci.DateTime)
	t := (jd - 2451545.0) / 36525.0 // Julian centuries since J2000.0

	// GMST in degrees (IAU 2006 formula, matches Vallado, p. 188, Eq. 3-47)
	gmst_deg := 280.46061837 +
		360.98564736629*(jd-2451545.0) +
		0.000387933*t*t -
		t*t*t/38710000.0

	gmst_deg = math.Mod(gmst_deg, 360.0)
	if gmst_deg < 0 {
		gmst_deg += 360.0
	}
	return gmst_deg * deg2rad // Convert to radians
}

// julianDateTime converts time.Time to Julian Date
func julianDateTime(t_utc time.Time) float64 {
	y := float64(t_utc.Year())
	m := float64(t_utc.Month())
	d := float64(t_utc.Day())
	h := float64(t_utc.Hour())
	min := float64(t_utc.Minute())
	s := float64(t_utc.Second())
	ns := float64(t_utc.Nanosecond())

	if m <= 2 {
		y--
		m += 12
	}

	A_jd := math.Floor(y / 100.0)
	B_jd := 2 - A_jd + math.Floor(A_jd/4.0)

	JD_day := math.Floor(365.25*(y+4716.0)) +
		math.Floor(30.6001*(m+1.0)) +
		d + B_jd - 1524.5

	dayFrac := (h + min/60.0 + s/3600.0 + ns/3.6e12) / 24.0
	return JD_day + dayFrac
}

// wrapLongitude wraps longitude to range -π to π
func wrapLongitude(lon float64) float64 {
	lon = math.Mod(lon, twoPi)
	if lon > math.Pi {
		lon -= twoPi
	} else if lon < -math.Pi {
		lon += twoPi
	}
	return lon
}

func gmsTime(t time.Time) float64 {
	jd := julianDateTime(t)
	centuries := (jd - 2451545.0) / 36525.0
	gmst_secs_at_0h := 24110.54841 + centuries*(8640184.812866+
		centuries*(0.093104-centuries*6.2e-6))
	dayFraction := jd - math.Floor(jd+0.5)
	if dayFraction < 0 {
		dayFraction += 1.0
	}
	ut_seconds := dayFraction * 86400.0
	gmst_total_seconds := gmst_secs_at_0h + ut_seconds*1.00273790935
	gmst_rad := math.Mod(gmst_total_seconds*(twoPi/86164.0905), twoPi)
	if gmst_rad < 0 {
		gmst_rad += twoPi
	}
	return gmst_rad
}
