package spg4

import "math"

// GetPosition calculates satellite position and velocity at given minutes since epoch
func (tle *TLE) GetPosition(tsince float64) (*StateVector, error) {
	elem, err := tle.Initialize()
	if err != nil {
		return nil, err
	}

	// Recover semi-major axis and mean motion
	a := elem.a
	e := elem.ecc
	incl := elem.incl
	omega := elem.omega
	raan := elem.raan
	m := elem.m
	n := elem.n

	// Update for secular gravity and atmospheric drag
	cosio := math.Cos(incl)

	mm := m + n*tsince
	omega_new := omega + (-3.5*ck2*cosio*tsince)/a/a
	raan_new := raan + (-3.5*ck2*cosio*tsince)/a/a

	// Solve Kepler's equation
	u := mm
	for i := 0; i < 10; i++ {
		sin_u := math.Sin(u)
		cos_u := math.Cos(u)
		du := (u - e*sin_u - mm) / (1.0 - e*cos_u)
		u -= du
		if math.Abs(du) < 1e-12 {
			break
		}
	}

	// Short-period periodic elements
	sin_u := math.Sin(u)
	cos_u := math.Cos(u)
	r := a * (1.0 - e*cos_u)
	v := math.Sqrt(398600.8/a) * (e*sin_u/r + 1.0/math.Sqrt(a))

	// Position and velocity in orbital plane
	x := r * cos_u
	y := r * sin_u
	xdot := -sin_u * v
	ydot := cos_u * v

	// Orientation angles for coordinate transformation
	cosOmega := math.Cos(omega_new)
	sinOmega := math.Sin(omega_new)
	sinRaan := math.Sin(raan_new)
	cosRaan := math.Cos(raan_new)
	sinI := math.Sin(incl)
	cosI := math.Cos(incl)

	// Earth-fixed coordinates
	state := &StateVector{
		X: x*(cosOmega*cosRaan-sinOmega*sinRaan*cosI) -
			y*(sinOmega*cosRaan+cosOmega*sinRaan*cosI),
		Y: x*(cosOmega*sinRaan+sinOmega*cosRaan*cosI) +
			y*(-sinOmega*sinRaan+cosOmega*cosRaan*cosI),
		Z: x*(sinOmega*sinI) + y*(cosOmega*sinI),
		VX: xdot*(cosOmega*cosRaan-sinOmega*sinRaan*cosI) -
			ydot*(sinOmega*cosRaan+cosOmega*sinRaan*cosI),
		VY: xdot*(cosOmega*sinRaan+sinOmega*cosRaan*cosI) +
			ydot*(-sinOmega*sinRaan+cosOmega*cosRaan*cosI),
		VZ: xdot*(sinOmega*sinI) + ydot*(cosOmega*sinI),
	}

	return state, nil
}
