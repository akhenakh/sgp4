package spg4

import "math"

// Initialize converts TLE data into orbital elements needed for SGP4
func (tle *TLE) Initialize() (*OrbitalElements, error) {
	// Convert degrees to radians
	incl := tle.Inclination * deg2rad
	omega := tle.ArgOfPerigee * deg2rad
	raan := tle.RightAscension * deg2rad
	m := tle.MeanAnomaly * deg2rad
	n0 := tle.MeanMotion * twoPi / minutesPerDay
	ecc := tle.Eccentricity

	// Initialize orbital elements
	elem := &OrbitalElements{
		ecc:   ecc,
		incl:  incl,
		omega: omega,
		raan:  raan,
		m:     m,
		n:     n0,
		bstar: tle.Bstar,
	}

	// Recover original mean motion (n0) and semimajor axis (a0)
	a1 := math.Pow(xke/n0, 2.0/3.0)
	cosio := math.Cos(incl)
	theta2 := cosio * cosio
	x3thm1 := 3.0*theta2 - 1.0
	eosq := ecc * ecc
	betao2 := 1.0 - eosq
	betao := math.Sqrt(betao2)

	// Corrections for mean motion
	del1 := 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2)
	ao := a1 * (1.0 - del1*(0.5*2.0/3.0+del1*(1.0+134.0/81.0*del1)))
	delo := 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2)
	elem.a = ao

	// Final mean motion
	elem.n = n0 / (1.0 + delo)

	return elem, nil
}
