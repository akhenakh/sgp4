package sgp4

import "math"

// Initialize converts TLE data into orbital elements and SGP4 constants
func (tle *TLE) Initialize() (*OrbitalElements, error) {
	elem := &OrbitalElements{
		ecc:   tle.Eccentricity,
		incl:  tle.Inclination * deg2rad,
		omega: tle.ArgOfPerigee * deg2rad,
		raan:  tle.RightAscension * deg2rad,
		m:     tle.MeanAnomaly * deg2rad,
		n:     tle.MeanMotion * twoPi / minutesPerDay, // This is TLE n, not no_kozai yet
		bstar: tle.Bstar,
	}

	// Recover original mean motion (no_kozai) and semimajor axis (aodp)
	// This part is from OrbitalElements constructor in libsgp4/OrbitalElements.cc
	// and SGP4::Initialise() in libsgp4/SGP4.cc
	a1 := math.Pow(xke/elem.n, 2.0/3.0) // elem.n is n_tle here
	elem.cosio = math.Cos(elem.incl)
	elem.sinio = math.Sin(elem.incl)
	theta2 := elem.cosio * elem.cosio
	elem.x3thm1 = 3.0*theta2 - 1.0
	elem.x1mth2 = 1.0 - theta2
	elem.x7thm1 = 7.0*theta2 - 1.0 // Used in CalculateFinalPositionVelocity

	eosq := elem.ecc * elem.ecc
	betao2 := 1.0 - eosq
	betao := math.Sqrt(betao2)

	// Similar to OrbitalElements constructor in libsgp4
	temp_init := (1.5 * ck2) * elem.x3thm1 / (betao * betao2)
	del1 := temp_init / (a1 * a1)
	a0 := a1 * (1.0 - del1*(1.0/3.0+del1*(1.0+del1*134.0/81.0)))
	del0 := temp_init / (a0 * a0)

	no_kozai := elem.n / (1.0 + del0) // This is elem.n in SGP4 terms (recovered_mean_motion_)
	aodp := a0 / (1.0 - del0)         // This is elem.a in SGP4 terms (recovered_semi_major_axis_)

	elem.n = no_kozai // Store Kozai mean motion
	elem.a = aodp     // Store SGP4 semi-major axis (in Earth Radii)

	// SGP4 Initialization specific calculations (from SGP4::Initialise)
	// Check for perigee < 156 km (simple model flag is for < 220km)
	perigeeKm := (elem.a*(1.0-elem.ecc) - ae) * xkmper // Perigee height in km
	elem.isSimpleModel = false
	if perigeeKm < 220.0 {
		elem.isSimpleModel = true
	}

	s4 := s          // s is kS from constants.go, (ae * (1.0 + 78.0/xkmper))
	qoms24 := qoms2t // qoms2t is kQOMS2T from constants.go
	if perigeeKm < 156.0 {
		s4_val := perigeeKm - 78.0
		if perigeeKm < 98.0 {
			s4_val = 20.0
		}
		qoms24 = math.Pow((120.0-s4_val)*ae/xkmper, 4.0) // ae is 1.0
		s4 = s4_val/xkmper + ae
	}

	pinvsq := 1.0 / (elem.a * elem.a * betao2 * betao2)
	tsi := 1.0 / (elem.a - s4)
	elem.eta = elem.a * elem.ecc * tsi
	etasq := elem.eta * elem.eta
	eeta := elem.ecc * elem.eta
	psisq := math.Abs(1.0 - etasq) // Should be > 0 if not decayed
	if psisq == 0 {
		psisq = 1e-12
	} // Prevent division by zero if etasq is exactly 1
	coef := qoms24 * math.Pow(tsi, 4.0)
	coef1 := coef / math.Pow(psisq, 3.5)

	c2 := coef1 * elem.n * (elem.a*(1.0+1.5*etasq+eeta*(4.0+etasq)) +
		0.75*ck2*tsi/psisq*elem.x3thm1*(8.0+3.0*etasq*(8.0+etasq)))
	elem.c1 = elem.bstar * c2

	elem.c4 = 2.0 * elem.n * coef1 * elem.a * betao2 *
		(elem.eta*(2.0+0.5*etasq) + elem.ecc*(0.5+2.0*etasq) -
			2.0*ck2*tsi/(elem.a*psisq)*
				(-3.0*elem.x3thm1*(1.0-2.0*eeta+etasq*(1.5-0.5*eeta))+
					0.75*elem.x1mth2*(2.0*etasq-eeta*(1.0+etasq))*math.Cos(2.0*elem.omega)))

	theta4 := theta2 * theta2
	temp1 := 3.0 * ck2 * pinvsq * elem.n
	temp2 := temp1 * ck2 * pinvsq
	temp3 := 1.25 * ck4 * pinvsq * pinvsq * elem.n

	elem.xmdot = elem.n + 0.5*temp1*betao*elem.x3thm1 +
		0.0625*temp2*betao*(13.0-78.0*theta2+137.0*theta4)

	x1m5th := 1.0 - 5.0*theta2
	elem.omgdot = -0.5*temp1*x1m5th +
		0.0625*temp2*(7.0-114.0*theta2+395.0*theta4) +
		temp3*(3.0-36.0*theta2+49.0*theta4)

	xhdot1 := -temp1 * elem.cosio
	elem.xnodot = xhdot1 + (0.5*temp2*(4.0-19.0*theta2)+
		2.0*temp3*(3.0-7.0*theta2))*elem.cosio

	elem.xnodcf = 3.5 * betao2 * xhdot1 * elem.c1
	elem.t2cof = 1.5 * elem.c1

	// For xlcof, aycof (long period periodic coefficients for L')
	// These are from SGP4::RecomputeConstants in libsgp4, but use initial inclination
	if math.Abs(elem.cosio+1.0) > 1.5e-12 {
		elem.xlcof = 0.125 * a3ovk2 * elem.sinio * (3.0 + 5.0*elem.cosio) / (1.0 + elem.cosio)
	} else {
		elem.xlcof = 0.125 * a3ovk2 * elem.sinio * (3.0 + 5.0*elem.cosio) / 1.5e-12
	}
	elem.aycof = 0.25 * a3ovk2 * elem.sinio

	// NearSpace specific constants (assuming not deep space for ISS)
	// From SGP4::Initialise when use_deep_space_ is false
	var c3 float64
	if elem.ecc > 1.0e-4 {
		c3 = coef * tsi * a3ovk2 * elem.n * ae * elem.sinio / elem.ecc // ae=1
	}
	elem.c5 = 2.0 * coef1 * elem.a * betao2 * (1.0 + 2.75*(etasq+eeta) + eeta*etasq)
	elem.omgcof = elem.bstar * c3 * math.Cos(elem.omega)

	elem.xmcof = 0.0
	if elem.ecc > 1.0e-4 {
		elem.xmcof = -2.0 / 3.0 * coef * elem.bstar * ae / eeta // ae=1
	}
	elem.delmo = math.Pow(1.0+elem.eta*math.Cos(elem.m), 3.0)
	elem.sinmo = math.Sin(elem.m)

	if !elem.isSimpleModel {
		c1sq := elem.c1 * elem.c1
		elem.d2 = 4.0 * elem.a * tsi * c1sq
		dtemp := elem.d2 * tsi * elem.c1 / 3.0 // Renamed to avoid conflict
		elem.d3 = (17.0*elem.a + s4) * dtemp
		elem.d4 = 0.5 * dtemp * elem.a * tsi * (221.0*elem.a + 31.0*s4) * elem.c1
		elem.t3cof = elem.d2 + 2.0*c1sq
		elem.t4cof = 0.25 * (3.0*elem.d3 + elem.c1*(12.0*elem.d2+10.0*c1sq))
		elem.t5cof = 0.2 * (3.0*elem.d4 + 12.0*elem.c1*elem.d3 + 6.0*elem.d2*elem.d2 +
			15.0*c1sq*(2.0*elem.d2+c1sq))
	}
	return elem, nil
}
