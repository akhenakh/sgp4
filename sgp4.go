package sgp4

import (
	"fmt"
	"math"
	"time"
)

const MinElevationForPass = 10.0 // degrees

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
	currentRe := reSGP4 // Use SGP4-aligned constants
	currentF := fSGP4
	e2 := currentF * (2.0 - currentF)

	gmst := eci.GreenwichSiderealTime()
	x := eci.Position.X
	y := eci.Position.Y
	z := eci.Position.Z

	lon = math.Atan2(y, x) - gmst
	lon = wrapLongitude(lon)
	r := math.Sqrt(x*x + y*y)
	lat = math.Atan2(z, r)

	const maxIter = 10
	const tol = 1e-10
	var oldLat float64
	var c_iter float64

	for range maxIter {
		oldLat = lat
		sinLat := math.Sin(lat)
		if math.Abs(1.0-e2*sinLat*sinLat) < 1e-14 { // Avoid division by zero or sqrt of negative
			c_iter = 1.0 / math.Sqrt(1e-14)
		} else {
			c_iter = 1.0 / math.Sqrt(1.0-e2*sinLat*sinLat)
		}
		lat = math.Atan2(z+currentRe*c_iter*e2*sinLat, r)
		if math.Abs(lat-oldLat) < tol {
			break
		}
	}

	sinLat := math.Sin(lat)
	cosLat := math.Cos(lat)
	var N_val float64
	if math.Abs(1.0-e2*sinLat*sinLat) < 1e-14 {
		N_val = currentRe * (1.0 / math.Sqrt(1e-14))
	} else {
		N_val = currentRe * (1.0 / math.Sqrt(1.0-e2*sinLat*sinLat))
	}

	if math.Abs(cosLat) < 1e-10 {
		alt = math.Abs(z) - currentRe*math.Sqrt(1.0-e2)
	} else {
		alt = r/cosLat - N_val
	}

	lat = lat * rad2deg
	lon = lon * rad2deg
	return lat, lon, alt
}

// FindPosition propagates the TLE to the given time offset (tsince) in minutes
// FindPosition propagates the TLE to the given time offset (tsince) in minutes
func (tle *TLE) FindPosition(tsince float64) (Eci, error) {
	elems, err := tle.Initialize()
	if err != nil {
		return Eci{}, fmt.Errorf("SGP4 propagation error during initialization: %w", err)
	}

	// Retrieve initialized SGP4 elements and constants
	// These are mean elements at epoch, some are already perturbed by J2 (a, n)
	// Others are TLE raw values (ecc, m, omega, raan, incl)
	// And derived SGP4 constants (c1, c4, xmdot, etc.)

	// Secular effects (similar to SGP4::FindPositionSGP4 before CalculateFinalPositionVelocity)
	xmdf := elems.m + elems.xmdot*tsince
	omgadf := elems.omega + elems.omgdot*tsince
	xnoddf := elems.raan + elems.xnodot*tsince

	omega := omgadf // current argument of perigee
	xmp := xmdf     // current mean anomaly

	tsq := tsince * tsince
	xnode := xnoddf + elems.xnodcf*tsq // current RAAN
	tempa := 1.0 - elems.c1*tsince
	tempe := elems.bstar * elems.c4 * tsince
	templ := elems.t2cof * tsq

	if !elems.isSimpleModel {
		delomg := elems.omgcof * tsince
		delm_term := 0.0
		if elems.eta != 0.0 { // Avoid division by zero if eta is zero (circular orbit or specific init issue)
			delm_term = elems.xmcof * (math.Pow(1.0+elems.eta*math.Cos(xmdf), 3.0) - elems.delmo)
		}

		temp := delomg + delm_term
		xmp += temp
		omega -= temp

		tcube := tsq * tsince
		tfour := tsince * tcube
		tempa = tempa - elems.d2*tsq - elems.d3*tcube - elems.d4*tfour
		tempe += elems.bstar * elems.c5 * (math.Sin(xmp) - elems.sinmo)
		templ += elems.t3cof*tcube + tfour*(elems.t4cof+tsince*elems.t5cof)
	}

	a := elems.a * tempa * tempa              // current semi-major axis (ER)
	e := elems.ecc - tempe                    // current eccentricity
	xl := xmp + omega + xnode + elems.n*templ // current mean longitude (M + omega + Omega)

	// Ensure eccentricity is within sane bounds
	if e <= -0.001 { // Check current eccentricity 'e'
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonEccentricityTooLow, Value: e}
	} else if e < 1.0e-6 {
		e = 1.0e-6
	} else if e > (1.0 - 1.0e-6) { // Near parabolic
		// SGP4 usually caps this, but if we wanted to error:
		// return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonEccentricityTooHigh, Value: e}
		e = 1.0 - 1.0e-6
	}

	// Call CalculateFinalPositionVelocity equivalent
	// This is where short-period perturbations are applied.
	// Constants cosio, sinio, x3thm1, x1mth2, xlcof, aycof are from the *initial* inclination (elems.cosio etc)
	// x7thm1 was calculated during Initialize as well
	beta2 := 1.0 - e*e // current beta_sq, based on current 'e'
	if beta2 < 0.0 {   // Should not happen if e is capped correctly to < 1
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonBeta2Negative, Value: beta2, Message: fmt.Sprintf("based on current eccentricity e=%.6e", e)}
	}
	xn := xke / math.Pow(a, 1.5) // current mean motion (rad/min)

	// Long period periodics (LPP) affecting argument of latitude
	axn := e * math.Cos(omega) // LPP term for L' (e_k * cos(omega_k))
	temp11_lpp := 1.0 / (a * beta2)
	xll_lpp := temp11_lpp * elems.xlcof * axn
	aynl_lpp := temp11_lpp * elems.aycof
	xlt_lpp := xl + xll_lpp                 // L' = L + Lpp
	ayn_lpp := e*math.Sin(omega) + aynl_lpp // LPP term for L' (e_k * sin(omega_k))

	elsq := axn*axn + ayn_lpp*ayn_lpp // (e_k)^2, where e_k is the LPP-perturbed eccentricity vector magnitude
	if elsq >= 1.0 {
		// This can happen if perturbations drive eccentricity too high.
		// For robust code, might cap elsq or return error.
		// libsgp4 throws "Error: (elsq >= 1.0)"
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonPerturbedEccSqTooHigh, Value: elsq}
	}

	// Solve Kepler's equation for L' (eccentric longitude)
	// capu is M' = L' - Omega - omega (where Omega and omega are current nodek and omegak)
	// In SGP4, this is M_k = L' - Omega_k (argument of latitude from ascending node)
	capu := math.Mod(xlt_lpp-xnode, twoPi) // E_k_mean = L' - Omega_k
	epw := capu                            // Initial guess for eccentric anomaly (perturbed) E_k_ecc

	var sinepw, cosepw, ecose, esine float64
	max_newton_raphson := 1.25 * math.Abs(math.Sqrt(elsq)) // Cap on N-R step

	for i := range 10 {
		sinepw = math.Sin(epw)
		cosepw = math.Cos(epw)
		ecose = axn*cosepw + ayn_lpp*sinepw // e_k cos(E_k_ecc)
		esine = axn*sinepw - ayn_lpp*cosepw // e_k sin(E_k_ecc)
		f_kepler := capu - epw + esine      // Kepler's eq: M' - E_k_ecc + e_k sin(E_k_ecc) = 0
		if math.Abs(f_kepler) < 1.0e-12 {
			break
		}
		fdot_kepler := 1.0 - ecose // d(Kepler)/d(E_k_ecc)
		delta_epw := f_kepler / fdot_kepler
		if i == 0 { // First iteration, apply cap
			if delta_epw > max_newton_raphson {
				delta_epw = max_newton_raphson
			} else if delta_epw < -max_newton_raphson {
				delta_epw = -max_newton_raphson
			}
		}
		epw += delta_epw
	}

	// Short period preliminary quantities
	temp21_sp := max(1.0-elsq, 0.0) // Ensure non-negative for sqrt
	pl := a * temp21_sp             // semi-latus rectum p_k = a_k * (1 - (e_k)^2)
	if pl < 0.0 {                   // Should be caught by elsq >= 1.0 or a < 0 (if tempa becomes very negative)
		return Eci{}, &SGP4ModelLimitsError{Tsince: tsince, Reason: ReasonSemiLatusRectumNegative, Value: pl}
	}

	r_val := a * (1.0 - ecose) // distance from primary focus r_k = a_k * (1 - e_k cos(E_k_ecc))
	if r_val == 0.0 {
		r_val = 1e-9
	} // Avoid division by zero
	temp31_sp := 1.0 / r_val
	rdot_val := xke * math.Sqrt(a) * esine * temp31_sp // r_dot_k
	rfdot_val := xke * math.Sqrt(pl) * temp31_sp       // r_k * f_dot_k (f is true anomaly)

	temp32_sp := a * temp31_sp       // a_k / r_k
	betal_sp := math.Sqrt(temp21_sp) // beta_k = sqrt(1 - (e_k)^2)
	temp33_sp := 0.0
	if (1.0 + betal_sp) != 0.0 {
		temp33_sp = 1.0 / (1.0 + betal_sp)
	} else {
		// This case (betal_sp = -1) should not happen for real eccentricities
		temp33_sp = 1.0e12 // Avoid division by zero, effectively making next terms small
	}

	cosu_sp := temp32_sp * (cosepw - axn + ayn_lpp*esine*temp33_sp)
	sinu_sp := temp32_sp * (sinepw - ayn_lpp - axn*esine*temp33_sp)
	u_sp := math.Atan2(sinu_sp, cosu_sp)

	sin2u_sp := 2.0 * sinu_sp * cosu_sp
	cos2u_sp := 2.0*cosu_sp*cosu_sp - 1.0

	// Short period perturbations (SPP)
	// Constants x3thm1, x1mth2, cosio, sinio, x7thm1 are from initial elements
	temp41_spp := 0.0
	if pl != 0.0 {
		temp41_spp = 1.0 / pl
	} else {
		temp41_spp = 1.0e12
	}

	temp42_spp := ck2 * temp41_spp
	temp43_spp := temp42_spp * temp41_spp

	rk_spp := r_val*(1.0-1.5*temp43_spp*betal_sp*elems.x3thm1) + 0.5*temp42_spp*elems.x1mth2*cos2u_sp
	uk_spp := u_sp - 0.25*temp43_spp*elems.x7thm1*sin2u_sp
	xnodek_spp := xnode + 1.5*temp43_spp*elems.cosio*sin2u_sp
	xinck_spp := elems.incl + 1.5*temp43_spp*elems.cosio*elems.sinio*cos2u_sp
	rdotk_spp := rdot_val - xn*temp42_spp*elems.x1mth2*sin2u_sp
	rfdotk_spp := rfdot_val + xn*temp42_spp*(elems.x1mth2*cos2u_sp+1.5*elems.x3thm1)

	// Orientation vectors (using perturbed elements)
	sinuk := math.Sin(uk_spp)
	cosuk := math.Cos(uk_spp)
	sinik := math.Sin(xinck_spp)
	cosik := math.Cos(xinck_spp)
	sinnok := math.Sin(xnodek_spp)
	cosnok := math.Cos(xnodek_spp)

	xmx := -sinnok * cosik
	xmy := cosnok * cosik

	// Position in ECI (km)
	// ux, uy, uz are components of position unit vector in ECI
	ux := xmx*sinuk + cosnok*cosuk
	uy := xmy*sinuk + sinnok*cosuk
	uz := sinik * sinuk
	posX := rk_spp * ux * xkmper
	posY := rk_spp * uy * xkmper
	posZ := rk_spp * uz * xkmper

	// Velocity in ECI (km/s)
	// vx, vy, vz are components of (velocity / (r_k * f_dot_k)) unit vector, or similar
	vx_orient := xmx*cosuk - cosnok*sinuk
	vy_orient := xmy*cosuk - sinnok*sinuk
	vz_orient := sinik * cosuk

	// rdotk_spp is in ER/min, rfdotk_spp is in ER/min
	// Convert to km/s by multiplying by (xkmper / 60.0)
	vFactor := xkmper / 60.0
	velX := (rdotk_spp*ux + rfdotk_spp*vx_orient) * vFactor
	velY := (rdotk_spp*uy + rfdotk_spp*vy_orient) * vFactor
	velZ := (rdotk_spp*uz + rfdotk_spp*vz_orient) * vFactor

	// Check for decay (physical condition)
	if rk_spp < 1.0 {
		// Satellite has decayed. SGP4 docs say prediction is unreliable.
		return Eci{}, &SatelliteDecayedError{Tsince: tsince, Radius: rk_spp}
	}

	return Eci{
		DateTime: tle.EpochTime().Add(time.Duration(tsince) * time.Minute),
		Position: Vector{X: posX, Y: posY, Z: posZ},
		Velocity: Vector{X: velX, Y: velY, Z: velZ},
	}, nil
}

// GreenwichSiderealTime calculates the Greenwich Mean Sidereal Time
func (eci *Eci) GreenwichSiderealTime() float64 {
	jd := julianDateTime(eci.DateTime)
	t := (jd - 2451545.0) / 36525.0

	gmst_deg := 280.46061837 +
		360.98564736629*(jd-2451545.0) +
		0.000387933*t*t -
		t*t*t/38710000.0

	gmst_deg = math.Mod(gmst_deg, 360.0)
	if gmst_deg < 0 {
		gmst_deg += 360.0
	}
	return gmst_deg * deg2rad
}

// FindPositionAtTime propagates the TLE to a specific absolute time.
func (tle *TLE) FindPositionAtTime(t time.Time) (Eci, error) {
	// Calculate tsince (time since TLE epoch in minutes)
	tsince := t.Sub(tle.EpochTime()).Minutes()
	return tle.FindPosition(tsince)
}

// GeneratePasses predicts satellite passes over a ground station within a given time window.
// lat, lng are observer's geodetic latitude/longitude in degrees.
// alt is observer's altitude in meters above sea level.
// start, stop are the time window boundaries.
// stepSeconds is the time step for propagation in seconds.
func (tle *TLE) GeneratePasses(obsLat, obsLng, obsAltMeters float64, start, stop time.Time, stepSeconds int) ([]PassDetails, error) {
	if start.After(stop) {
		return nil, fmt.Errorf("start time must be before stop time")
	}
	if stepSeconds <= 0 {
		return nil, fmt.Errorf("stepSeconds must be positive")
	}

	observer := &Location{
		Latitude:  obsLat,
		Longitude: obsLng,
		Altitude:  obsAltMeters,
	}

	var passes []PassDetails
	var currentPass *PassDetails
	var prevObservation *Observation

	stepDuration := time.Duration(stepSeconds) * time.Second
	currentTime := start

	for !currentTime.After(stop) {
		tsince := currentTime.Sub(tle.EpochTime()).Minutes()

		eciState, err := tle.FindPosition(tsince)
		if err != nil {
			// If a propagation error occurs (e.g. decay or model breakdown),
			// end any current pass and skip this step.
			// fmt.Printf("Warning: Could not propagate for time %v: %v\n", currentTime, err)
			if currentPass != nil { // Finalize pass if it was ongoing and propagation fails
				// Use previous observation for LOS if possible
				if prevObservation != nil {
					currentPass.LOS = prevObservation.SatellitePos.Timestamp
					currentPass.LOSAzimuth = prevObservation.LookAngles.Azimuth
					currentPass.LOSObservation = *prevObservation
				} else { // Fallback if no previous point (pass started and ended at this error point)
					currentPass.LOS = currentTime // Approximate LOS time
					// LOS Azimuth and Observation might be inaccurate or unavailable here.
					// Use AOS details if nothing else.
					currentPass.LOSAzimuth = currentPass.AOSAzimuth
					currentPass.LOSObservation = currentPass.AOSObservation
				}
				currentPass.Duration = currentPass.LOS.Sub(currentPass.AOS)
				passes = append(passes, *currentPass)
				currentPass = nil
			}
			currentTime = currentTime.Add(stepDuration)
			prevObservation = nil // Reset prevObservation due to error
			continue
		}

		sv := &StateVector{
			X: eciState.Position.X, Y: eciState.Position.Y, Z: eciState.Position.Z,
			VX: eciState.Velocity.X, VY: eciState.Velocity.Y, VZ: eciState.Velocity.Z,
		}

		currentObservation, err := sv.GetLookAngle(observer, currentTime)
		if err != nil {
			// fmt.Printf("Warning: Could not get look angle for time %v: %v\n", currentTime, err)
			currentTime = currentTime.Add(stepDuration)
			prevObservation = nil
			continue
		}

		// Pass detection logic
		isCurrentlyVisible := currentObservation.LookAngles.Elevation >= MinElevationForPass

		if isCurrentlyVisible {
			if currentPass == nil { // Start of a new pass (AOS)
				currentPass = &PassDetails{
					AOS:              currentTime,
					AOSAzimuth:       currentObservation.LookAngles.Azimuth,
					MaxElevation:     currentObservation.LookAngles.Elevation, // Initialize with first point
					MaxElevationAz:   currentObservation.LookAngles.Azimuth,
					MaxElevationTime: currentTime,
					AOSObservation:   *currentObservation,
					MaxElObservation: *currentObservation, // Set initial MaxElObservation
					DataPoints:       []PassDataPoint{},   // Initialize the slice
				}
			}

			// Add data point to the current pass
			currentPass.DataPoints = append(currentPass.DataPoints, PassDataPoint{
				Timestamp: currentTime,
				Azimuth:   currentObservation.LookAngles.Azimuth,
				Elevation: currentObservation.LookAngles.Elevation,
				Range:     currentObservation.LookAngles.Range,
				RangeRate: currentObservation.LookAngles.RangeRate,
			})

			// Update max elevation for the current pass
			if currentObservation.LookAngles.Elevation > currentPass.MaxElevation {
				currentPass.MaxElevation = currentObservation.LookAngles.Elevation
				currentPass.MaxElevationAz = currentObservation.LookAngles.Azimuth
				currentPass.MaxElevationTime = currentTime
				currentPass.MaxElObservation = *currentObservation
			}
		} else { // Not currently visible
			if currentPass != nil { // End of the current pass (LOS)
				currentPass.LOS = currentTime
				currentPass.LOSAzimuth = currentObservation.LookAngles.Azimuth

				// Take LOS details from the last point that was visible (prevObservation) if available
				// or the current one if it just dropped below the horizon
				if prevObservation != nil && prevObservation.LookAngles.Elevation >= MinElevationForPass {
					currentPass.LOS = prevObservation.SatellitePos.Timestamp // LOS is end of previous visible step
					currentPass.LOSAzimuth = prevObservation.LookAngles.Azimuth
					currentPass.LOSObservation = *prevObservation
				} else {
					// Fallback: If somehow prevObservation isn't above horizon, use current (less ideal for LOS)
					currentPass.LOSObservation = *currentObservation
				}

				currentPass.Duration = currentPass.LOS.Sub(currentPass.AOS)
				passes = append(passes, *currentPass)
				currentPass = nil // Reset for the next pass
			}
		}

		prevObservation = currentObservation
		currentTime = currentTime.Add(stepDuration)
	}

	// If the loop ends and a pass is still ongoing (satellite is visible at 'stop' time)
	if currentPass != nil {
		currentPass.LOS = stop      // Pass extends to the end of the window
		if prevObservation != nil { // Use the last valid observation for LOS details
			currentPass.LOSAzimuth = prevObservation.LookAngles.Azimuth
			currentPass.LOSObservation = *prevObservation
		}
		// No need for 'else if currentObservation != nil' here, as currentObservation is out of scope.
		// If prevObservation is nil, it implies pass started and finished within a single step at stop time,
		// and AOSObservation would hold the only point.

		currentPass.Duration = currentPass.LOS.Sub(currentPass.AOS)
		passes = append(passes, *currentPass)
	}

	return passes, nil
}

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

func wrapLongitude(lon float64) float64 {
	lon = math.Mod(lon, twoPi)
	if lon > math.Pi {
		lon -= twoPi
	} else if lon < -math.Pi {
		lon += twoPi
	}
	return lon
}
