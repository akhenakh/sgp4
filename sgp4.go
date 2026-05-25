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
	lat = AcTan(z, r)

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
		lat = AcTan(z+currentRe*c_iter*e2*sinLat, r)
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
		ecose = axn*cosepw + ayn_lpp*sinepw
		esine = axn*sinepw - ayn_lpp*cosepw

		f_kepler := capu - epw + esine
		if math.Abs(f_kepler) < 1.0e-12 {
			break
		}

		fdot_kepler := 1.0 - ecose
		delta_epw := f_kepler / fdot_kepler

		if i == 0 { // First iteration, apply cap
			if delta_epw > max_newton_raphson {
				delta_epw = max_newton_raphson
			} else if delta_epw < -max_newton_raphson {
				delta_epw = -max_newton_raphson
			}
		} else {
			// Second-order Newton-Raphson correction (matches libsgp4)
			delta_epw = f_kepler / (fdot_kepler + 0.5*esine*delta_epw)
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
		DateTime: tle.EpochTime().Add(time.Duration(tsince * float64(time.Minute))),
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

	// Create a wrapper function to get elevation for binary search
	getElevationAtTime := func(t time.Time) (float64, error) {
		tsince := t.Sub(tle.EpochTime()).Minutes()
		eciState, err := tle.FindPosition(tsince)
		if err != nil {
			return 0.0, err
		}
		sv := &StateVector{
			X: eciState.Position.X, Y: eciState.Position.Y, Z: eciState.Position.Z,
			VX: eciState.Velocity.X, VY: eciState.Velocity.Y, VZ: eciState.Velocity.Z,
		}
		observation, err := sv.GetLookAngle(observer, t)
		if err != nil {
			return 0.0, err
		}
		return observation.LookAngles.Elevation, nil
	}

	// findCrossingPoint finds the exact second the satellite crosses the horizon (0° elevation)
	// between time1 and time2. `findingAOS` should be true for AOS, false for LOS.
	findCrossingPoint := func(time1, time2 time.Time, findingAOS bool) (time.Time, error) {
		running := true
		cnt := 0
		middleTime := time1

		// Binary search for crossing point
		for running && cnt < 16 { // Limit iterations
			cnt++
			duration := time2.Sub(time1)
			middleTime = time1.Add(duration / 2)

			elevation, err := getElevationAtTime(middleTime)
			if err != nil {
				return time1, err // Return start time on error
			}

			if elevation > 0.0 {
				// Satellite is above horizon
				if findingAOS {
					time2 = middleTime
				} else {
					time1 = middleTime
				}
			} else {
				// Satellite is below horizon
				if findingAOS {
					time1 = middleTime
				} else {
					time2 = middleTime
				}
			}

			// Check if times are within 1 second
			if time2.Sub(time1).Seconds() < 1.0 {
				running = false
				// Truncate to whole seconds
				middleTime = time.Date(
					middleTime.Year(), middleTime.Month(), middleTime.Day(),
					middleTime.Hour(), middleTime.Minute(), middleTime.Second(),
					0, middleTime.Location(),
				)
				// Step back/forward 1 second into the pass
				if findingAOS {
					middleTime = middleTime.Add(time.Second)
				} else {
					middleTime = middleTime.Add(-time.Second)
				}
			}
		}

		// Final refinement: step back/forward until just below horizon
		running = true
		cnt = 0
		for running && cnt < 6 {
			cnt++
			elevation, err := getElevationAtTime(middleTime)
			if err != nil {
				break
			}

			if elevation > 0.0 {
				if findingAOS {
					middleTime = middleTime.Add(-time.Second)
				} else {
					middleTime = middleTime.Add(time.Second)
				}
			} else {
				running = false
			}
		}

		return middleTime, nil
	}

	// findMaxElevation finds the maximum elevation between aos and los times.
	findMaxElevation := func(aos, los time.Time) (maxEl float64, maxElTime time.Time, err error) {
		maxEl = -999999.9
		maxElTime = aos

		timeStep := los.Sub(aos) / 9 // Initial coarse step
		currentTime := aos

		for {
			running := true
			newMaxEl := -999999.9
			newMaxElTime := aos

			// Scan with current step size
			for running && currentTime.Before(los) {
				elevation, err := getElevationAtTime(currentTime)
				if err != nil {
					return maxEl, maxElTime, err
				}

				if elevation > newMaxEl {
					newMaxEl = elevation
					newMaxElTime = currentTime
					// Move to next step
					nextTime := currentTime.Add(timeStep)
					if nextTime.After(los) {
						nextTime = los
					}
					currentTime = nextTime
				} else {
					// Elevation is decreasing, stop this scan
					running = false
				}
			}

			// If max elevation didn't improve, we're done
			if newMaxEl <= maxEl {
				break
			}

			// Update best guess
			maxEl = newMaxEl
			maxElTime = newMaxElTime

			// Refine search window and step size
			time1 := currentTime.Add(-2 * timeStep) // Start 2 steps back
			if time1.Before(aos) {
				time1 = aos
			}
			time2 := currentTime // End at current time
			if time2.After(los) {
				time2 = los
			}

			currentTime = time1
			timeStep = time2.Sub(time1) / 9

			// Stop if step is less than 1 second
			if timeStep < time.Second {
				break
			}
		}

		return maxEl, maxElTime, nil
	}

	// Main pass generation loop
	currentTime := start
	previousTime := start
	foundAOS := false

	for currentTime.Before(stop) {
		endOfPass := false

		// Get current elevation
		elevation, err := getElevationAtTime(currentTime)
		if err != nil {
			// Skip this step on propagation error
			previousTime = currentTime
			currentTime = currentTime.Add(time.Duration(stepSeconds) * time.Second)
			if currentTime.After(stop) {
				currentTime = stop
			}
			continue
		}

		if !foundAOS && elevation > 0.0 {
			// Satellite has risen above horizon
			var aosTime time.Time
			if currentTime.Equal(start) {
				// Already above horizon at start
				aosTime = start
			} else {
				// Find exact AOS between previousTime and currentTime
				aos, err := findCrossingPoint(previousTime, currentTime, true)
				if err != nil {
					// Fallback to current time if error
					aosTime = currentTime
				} else {
					aosTime = aos
				}
			}
			foundAOS = true

			// Initialize new pass
			currentPass = &PassDetails{
				AOS: aosTime,
			}

			// Get AOS observation details
			tsinceAOS := aosTime.Sub(tle.EpochTime()).Minutes()
			eciAOS, _ := tle.FindPosition(tsinceAOS) // Ignore error, we just got elevation
			svAOS := &StateVector{X: eciAOS.Position.X, Y: eciAOS.Position.Y, Z: eciAOS.Position.Z}
			aosObs, _ := svAOS.GetLookAngle(observer, aosTime)
			currentPass.AOSAzimuth = aosObs.LookAngles.Azimuth
			currentPass.AOSObservation = *aosObs
		} else if foundAOS && elevation < 0.0 {
			// Satellite has set below horizon
			foundAOS = false
			endOfPass = true

			// Find exact LOS
			los, err := findCrossingPoint(previousTime, currentTime, false)
			if err != nil {
				los = currentTime // Fallback
			}

			// Finalize the pass
			currentPass.LOS = los

			// Get LOS observation details
			tsinceLOS := los.Sub(tle.EpochTime()).Minutes()
			eciLOS, _ := tle.FindPosition(tsinceLOS)
			svLOS := &StateVector{X: eciLOS.Position.X, Y: eciLOS.Position.Y, Z: eciLOS.Position.Z}
			losObs, _ := svLOS.GetLookAngle(observer, los)
			currentPass.LOSAzimuth = losObs.LookAngles.Azimuth
			currentPass.LOSObservation = *losObs

			// Find max elevation within the pass
			maxEl, maxElTime, err := findMaxElevation(currentPass.AOS, currentPass.LOS)
			if err == nil {
				currentPass.MaxElevation = maxEl
				currentPass.MaxElevationTime = maxElTime

				// Get observation at max elevation
				tsinceMax := maxElTime.Sub(tle.EpochTime()).Minutes()
				eciMax, _ := tle.FindPosition(tsinceMax)
				svMax := &StateVector{X: eciMax.Position.X, Y: eciMax.Position.Y, Z: eciMax.Position.Z}
				maxElObs, _ := svMax.GetLookAngle(observer, maxElTime)
				currentPass.MaxElevationAz = maxElObs.LookAngles.Azimuth
				currentPass.MaxElObservation = *maxElObs
			}

			currentPass.Duration = currentPass.LOS.Sub(currentPass.AOS)
			passes = append(passes, *currentPass)
			currentPass = nil
		}

		// Update times
		previousTime = currentTime
		if endOfPass {
			// Jump ahead after a pass to avoid unnecessary calculations
			currentTime = currentTime.Add(30 * time.Minute)
		} else {
			currentTime = currentTime.Add(time.Duration(stepSeconds) * time.Second)
		}

		if currentTime.After(stop) {
			currentTime = stop
		}
	}

	// Handle case where satellite is still visible at end time
	if foundAOS && currentPass != nil {
		currentPass.LOS = stop

		// Get LOS observation details (at stop time)
		tsinceLOS := stop.Sub(tle.EpochTime()).Minutes()
		eciLOS, _ := tle.FindPosition(tsinceLOS)
		svLOS := &StateVector{X: eciLOS.Position.X, Y: eciLOS.Position.Y, Z: eciLOS.Position.Z}
		losObs, _ := svLOS.GetLookAngle(observer, stop)
		currentPass.LOSAzimuth = losObs.LookAngles.Azimuth
		currentPass.LOSObservation = *losObs

		// Find max elevation within the pass
		maxEl, maxElTime, err := findMaxElevation(currentPass.AOS, currentPass.LOS)
		if err == nil {
			currentPass.MaxElevation = maxEl
			currentPass.MaxElevationTime = maxElTime

			// Get observation at max elevation
			tsinceMax := maxElTime.Sub(tle.EpochTime()).Minutes()
			eciMax, _ := tle.FindPosition(tsinceMax)
			svMax := &StateVector{X: eciMax.Position.X, Y: eciMax.Position.Y, Z: eciMax.Position.Z}
			maxElObs, _ := svMax.GetLookAngle(observer, maxElTime)
			currentPass.MaxElevationAz = maxElObs.LookAngles.Azimuth
			currentPass.MaxElObservation = *maxElObs
		}

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

// AcTan calculates the arctangent in the correct quadrant based on sin and cos values.
// It replicates the behavior of the C++ libsgp4 Util::AcTan function.
func AcTan(sinx, cosx float64) float64 {
	// Handle the case where cosx is effectively zero to avoid division by zero
	// and determine the correct quadrant based on the sign of sinx.
	if math.Abs(cosx) < 1e-14 { // Using a small epsilon, similar to other checks
		if sinx > 0.0 {
			return math.Pi / 2.0
		} else {
			// Note: C++ returns 3*PI/2. This is equivalent to -PI/2 in terms of angle,
			// but let's stick exactly to the C++ logic for maximum fidelity.
			// If the downstream code expects [-PI, PI] range, you might consider
			// returning -math.Pi / 2.0 instead, but C++ uses [0, 2*PI).
			return 3.0 * math.Pi / 2.0
		}
	} else {
		// If cosx is positive, the angle is in quadrants I or IV.
		// math.Atan will return the correct angle in [-PI/2, PI/2].
		if cosx > 0.0 {
			return math.Atan(sinx / cosx)
		} else {
			// If cosx is negative, the angle is in quadrants II or III.
			// math.Atan will return an angle in [-PI/2, PI/2], but the actual
			// angle is in [PI/2, 3*PI/2]. Adding PI corrects the quadrant.
			// This also correctly handles the case where sinx is zero and cosx is negative
			// (angle is PI).
			return math.Pi + math.Atan(sinx/cosx) // Adding PI shifts to correct quadrant
		}
	}
}
