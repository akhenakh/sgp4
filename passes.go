package spg4

import (
	"fmt"
	"math"
	"time"
)

// Pass represents a single pass of a satellite over a ground station
type Pass struct {
	AOS       time.Time          // Time of Acquisition of Signal (when satellite rises)
	TCA       time.Time          // Time of Closest Approach (max elevation)
	LOS       time.Time          // Time of Loss of Signal (when satellite sets)
	MaxEl     float64            // Maximum elevation during pass (degrees)
	MaxAzEl   *TopocentricCoords // Az/El at max elevation point
	StartAzEl *TopocentricCoords // Az/El at start of pass
	EndAzEl   *TopocentricCoords // Az/El at end of pass
}

// minutesSince calculates minutes elapsed between two time instances
func minutesSince(start, end time.Time) float64 {
	return end.Sub(start).Minutes()
}

// findCrossing uses binary search to find the exact time when elevation crosses the
// specified threshold. rising=true for AOS, false for LOS
func findCrossing(tle *TLE, loc *Location, t1, t2 float64, threshold float64, rising bool) float64 {
	const tolerance = 0.001 // Tolerance in degrees
	for i := 0; i < 50 && math.Abs(t2-t1) > 0.0001; i++ {
		tmid := (t1 + t2) / 2
		pos, err := tle.GetPosition(tmid)
		if err != nil {
			return tmid // Return mid point in case of error
		}

		el := pos.GetLookAngle(loc).Elevation
		if rising {
			if el < threshold {
				t1 = tmid
			} else {
				t2 = tmid
			}
		} else {
			if el > threshold {
				t1 = tmid
			} else {
				t2 = tmid
			}
		}

		if math.Abs(el-threshold) < tolerance {
			return tmid
		}
	}
	return (t1 + t2) / 2
}

// findMaxElevation finds the time of maximum elevation during a pass
// and returns the maximum elevation and its time
func findMaxElevation(tle *TLE, loc *Location, start, end float64) (maxEl float64, maxTime float64) {
	const steps = 100
	dt := (end - start) / steps

	maxEl = -90.0
	maxTime = start

	// First pass: coarse search
	for t := start; t <= end; t += dt {
		pos, err := tle.GetPosition(t)
		if err != nil {
			continue
		}
		el := pos.GetLookAngle(loc).Elevation
		if el > maxEl {
			maxEl = el
			maxTime = t
		}
	}

	// Second pass: fine search around maximum
	if maxTime > start && maxTime < end {
		fineStart := maxTime - dt
		fineEnd := maxTime + dt
		fineDt := dt / 10

		for t := fineStart; t <= fineEnd; t += fineDt {
			pos, err := tle.GetPosition(t)
			if err != nil {
				continue
			}
			el := pos.GetLookAngle(loc).Elevation
			if el > maxEl {
				maxEl = el
				maxTime = t
			}
		}
	}

	return maxEl, maxTime
}

// getAzElAtTime calculates azimuth and elevation at a specific time
func getAzElAtTime(tle *TLE, loc *Location, t float64) *TopocentricCoords {
	pos, err := tle.GetPosition(t)
	if err != nil {
		return &TopocentricCoords{} // Return empty coords in case of error
	}
	return pos.GetLookAngle(loc)
}

// FindPasses calculates all visible passes of the satellite over the given location
// within the specified time window. A pass is considered visible if it reaches an
// elevation above minElevation degrees.
func (tle *TLE) FindPasses(loc *Location, start time.Time, end time.Time, minElevation float64) ([]*Pass, error) {
	const (
		stepSize   = 1.0  // Step size in minutes for coarse search
		fineStep   = 0.1  // Step size in minutes for fine search
		maxPasses  = 100  // Safety limit for number of passes
		marginMins = 10.0 // Time margin in minutes for finding exact AOS/LOS
	)

	var passes []*Pass
	var inPass bool
	var currentPass *Pass
	var lastEl float64
	var aosTime float64 // Store AOS time for the current pass

	// Convert times to minutes since epoch
	epoch := tle.EpochTime()
	startMins := minutesSince(epoch, start)
	endMins := minutesSince(epoch, end)

	// Coarse search for passes
	for t := startMins; t <= endMins; t += stepSize {
		pos, err := tle.GetPosition(t)
		if err != nil {
			return nil, fmt.Errorf("position calculation failed at t=%f: %w", t, err)
		}

		azEl := pos.GetLookAngle(loc)
		el := azEl.Elevation

		// Detect pass start (rising above horizon + minElevation)
		if !inPass && lastEl < minElevation && el >= minElevation {
			currentPass = &Pass{}
			// Fine search for exact AOS
			aosTime = findCrossing(tle, loc, t-marginMins, t+marginMins, minElevation, true)
			currentPass.AOS = epoch.Add(time.Duration(aosTime * float64(time.Minute)))
			currentPass.StartAzEl = getAzElAtTime(tle, loc, aosTime)
			inPass = true
		}

		// Detect pass end (setting below minElevation)
		if inPass && lastEl >= minElevation && el < minElevation {
			// Fine search for exact LOS
			losTime := findCrossing(tle, loc, t-marginMins, t+marginMins, minElevation, false)
			currentPass.LOS = epoch.Add(time.Duration(losTime * float64(time.Minute)))
			currentPass.EndAzEl = getAzElAtTime(tle, loc, losTime)

			// Find maximum elevation during pass
			maxEl, tca := findMaxElevation(tle, loc, aosTime, losTime)
			currentPass.MaxEl = maxEl
			currentPass.TCA = epoch.Add(time.Duration(tca * float64(time.Minute)))
			currentPass.MaxAzEl = getAzElAtTime(tle, loc, tca)

			passes = append(passes, currentPass)
			inPass = false

			if len(passes) >= maxPasses {
				break
			}
		}

		lastEl = el
	}

	return passes, nil
}
