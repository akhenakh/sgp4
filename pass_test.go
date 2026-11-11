package sgp4

import (
	"errors"
	"math"
	"testing"
	"time"
)

const passTestTLE = `ISS (ZARYA)
1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994
2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.49587957510533`

// firstTallPass returns the first pass in a 24 h window whose maximum elevation
// is high enough for a meaningful range-rate sign test.
func firstTallPass(t *testing.T) (*TLE, PassDetails) {
	t.Helper()
	tle, err := ParseTLE(passTestTLE)
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	start := tle.EpochTime()
	passes, err := tle.GeneratePasses(48.85, 2.35, 50, start, start.Add(24*time.Hour), 60)
	if err != nil {
		t.Fatalf("GeneratePasses: %v", err)
	}
	for _, p := range passes {
		if p.MaxElevation > 20 {
			return tle, p
		}
	}
	t.Fatalf("no pass above 20 degrees in the window (%d passes)", len(passes))
	return nil, PassDetails{}
}

func TestGeneratePassesPopulatesDataPoints(t *testing.T) {
	_, p := firstTallPass(t)
	if len(p.DataPoints) < 2 {
		t.Fatalf("pass has %d data points, want at least 2", len(p.DataPoints))
	}
	if !p.DataPoints[0].Timestamp.Equal(p.AOS) {
		t.Errorf("first data point at %v, want AOS %v", p.DataPoints[0].Timestamp, p.AOS)
	}
	if last := p.DataPoints[len(p.DataPoints)-1]; !last.Timestamp.Equal(p.LOS) {
		t.Errorf("last data point at %v, want LOS %v", last.Timestamp, p.LOS)
	}
	for i := 1; i < len(p.DataPoints); i++ {
		if !p.DataPoints[i].Timestamp.After(p.DataPoints[i-1].Timestamp) {
			t.Fatalf("data points are not in time order at index %d", i)
		}
	}
	if svg := p.GeneratePassPolarSVG(); len(svg) == 0 {
		t.Error("GeneratePassPolarSVG returned an empty plot")
	}
}

// TestRangeRateIsPhysical guards the units of the observer velocity: with the
// rotations-per-day ratio used as an angular velocity, range rate came out in
// the thousands of km/s.
func TestRangeRateIsPhysical(t *testing.T) {
	_, p := firstTallPass(t)

	// A near-Earth satellite moves at under 8 km/s and the observer at under
	// 0.5 km/s, so no range rate can exceed their sum.
	const maxRate = 8.5
	for i, dp := range p.DataPoints {
		if math.IsNaN(dp.RangeRate) || math.IsInf(dp.RangeRate, 0) {
			t.Fatalf("data point %d has a non-finite range rate: %v", i, dp.RangeRate)
		}
		if math.Abs(dp.RangeRate) > maxRate {
			t.Errorf("data point %d range rate %.3f km/s exceeds %.1f km/s", i, dp.RangeRate, maxRate)
		}
	}

	// Approaching at AOS, receding at LOS. This also proves the satellite
	// velocity was carried into the pass observations.
	if p.AOSObservation.LookAngles.RangeRate >= 0 {
		t.Errorf("AOS range rate %.3f, want negative (approaching)", p.AOSObservation.LookAngles.RangeRate)
	}
	if p.LOSObservation.LookAngles.RangeRate <= 0 {
		t.Errorf("LOS range rate %.3f, want positive (receding)", p.LOSObservation.LookAngles.RangeRate)
	}
	if math.Abs(p.MaxElObservation.LookAngles.RangeRate) > 1 {
		t.Errorf("range rate at maximum elevation %.3f, want near zero", p.MaxElObservation.LookAngles.RangeRate)
	}
}

func TestGeneratePassesRejectsInvalidObserver(t *testing.T) {
	tle, err := ParseTLE(passTestTLE)
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}
	start := tle.EpochTime()
	stop := start.Add(time.Hour)

	tests := []struct {
		name string
		lat  float64
		lng  float64
		alt  float64
		want error
	}{
		{"latitude above 90", 91, 0, 0, ErrInvalidLocationLatitude},
		{"latitude below -90", -91, 0, 0, ErrInvalidLocationLatitude},
		{"NaN latitude", math.NaN(), 0, 0, ErrInvalidLocationLatitude},
		{"NaN longitude", 0, math.NaN(), 0, ErrInvalidLocationLongitude},
		{"NaN altitude", 0, 0, math.NaN(), ErrInvalidLocationAltitude},
		{"infinite altitude", 0, 0, math.Inf(1), ErrInvalidLocationAltitude},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			passes, err := tle.GeneratePasses(tt.lat, tt.lng, tt.alt, start, stop, 60)
			if !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
			if passes != nil {
				t.Errorf("returned %d passes alongside the error", len(passes))
			}
		})
	}
}

func TestGeneratePassPolarSVGStable(t *testing.T) {
	_, p := firstTallPass(t)
	first := p.GeneratePassPolarSVG()
	second := p.GeneratePassPolarSVG()
	if first != second {
		t.Error("GeneratePassPolarSVG output changes between calls")
	}
}
