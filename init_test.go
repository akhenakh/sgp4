package sgp4

import (
	"errors"
	"math"
	"testing"
)

// TestInitializeRejectsDegenerateElements covers element sets that are
// constructible (hand-built or converted from OMM) but cannot be propagated:
// they used to return NaN positions with a nil error.
func TestInitializeRejectsDegenerateElements(t *testing.T) {
	base, err := ParseTLE(passTestTLE)
	if err != nil {
		t.Fatalf("ParseTLE: %v", err)
	}

	tests := []struct {
		name   string
		mutate func(*TLE)
		reason SGP4ModelLimitsErrorReason
	}{
		{"zero mean motion", func(x *TLE) { x.MeanMotion = 0 }, ReasonInvalidMeanMotion},
		{"negative mean motion", func(x *TLE) { x.MeanMotion = -1 }, ReasonInvalidMeanMotion},
		{"NaN mean motion", func(x *TLE) { x.MeanMotion = math.NaN() }, ReasonInvalidMeanMotion},
		{"infinite mean motion", func(x *TLE) { x.MeanMotion = math.Inf(1) }, ReasonInvalidMeanMotion},
		{"negative eccentricity", func(x *TLE) { x.Eccentricity = -0.1 }, ReasonInvalidEccentricity},
		{"eccentricity of 1", func(x *TLE) { x.Eccentricity = 1 }, ReasonInvalidEccentricity},
		{"eccentricity above 1", func(x *TLE) { x.Eccentricity = 1.5 }, ReasonInvalidEccentricity},
		{"NaN eccentricity", func(x *TLE) { x.Eccentricity = math.NaN() }, ReasonInvalidEccentricity},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tle := *base
			tt.mutate(&tle)

			eci, err := tle.FindPosition(0)
			if err == nil {
				t.Fatalf("FindPosition returned %v with no error", eci.Position)
			}
			var limits *SGP4ModelLimitsError
			if !errors.As(err, &limits) {
				t.Fatalf("error = %v, want SGP4ModelLimitsError", err)
			}
			if limits.Reason != tt.reason {
				t.Errorf("reason = %q, want %q", limits.Reason, tt.reason)
			}
		})
	}
}
