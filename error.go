package sgp4

import (
	"fmt"
)

// SatelliteDecayedError is returned when the SGP4 model predicts the satellite has decayed.
type SatelliteDecayedError struct {
	Tsince float64 // Time since epoch in minutes when decay was detected
	Radius float64 // Final orbital radius in Earth Radii that triggered decay
}

// Error returns the error message for SatelliteDecayedError.
func (e *SatelliteDecayedError) Error() string {
	return fmt.Sprintf("SGP4: satellite has decayed (at tsince %.2f min, final orbital radius %.4f < 1.0 Earth radii)", e.Tsince, e.Radius)
}

// SGP4ModelLimitsErrorReason defines the specific reason for the model limit violation.
type SGP4ModelLimitsErrorReason string

const (
	ReasonEccentricityTooLow      SGP4ModelLimitsErrorReason = "eccentricity too low (<= -0.001)"
	ReasonEccentricityTooHigh     SGP4ModelLimitsErrorReason = "eccentricity too high (>= 1.0 - 1e-6)" // Though SGP4 tries to cap this
	ReasonPerturbedEccSqTooHigh   SGP4ModelLimitsErrorReason = "perturbed eccentricity squared (elsq) >= 1.0"
	ReasonBeta2Negative           SGP4ModelLimitsErrorReason = "beta2 (1-e^2) negative"
	ReasonSemiLatusRectumNegative SGP4ModelLimitsErrorReason = "semi-latus rectum (pl) negative"
)

// SGP4ModelLimitsError is returned when SGP4 internal mathematical limits are exceeded,
// often due to extreme orbital parameters (e.g., high drag).
type SGP4ModelLimitsError struct {
	Tsince  float64                    // Time since epoch in minutes when the limit was hit
	Reason  SGP4ModelLimitsErrorReason // The specific limit that was violated
	Value   float64                    // The value that caused the limit violation
	Message string                     // Additional message if any
}

// Error returns the error message for SGP4ModelLimitsError.
func (e *SGP4ModelLimitsError) Error() string {
	return fmt.Sprintf("SGP4 model limits exceeded at tsince %.2f min: %s (value: %.6e). %s", e.Tsince, e.Reason, e.Value, e.Message)
}

// SGPError defines a custom error type for SGP4 related errors.
type SGPError struct {
	msg string
}

func (e *SGPError) Error() string {
	return e.msg
}
