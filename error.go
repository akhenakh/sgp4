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
	ReasonDeepSpaceUnsupported    SGP4ModelLimitsErrorReason = "deep-space orbit (period >= 225 min) not supported: SDP4 is not implemented"
	ReasonInvalidMeanMotion       SGP4ModelLimitsErrorReason = "mean motion is not a finite positive number"
	ReasonInvalidEccentricity     SGP4ModelLimitsErrorReason = "eccentricity is outside [0, 1)"
	ReasonKeplerNotConverged      SGP4ModelLimitsErrorReason = "Kepler's equation did not converge"
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
	msg := fmt.Sprintf("SGP4 model limits exceeded at tsince %.2f min: %s (value: %.6e).", e.Tsince, e.Reason, e.Value)
	if e.Message != "" {
		msg += " " + e.Message
	}
	return msg
}

// SGPError defines a custom error type for SGP4 related errors.
type SGPError struct {
	msg string
}

func (e *SGPError) Error() string {
	return e.msg
}

// UnsupportedEphemerisTypeError is returned when an element set declares an
// ephemeris type this package cannot propagate. Only SGP4 element sets (types 0
// and 2) are supported: type 3 is PPT3, type 4 is SGP4-XP and type 6 is Special
// Perturbations, each a different theory, so propagating them with SGP4 would
// return a silently wrong position.
type UnsupportedEphemerisTypeError struct {
	Type int    // EPHEMERIS_TYPE from the element set
	Name string // Human-readable theory name, empty when unknown
}

// Error returns the error message for UnsupportedEphemerisTypeError.
func (e *UnsupportedEphemerisTypeError) Error() string {
	if e.Name != "" {
		return fmt.Sprintf("sgp4: ephemeris type %d (%s) is not supported: only SGP4 element sets (types 0 and 2) can be propagated", e.Type, e.Name)
	}
	return fmt.Sprintf("sgp4: ephemeris type %d is not supported: only SGP4 element sets (types 0 and 2) can be propagated", e.Type)
}

// ephemerisTypeName returns the CCSDS name of an EPHEMERIS_TYPE value.
func ephemerisTypeName(t int) string {
	switch t {
	case 0:
		return "SGP"
	case 2:
		return "SGP4"
	case 3:
		return "PPT3"
	case 4:
		return "SGP4-XP"
	case 6:
		return "Special Perturbations"
	}
	return ""
}
