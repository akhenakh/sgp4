package sgp4

import (
	"strings"
	"testing"
)

func TestSGP4ModelLimitsErrorMessage(t *testing.T) {
	e := &SGP4ModelLimitsError{Reason: ReasonBeta2Negative, Value: -1e-9}
	msg := e.Error()
	if strings.HasSuffix(msg, " ") || strings.HasSuffix(msg, ". ") {
		t.Errorf("message has trailing whitespace: %q", msg)
	}
	if !strings.Contains(msg, string(ReasonBeta2Negative)) {
		t.Errorf("message does not name the reason: %q", msg)
	}

	withDetail := &SGP4ModelLimitsError{Reason: ReasonInvalidMeanMotion, Value: 0, Message: "mean motion must be positive"}
	if got := withDetail.Error(); !strings.HasSuffix(got, "mean motion must be positive") {
		t.Errorf("detail missing from message: %q", got)
	}
}
