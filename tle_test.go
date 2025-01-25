package spg4

import (
	"math"
	"testing"
)

func TestParseTLE(t *testing.T) {
	issTLE := `1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse ISS TLE: %v", err)
	}

	// Test Line 1 fields
	tests := []struct {
		name    string
		got     interface{}
		want    interface{}
		epsilon float64
		compare func(got, want interface{}, epsilon float64) bool
	}{
		{"Satellite Number", tle.SatelliteNumber, 25544, 0, compareExact},
		{"Classification", string(tle.Classification), "U", 0, compareExact},
		{"International ID", tle.International, "98067A", 0, compareExact},
		{"Epoch Year", tle.EpochYear, 2025, 0, compareExact},
		{"Epoch Day", tle.EpochDay, 25.00048859, 1e-8, compareFloat},
		{"Mean Motion Dot", tle.MeanMotionDot, 0.00033214, 1e-8, compareFloat},
		{"B* Drag Term", tle.Bstar, 0.00057704, 1e-8, compareFloat},

		// Line 2 fields
		{"Inclination", tle.Inclination, 51.6377, 1e-4, compareFloat},
		{"Right Ascension", tle.RightAscension, 296.2827, 1e-4, compareFloat},
		{"Eccentricity", tle.Eccentricity, 0.0003104, 1e-7, compareFloat},
		{"Argument of Perigee", tle.ArgOfPerigee, 141.8447, 1e-4, compareFloat},
		{"Mean Anomaly", tle.MeanAnomaly, 313.9175, 1e-4, compareFloat},
		{"Mean Motion", tle.MeanMotion, 15.50506992, 1e-8, compareFloat},
		{"Revolution Number", tle.RevolutionNumber, 49295, 0, compareExact},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if !tt.compare(tt.got, tt.want, tt.epsilon) {
				t.Errorf("%s = %v, want %v", tt.name, tt.got, tt.want)
			}
		})
	}
}

func compareExact(got, want interface{}, _ float64) bool {
	return got == want
}

func compareFloat(got, want interface{}, epsilon float64) bool {
	g, ok1 := got.(float64)
	w, ok2 := want.(float64)
	if !ok1 || !ok2 {
		return false
	}
	return math.Abs(g-w) < epsilon
}

func TestInvalidTLE(t *testing.T) {
	tests := []struct {
		name    string
		tle     string
		wantErr bool
	}{
		{
			name:    "Empty input",
			tle:     "",
			wantErr: true,
		},
		{
			name:    "Single line",
			tle:     "1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996",
			wantErr: true,
		},
		{
			name: "Invalid line length",
			tle: `1 25544U 98067A   25025.00048859
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`,
			wantErr: true,
		},
		{
			name: "Invalid line numbers",
			tle: `3 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`,
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			_, err := ParseTLE(tt.tle)
			if (err != nil) != tt.wantErr {
				t.Errorf("ParseTLE() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}
