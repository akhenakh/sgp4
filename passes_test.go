package spg4

import (
	"testing"
	"time"
)

func TestFindPasses(t *testing.T) {
	// ISS TLE from a known epoch
	issTLE := `1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	// Test location: McDonald Observatory
	loc := &Location{
		Latitude:  30.6715,
		Longitude: -104.0227,
		Altitude:  2070, // meters
	}

	epoch := tle.EpochTime()
	start := epoch
	end := epoch.Add(24 * time.Hour)

	passes, err := tle.FindPasses(loc, start, end, 10.0) // 10° minimum elevation
	if err != nil {
		t.Fatalf("FindPasses failed: %v", err)
	}

	// Basic sanity checks for pass prediction
	if len(passes) == 0 {
		t.Error("No passes found in 24-hour period")
	}

	for i, pass := range passes {
		t.Run("Pass validity checks", func(t *testing.T) {
			checkPass(t, pass, i)
		})
	}
}

func checkPass(t *testing.T, pass *Pass, index int) {
	// Time sequence check
	if !pass.AOS.Before(pass.TCA) {
		t.Errorf("Pass %d: AOS (%v) not before TCA (%v)", index, pass.AOS, pass.TCA)
	}
	if !pass.TCA.Before(pass.LOS) {
		t.Errorf("Pass %d: TCA (%v) not before LOS (%v)", index, pass.TCA, pass.LOS)
	}

	// Elevation checks
	if pass.MaxEl <= 0 || pass.MaxEl > 90 {
		t.Errorf("Pass %d: Invalid maximum elevation: %.2f", index, pass.MaxEl)
	}
	if pass.MaxAzEl.Elevation != pass.MaxEl {
		t.Errorf("Pass %d: MaxEl (%.2f) doesn't match MaxAzEl elevation (%.2f)",
			index, pass.MaxEl, pass.MaxAzEl.Elevation)
	}

	// Azimuth checks
	checkAzimuth := func(az float64, point string) {
		if az < 0 || az >= 360 {
			t.Errorf("Pass %d: Invalid %s azimuth: %.2f", index, point, az)
		}
	}
	checkAzimuth(pass.StartAzEl.Azimuth, "start")
	checkAzimuth(pass.MaxAzEl.Azimuth, "max")
	checkAzimuth(pass.EndAzEl.Azimuth, "end")

	// Range checks
	checkRange := func(rng float64, point string) {
		if rng <= 0 {
			t.Errorf("Pass %d: Invalid %s range: %.2f", index, point, rng)
		}
	}
	checkRange(pass.StartAzEl.Range, "start")
	checkRange(pass.MaxAzEl.Range, "max")
	checkRange(pass.EndAzEl.Range, "end")
}

func TestFindPassesEdgeCases(t *testing.T) {
	// ISS TLE
	issTLE := `1 25544U 98067A   25025.00048859  .00033214  00000+0  57704-3 0  9996
2 25544  51.6377 296.2827 0003104 141.8447 313.9175 15.50506992492954`

	tle, err := ParseTLE(issTLE)
	if err != nil {
		t.Fatalf("Failed to parse TLE: %v", err)
	}

	epoch := tle.EpochTime()

	tests := []struct {
		name     string
		loc      *Location
		start    time.Time
		end      time.Time
		minEl    float64
		wantPass bool
	}{
		{
			name: "Very short time window",
			loc: &Location{
				Latitude:  30.6715,
				Longitude: -104.0227,
				Altitude:  2070,
			},
			start:    epoch,
			end:      epoch.Add(10 * time.Minute),
			minEl:    10.0,
			wantPass: false, // May be true depending on exact timing
		},
		{
			name: "High minimum elevation",
			loc: &Location{
				Latitude:  30.6715,
				Longitude: -104.0227,
				Altitude:  2070,
			},
			start:    epoch,
			end:      epoch.Add(24 * time.Hour),
			minEl:    85.0,
			wantPass: false,
		},
		{
			name: "Location near pole",
			loc: &Location{
				Latitude:  89.9,
				Longitude: 0.0,
				Altitude:  0,
			},
			start:    epoch,
			end:      epoch.Add(24 * time.Hour),
			minEl:    10.0,
			wantPass: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			passes, err := tle.FindPasses(tt.loc, tt.start, tt.end, tt.minEl)
			if err != nil {
				t.Fatalf("FindPasses failed: %v", err)
			}

			if tt.wantPass && len(passes) == 0 {
				t.Error("Expected at least one pass, got none")
			}

			for _, pass := range passes {
				if pass.MaxEl < tt.minEl {
					t.Errorf("Found pass with max elevation %.2f below minimum %.2f",
						pass.MaxEl, tt.minEl)
				}
			}
		})
	}
}
