package sgp4

import (
	"math"
	"testing"
)

// TestValladoSGP4Verification checks FindPosition against Vallado's SGP4-VER suite
// (AIAA 2006-6753, "Revisiting Spacetrack Report #3") — the canonical SGP4 verification
// vectors. Reference TEME position (km) / velocity (km/s) are the published tcppver.out
// values. These are near-earth cases (period < 225 min), which this library supports;
// the deep-space cases in SGP4-VER need the SDP4 model.
//
// Observed agreement is ~1e-8 km / ~1e-9 km/s; the tolerances below are deliberately
// loose to stay robust across platforms while still catching any real (km-scale)
// propagation regression.
func TestValladoSGP4Verification(t *testing.T) {
	type point struct{ tsince, x, y, z, vx, vy, vz float64 }
	cases := []struct {
		name   string
		l1, l2 string
		pts    []point
	}{
		{
			name: "NORAD 5 (TEME example)",
			l1:   "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
			l2:   "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
			pts: []point{
				{0.0, 7022.46529266, -1400.08296755, 0.03995155, 1.893841015, 6.405893759, 4.534807250},
				{360.0, -7154.03120202, -3783.17682504, -3536.19412294, 4.741887409, -4.151817765, -2.093935425},
				{720.0, -7134.59340119, 6531.68641334, 3260.27186483, -4.113793027, -2.911922039, -2.557327851},
				{1440.0, -938.55923943, -6268.18748831, -4294.02924751, 7.536105209, -0.427127707, 0.989878080},
			},
		},
		{
			name: "NORAD 28057 (sun-sync LEO)",
			l1:   "1 28057U 03049A   06177.78615833  .00000060  00000-0  35940-4 0  1836",
			l2:   "2 28057  98.4283 247.6961 0000884  88.1964 271.9322 14.35478080140550",
			pts: []point{
				{0.0, -2715.28237486, -6619.26436889, -0.01341443, -1.008587273, 0.422782003, 7.385272942},
				{120.0, -1816.87920942, -1835.78762132, 6661.07926465, 2.325140071, 6.655669329, 2.463394512},
				{240.0, 1483.17364291, 5395.21248786, 4448.65907172, 2.560540387, 4.039025766, -5.736648561},
				{360.0, 2801.25607157, 5455.03931333, -3692.12865695, -0.595095864, -3.951923117, -6.298799125},
			},
		},
		{
			name: "NORAD 88888 (Spacetrack Report #3)",
			l1:   "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87",
			l2:   "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058",
			pts: []point{
				{0.0, 2328.96975262, -5995.22051338, 1719.97297192, 2.912073281, -0.983417956, -7.090816210},
				{120.0, 1020.69234558, 2286.56260634, -6191.55565927, -3.746543902, 6.467532721, 1.827985678},
				{240.0, -3226.54349155, 3503.70977525, 4532.80979343, 1.000992116, -5.788042888, 5.162585826},
				{360.0, 2456.10706533, -6071.93855503, 1222.89768554, 2.679390040, -0.448290811, -7.228792155},
			},
		},
	}
	const posTol, velTol = 1e-4, 1e-6 // km, km/s

	mag := func(a, b, c float64) float64 { return math.Sqrt(a*a + b*b + c*c) }
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			tle, err := ParseTLELines([]string{c.l1, c.l2})
			if err != nil {
				t.Fatalf("ParseTLELines: %v", err)
			}
			for _, p := range c.pts {
				eci, err := tle.FindPosition(p.tsince)
				if err != nil {
					t.Fatalf("FindPosition(%.0f min): %v", p.tsince, err)
				}
				if dp := mag(eci.Position.X-p.x, eci.Position.Y-p.y, eci.Position.Z-p.z); dp > posTol {
					t.Errorf("t=%.0f min: TEME position off by %.3e km (tol %.0e)", p.tsince, dp, posTol)
				}
				if dv := mag(eci.Velocity.X-p.vx, eci.Velocity.Y-p.vy, eci.Velocity.Z-p.vz); dv > velTol {
					t.Errorf("t=%.0f min: TEME velocity off by %.3e km/s (tol %.0e)", p.tsince, dv, velTol)
				}
			}
		})
	}
}
