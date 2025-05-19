package main

import (
	"fmt"
	"log"
	"os"
	"time"

	"github.com/akhenakh/sgp4"
)

func main() {
	tleStr := `ISS (ZARYA)
1 25544U 98067A   25138.37048074  .00007749  00000+0  14567-3 0  9994
2 25544  51.6369  94.7823 0002558 120.7586  15.7840 15.49587957510533`

	tle, err := sgp4.ParseTLE(tleStr)
	if err != nil {
		log.Fatalf("Failed to parse TLE: %v", err)
	}

	observerLat := 46.829853 // Example: Quebec
	observerLng := -71.254028
	observerAltM := 0.0 // meters

	startTime := time.Now().UTC().Add(48 * time.Hour)
	stopTime := startTime.Add(48 * time.Hour) // Predict for the next 24 hours
	stepSeconds := 10                         // Propagation step in seconds

	passes, err := tle.GeneratePasses(observerLat, observerLng, observerAltM, startTime, stopTime, stepSeconds)
	if err != nil {
		log.Fatalf("Error generating passes: %v", err)
	}

	fmt.Printf("Predicted Passes for %s over Lat:%.2f Lon:%.2f:\n", tle.Name, observerLat, observerLng)
	if len(passes) == 0 {
		fmt.Println("No passes found in the given time window.")
	} else {
		for i, pass := range passes {
			fmt.Printf("Pass %d:\n", i+1)
			fmt.Printf("  AOS: %s (Az: %.1f), El: %.1f\n", pass.AOS.Local(), pass.AOSAzimuth, pass.AOSObservation.LookAngles.Elevation)
			fmt.Printf("  Max Elevation: %.1f° (Az: %.1f° at %s)\n", pass.MaxElevation, pass.MaxElevationAz, pass.MaxElevationTime.Local())
			fmt.Printf("  LOS: %s (Az: %.1f), El: %.1f\n", pass.LOS.Local(), pass.LOSAzimuth, pass.LOSObservation.LookAngles.Elevation)
			fmt.Printf("  Duration: %v\n", pass.Duration.Truncate(time.Second))

			// Generate and save SVG for the first pass found
			if i == 0 && len(pass.DataPoints) > 1 {
				svgContent := pass.GeneratePassPolarSVG()
				fileName := fmt.Sprintf("pass_%d_polar_plot.svg", i+1)
				err := os.WriteFile(fileName, []byte(svgContent), 0644)
				if err != nil {
					log.Printf("Error writing SVG to file: %v", err)
				} else {
					fmt.Printf("  Polar plot saved to %s\n", fileName)
				}
			}
		}
	}
}
