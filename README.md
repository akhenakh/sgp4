# Go SGP4 Satellite Propagation Library

This Go package provides functionality to parse Two-Line Element (TLE) sets and propagate satellite orbits using an SGP4-compatible model. The goal is to closely match the reference SGP4 propagator used by NORAD and Space-Track for near-Earth objects.

This an AI assisted port of the spg4 library from **Daniel Warner** ([github.com/dnwrnr/sgp4/](https://github.com/dnwrnr/sgp4/)).

## Features

*   **TLE Parsing:** Parses standard TLE format (2-line or 3-line with satellite name) into structured data. Includes checksum validation.
*   **SGP4 Propagation:**
    *   Initializes orbital elements and internal constants according to SGP4 methodology.
    *   Propagates satellite ECI (Earth-Centered Inertial) position and velocity over time:
        *   `FindPosition(tsinceMinutes float64)`: Propagates relative to TLE epoch.
        *   `FindPositionAtTime(t time.Time)`: Propagates to an absolute UTC time.
    *   Includes secular effects of J2, J3, J4 gravitational harmonics.
    *   Implements atmospheric drag model based on the B\* term.
    *   Applies short-period periodic perturbations due to Earth's oblateness (J2).
*   **Orbital Characteristics:**
    *   `IsGeostationary() bool`: Heuristically determines if a satellite is likely geostationary based on its TLE elements.
*   **Coordinate Transformations:**
    *   Converts ECI coordinates to Geodetic (latitude, longitude, altitude) using an Earth model aligned with common SGP4 implementations.
    *   Calculates topocentric look angles (azimuth, elevation, range, range rate) from an observer to the satellite.
*   **Pass Prediction:**
    *   Generates satellite pass predictions (AOS, Max Elevation, LOS times and corresponding look angles) over a specified ground station and time window.

## Accuracy and Reference

Significant portions of the SGP4 mathematical model and constant initialization logic are based on well-established C++ implementations, such as the one by **Daniel Warner** ([github.com/dnwrnr/sgp4/](https://github.com/dnwrnr/sgp4/)). This Go implementation strives for compatibility and comparable results with such reference libraries.

**Note:** This library currently focuses on near-Earth propagations (SGP4) and does not yet implement the deep-space corrections (SDP4) for objects with periods greater than 225 minutes.

## Basic Usage

### 1. Parsing a TLE

```go
package main

import (
	"fmt"
	"log"
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
	fmt.Printf("Successfully parsed TLE for: %s\n", tle.Name)
	fmt.Printf("Epoch Time: %v\n", tle.EpochTime())
}
```

### 2. Propagating Satellite Position

```go
// Propagate to a specific time
targetTime := tle.EpochTime().Add(60 * time.Minute)
eciState, err := tle.FindPositionAtTime(targetTime)
// Or propagate by minutes from epoch:
// eciState, err := tle.FindPosition(60.0)
if err != nil {
	log.Fatalf("Failed to propagate position: %v", err)
}

fmt.Printf("ECI Position at T+%.0f min (X,Y,Z km): %.3f, %.3f, %.3f\n",
	tsince, eciState.Position.X, eciState.Position.Y, eciState.Position.Z)
fmt.Printf("ECI Velocity at T+%.0f min (VX,VY,VZ km/s): %.3f, %.3f, %.3f\n",
	tsince, eciState.Velocity.X, eciState.Velocity.Y, eciState.Velocity.Z)

// Convert to Geodetic
lat, lon, alt := eciState.ToGeodetic()
fmt.Printf("Geodetic at T+%.0f min (Lat,Lon,Alt km): %.3f deg, %.3f deg, %.3f km\n",
	tsince, lat, lon, alt)

// validate is geostationary
if tle.IsGeostationary() {
    fmt.Printf("%s is likely a geostationary satellite.\n", tle.Name)
} else {
    fmt.Printf("%s is not classified as geostationary by this check.\n", tle.Name)
}
```

### 3. Calculating Look Angles

```go
// ... (inside main or another function, after getting eciState) ...
observer := &sgp4.Location{
	Latitude:  40.0,  // degrees North
	Longitude: -75.0, // degrees West
	Altitude:  100.0, // meters above sea level
}

// Construct StateVector for GetLookAngle (or modify GetLookAngle to take Eci.Position/Velocity)
sv := &sgp4.StateVector{
    X: eciState.Position.X, Y: eciState.Position.Y, Z: eciState.Position.Z,
    VX: eciState.Velocity.X, VY: eciState.Velocity.Y, VZ: eciState.Velocity.Z,
}

observation, err := sv.GetLookAngle(observer, eciState.DateTime)
if err != nil {
	log.Fatalf("Failed to get look angles: %v", err)
}

fmt.Printf("Look Angles (Az,El,Range km,RangeRate km/s): %.1f deg, %.1f deg, %.1f km, %.2f km/s\n",
	observation.LookAngles.Azimuth,
	observation.LookAngles.Elevation,
	observation.LookAngles.Range,
	observation.LookAngles.RangeRate)
```

### 4. Generating Pass Predictions

```go
// ... (inside main or another function) ...
startTime := time.Now().UTC()
stopTime := startTime.Add(24 * time.Hour) // Predict for the next 24 hours
stepSeconds := 30                         // Propagation step in seconds

passes, err := tle.GeneratePasses(observer.Latitude, observer.Longitude, observer.Altitude, startTime, stopTime, stepSeconds)
if err != nil {
	log.Fatalf("Error generating passes: %v", err)
}

fmt.Printf("\nPredicted Passes for %s over Lat:%.2f Lon:%.2f:\n", tle.Name, observer.Latitude, observer.Longitude)
if len(passes) == 0 {
	fmt.Println("No passes found in the given time window.")
}
for i, pass := range passes {
	fmt.Printf("Pass %d:\n", i+1)
	fmt.Printf("  AOS: %s (Az: %.1f°)\n", pass.AOS.Local(), pass.AOSAzimuth)
	fmt.Printf("  Max Elevation: %.1f° (Az: %.1f° at %s)\n", pass.MaxElevation, pass.MaxElevationAz, pass.MaxElevationTime.Local())
	fmt.Printf("  LOS: %s (Az: %.1f°)\n", pass.LOS.Local(), pass.LOSAzimuth)
	fmt.Printf("  Duration: %v\n", pass.Duration.Truncate(time.Second))
}
```

## License

Apache 2.0


## TODO
- reuse the same mechanism to find the exact AOC a sampling more than in the lookup
- make the min elevation a parameter
