# Go SGP4 Satellite Propagation Library


A Go (Golang) library for satellite propagation using the SGP4 model.

This library allows you to predict the position and velocity of Earth-orbiting satellites given their orbital elements, typically from Two-Line Element (TLE) sets or Orbit Mean-elements Message (OMM) JSON data.


This an AI assisted port of the spg4 library from **Daniel Warner** ([github.com/dnwrnr/sgp4/](https://github.com/dnwrnr/sgp4/)).

## Features

*   **TLE Parsing:** Parses standard TLE format (2-line or 3-line with satellite name) into structured data. Catalog numbers use the Alpha-5 encoding above 99999 (`A0000` is 100000, `Z9999` is 339999). Includes checksum validation.
*   **TLE Writing:** `TLELines` renders lines 1 and 2 with computed checksums, and `TLEString` the three-line form with the name. Values the fixed-column layout cannot hold (a catalog number above 339999, an out-of-range epoch, angle or eccentricity) return an error instead of being truncated.
*   **OMM / GP data:** Parses OMM JSON (`ParseOMMs`), OMM/General-Perturbations CSV (`ParseOMMsCSV`), or either (`ParseOMMsReader`), and converts to `TLE` for the propagator (`OMM.ToTLE`). `NORAD_CAT_ID` is a plain integer and holds up to nine digits, beyond the TLE ceiling of 339999. `TLE.ToOMM` converts back, and `OMM.ToTLE` tolerates a missing/invalid international designator.
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

**Note:** This library currently focuses on near-Earth propagations (SGP4) and does not yet implement the deep-space corrections (SDP4) for objects with periods greater than 225 minutes. Element sets that declare another theory -- `EPHEMERIS_TYPE` 3 (PPT3), 4 (SGP4-XP) or 6 (Special Perturbations) -- are refused with an `UnsupportedEphemerisTypeError` instead of being propagated with SGP4.

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

### 2. Parsing OMM / General-Perturbations data

The OMM JSON and CelesTrak GP CSV formats carry the catalog number as an integer, so they are not limited to the five-character TLE field:

```go
omms, err := sgp4.ParseOMMsReader(resp.Body) // JSON or CSV, auto-detected
if err != nil {
	log.Fatal(err)
}
for i := range omms {
	tle, err := omms[i].ToTLE() // ready for FindPosition / GeneratePasses
	if err != nil {
		continue
	}
	fmt.Println(tle.Name, tle.SatelliteNumber)
}
```

`ParseOMMsCSV` reads a header row and maps columns by name, so column order does not matter. Convert the other way with `tle.ToOMM()`.

### 3. Writing TLE text

```go
lines, err := tle.TLELines() // {line 1, line 2}, checksums computed
if err != nil {
	log.Fatal(err) // a value the layout cannot hold, e.g. catalog number > 339999
}
fmt.Println(lines[0])
fmt.Println(lines[1])

text, err := tle.TLEString() // three-line form: name padded to 24 columns, then the two lines
fmt.Println(text)
```

### 4. Propagating Satellite Position

```go
// Propagate to a specific time
tsince := 60.0 // minutes since the TLE epoch
eciState, err := tle.FindPosition(tsince)
// Or propagate to an absolute time:
// eciState, err := tle.FindPositionAtTime(tle.EpochTime().Add(time.Hour))
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

### 5. Calculating Look Angles

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

### 6. Generating Pass Predictions

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

## Catalog numbers

A catalog number appears in three forms:

* in TLE text, five characters: digits below 100000, the Alpha-5 letter form from 100000 to 339999 (`A0000`-`Z9999`, with `I` and `O` unused);
* in OMM/GP JSON and CSV, a plain integer, up to nine digits;
* in `TLE.SatelliteNumber`, the decoded integer.

`ParseTLE` decodes Alpha-5 fields, and `ParseAlpha5` / `FormatAlpha5` convert between a five-character field and the number directly. Objects above 339999 have no TLE representation, so they have to be carried as OMM/GP data; writing such a number with `TLELines` returns an error.

## Use This Package

- This is the reason for this package, I was maintaining an iOS (ObjC) app called [SatSaT](https://satsat.inair.space/) but cost was too high, so I created a [web app](https://websat.inair.space/) with the same ideas.

## License

Apache 2.0


## TODO
- reuse the same mechanism to find the exact AOC a sampling more than in the lookup
