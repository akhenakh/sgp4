package sgp4

import (
	"math"
	"time"
)

// ErrLocationNil is returned when the location is nil.
var ErrLocationNil = &SGPError{"location cannot be nil"}

// ErrInvalidLocationLatitude is returned for invalid latitude.
var ErrInvalidLocationLatitude = &SGPError{"invalid location latitude"}

// Location represents a ground station or observation point on Earth
type Location struct {
	Latitude  float64 // Latitude in degrees (North positive)
	Longitude float64 // Longitude in degrees (East positive)
	Altitude  float64 // Altitude in METERS above sea level (as in original)
}

// TopocentricCoords represents the position of a satellite relative to
// an observation point in a local horizontal coordinate system:
type TopocentricCoords struct {
	Azimuth   float64 // degrees clockwise from true North (0° to 360°)
	Elevation float64 // degrees above the local horizon (-90° to 90°)
	Range     float64 // distance in kilometers from observer to satellite
	RangeRate float64 // rate of change of range in km/s (positive = moving away)
}

// SatellitePosition defines the geodetic position of the satellite.
type SatellitePosition struct {
	Latitude  float64   // Satellite latitude in degrees
	Longitude float64   // Satellite longitude in degrees
	Altitude  float64   // Satellite altitude in km above the ellipsoid
	Timestamp time.Time // Timestamp of this position
}

// Observation combines satellite position and look angles from a ground station.
type Observation struct {
	SatellitePos SatellitePosition // Geodetic position of the satellite
	LookAngles   TopocentricCoords // Look angles from the observer to the satellite
}

// PassDataPoint stores the calculated data for a single point during a pass.
type PassDataPoint struct {
	Timestamp time.Time
	Azimuth   float64 // degrees
	Elevation float64 // degrees
	Range     float64 // km
	RangeRate float64 // km/s
}

// PassDetails stores information about a single satellite pass over a ground station.
type PassDetails struct {
	AOS              time.Time       // Acquisition of Signal time
	LOS              time.Time       // Loss of Signal time
	AOSAzimuth       float64         // Azimuth at AOS (degrees)
	LOSAzimuth       float64         // Azimuth at LOS (degrees)
	MaxElevation     float64         // Maximum elevation during the pass (degrees)
	MaxElevationAz   float64         // Azimuth at maximum elevation (degrees)
	MaxElevationTime time.Time       // Time of maximum elevation
	AOSObservation   Observation     // Observation details at AOS
	LOSObservation   Observation     // Observation details at LOS
	MaxElObservation Observation     // Observation details at Max Elevation
	Duration         time.Duration   // Duration of the pass
	DataPoints       []PassDataPoint // Slice of data points for plotting the pass path
}

// GetLookAngle calculates the topocentric coordinates (azimuth, elevation, range)
// of a satellite relative to the given observation location
func (sv *StateVector) GetLookAngle(loc *Location, currentDateTime time.Time) (*Observation, error) {
	if loc == nil {
		return nil, ErrLocationNil
	}

	// Validate latitude range
	if loc.Latitude < -90 || loc.Latitude > 90 {
		return nil, ErrInvalidLocationLatitude
	}
	// Longitude will be wrapped by math.Sin/Cos, altitude can be anything.

	// Convert observer lat/lon to radians
	obsLatRad := loc.Latitude * deg2rad
	obsLonRad := loc.Longitude * deg2rad

	// Calculate observer position in Earth-fixed coordinates (ECEF)
	// This uses WGS-84 constants from the constants.go file for observer ECEF
	// If you want to use SGP4 specific constants for observer, change reObs and fObs here.
	reObs := reWGS84 // Or reSGP4 if you want SGP4 ellipsoid for observer
	fObs := fWGS84   // Or fSGP4

	altKm := loc.Altitude / 1000.0 // Observer altitude in km

	sinObsLat := math.Sin(obsLatRad)
	cosObsLat := math.Cos(obsLatRad)

	// Radius of curvature in prime vertical (N) for observer
	var N_obs float64
	e2Obs := fObs * (2.0 - fObs)
	if math.Abs(1.0-e2Obs*sinObsLat*sinObsLat) < 1e-14 {
		N_obs = reObs / math.Sqrt(1e-14)
	} else {
		N_obs = reObs / math.Sqrt(1.0-e2Obs*sinObsLat*sinObsLat)
	}

	// Observer's ECEF coordinates
	obsXecef := (N_obs + altKm) * cosObsLat * math.Cos(obsLonRad)
	obsYecef := (N_obs + altKm) * cosObsLat * math.Sin(obsLonRad)
	obsZecef := (N_obs*(1.0-e2Obs) + altKm) * sinObsLat

	// GMST for rotating observer ECEF to ECI
	tempEciForGmst := Eci{DateTime: currentDateTime}
	gmst := tempEciForGmst.GreenwichSiderealTime()

	// Rotate observer ECEF to ECI
	cosGmst := math.Cos(gmst)
	sinGmst := math.Sin(gmst)
	obsXeci := obsXecef*cosGmst - obsYecef*sinGmst
	obsYeci := obsXecef*sinGmst + obsYecef*cosGmst
	obsZeci := obsZecef // Z-axis aligns for ECEF and ECI using GMST rotation

	// Vector from observer (ECI) to satellite (ECI)
	rx := sv.X - obsXeci
	ry := sv.Y - obsYeci
	rz := sv.Z - obsZeci

	range_ := math.Sqrt(rx*rx + ry*ry + rz*rz)
	if range_ == 0 {
		range_ = 1e-9
	} // Avoid division by zero

	// To convert to topocentric-horizon (SEZ or ENU), rotate this ECI range vector.
	// The theta for SEZ transform is Local Sidereal Time (LST = GMST + East Longitude)
	lst := gmst + obsLonRad

	sinLatObs := math.Sin(obsLatRad) // Observer's geodetic latitude
	cosLatObs := math.Cos(obsLatRad)
	sinLstObs := math.Sin(lst)
	cosLstObs := math.Cos(lst)

	topS := sinLatObs*cosLstObs*rx + sinLatObs*sinLstObs*ry - cosLatObs*rz
	topE := -sinLstObs*rx + cosLstObs*ry
	topZ := cosLatObs*cosLstObs*rx + cosLatObs*sinLstObs*ry + sinLatObs*rz

	// Azimuth and Elevation
	elRad := math.Asin(topZ / range_)
	azRad := math.Atan2(topE, topS) // Azimuth from North, clockwise: atan2(E, N) often used.
	// Here, topS is "South component", topE is "East component".
	// Azimuth from North: atan2(E,N).
	// Azimuth from South, clockwise: atan2(E,S).

	azimuth := azRad * rad2deg
	if azimuth < 0.0 {
		azimuth += 360.0
	}
	elevation := elRad * rad2deg

	// Range Rate calculation
	// Observer velocity in ECI (due to Earth rotation)
	// omega_earth_rad_sec = 7.2921150e-5 (from constants.go we)
	obsVXeci := -we * obsYeci // vx = -omega_earth * y_eci
	obsVYeci := we * obsXeci  // vy =  omega_earth * x_eci
	obsVZeci := 0.0           // vz = 0

	deltaVX := sv.VX - obsVXeci
	deltaVY := sv.VY - obsVYeci
	deltaVZ := sv.VZ - obsVZeci

	rangeRate := (rx*deltaVX + ry*deltaVY + rz*deltaVZ) / range_

	// Geodetic position of the satellite (already ECI, convert to Geodetic)
	satEci := Eci{DateTime: currentDateTime, Position: Vector{X: sv.X, Y: sv.Y, Z: sv.Z}}
	satLat, satLon, satAltKm := satEci.ToGeodetic()

	return &Observation{
		SatellitePos: SatellitePosition{
			Latitude:  satLat,
			Longitude: satLon,
			Altitude:  satAltKm,
			Timestamp: currentDateTime,
		},
		LookAngles: TopocentricCoords{
			Azimuth:   azimuth,
			Elevation: elevation,
			Range:     range_,
			RangeRate: rangeRate,
		},
	}, nil
}
