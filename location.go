package sgp4

import (
	"math"
)

// Location represents a ground station or observation point on Earth
type Location struct {
	Latitude  float64 // Latitude in degrees (North positive)
	Longitude float64 // Longitude in degrees (East positive)
	Altitude  float64 // Altitude in meters above sea level
}

// TopocentricCoords represents the position of a satellite relative to
// an observation point in a local horizontal coordinate system:
type TopocentricCoords struct {
	Azimuth   float64 // degrees clockwise from true North (0° to 360°)
	Elevation float64 // degrees above the local horizon (0° to 90°)
	Range     float64 // distance in kilometers from observer to satellite
	RangeRate float64 // rate of change of range in km/s (positive = moving away)
}

// GetLookAngle calculates the topocentric coordinates (azimuth, elevation, range)
// of a satellite relative to the given observation location
func (sv *StateVector) GetLookAngle(loc *Location) *TopocentricCoords {
	if loc == nil {
		return nil
	}

	// Validate latitude range
	if loc.Latitude < -90 || loc.Latitude > 90 {
		return nil
	}

	// Convert observer lat/lon to radians
	latRad := loc.Latitude * deg2rad
	lonRad := loc.Longitude * deg2rad

	// Calculate observer position in Earth-fixed coordinates
	altKm := loc.Altitude / 1000.0

	// WGS-84 Earth radius calculation
	sinLat := math.Sin(latRad)
	cosLat := math.Cos(latRad)
	sinLon := math.Sin(lonRad)
	cosLon := math.Cos(lonRad)

	// Calculate Earth radius at observer's latitude
	C := 1.0 / math.Sqrt(1.0+f*(f-2)*sinLat*sinLat)
	S := (1.0 - f) * (1.0 - f) * C

	// Observer's position vector
	obsX := (re*C + altKm) * cosLat * cosLon
	obsY := (re*C + altKm) * cosLat * sinLon
	obsZ := (re*S + altKm) * sinLat

	// Vector from observer to satellite (topocentric)
	rx := sv.X - obsX
	ry := sv.Y - obsY
	rz := sv.Z - obsZ

	// Convert to topocentric-horizon coordinates (SEZ)
	top_s := sinLat*cosLon*rx + sinLat*sinLon*ry - cosLat*rz
	top_e := -sinLon*rx + cosLon*ry
	top_z := cosLat*cosLon*rx + cosLat*sinLon*ry + sinLat*rz

	// Calculate look angles
	range_ := math.Sqrt(rx*rx + ry*ry + rz*rz)

	// Range rate calculation using velocity components
	rdotx := sv.VX - (-we * sv.Y) // Account for Earth rotation
	rdoty := sv.VY + (we * sv.X)
	rdotz := sv.VZ

	rangeRate := (rx*rdotx + ry*rdoty + rz*rdotz) / range_

	// Calculate azimuth and elevation
	elevation := math.Asin(top_z/range_) * rad2deg

	azimuth := math.Atan2(top_e, top_s) * rad2deg
	if azimuth < 0 {
		azimuth += 360.0
	}

	return &TopocentricCoords{
		Azimuth:   azimuth,
		Elevation: elevation,
		Range:     range_,
		RangeRate: rangeRate,
	}
}
