package spg4

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
// an observation point in a local horizontal coordinate system
type TopocentricCoords struct {
	Azimuth   float64 // Azimuth angle in degrees (clockwise from North)
	Elevation float64 // Elevation angle in degrees above horizon
	Range     float64 // Slant range in kilometers
	RangeRate float64 // Range rate in kilometers per second
}

// GetLookAngle calculates the topocentric coordinates (azimuth, elevation, range)
// of a satellite relative to the given observation location
func (sv *StateVector) GetLookAngle(loc *Location) *TopocentricCoords {
	// Convert observer lat/lon to radians
	latRad := loc.Latitude * deg2rad
	lonRad := loc.Longitude * deg2rad

	// Calculate observer position in Earth-fixed coordinates
	// Convert altitude from meters to kilometers
	altKm := loc.Altitude / 1000.0
	cosLat := math.Cos(latRad)
	sinLat := math.Sin(latRad)
	cosLon := math.Cos(lonRad)
	sinLon := math.Sin(lonRad)

	// Earth radius at observation point using WGS-84 ellipsoid
	c := 1.0 / math.Sqrt(1.0+f*(f-2.0)*sinLat*sinLat)
	s := (1.0 - f) * (1.0 - f) * c
	obsX := (re*c + altKm) * cosLat * cosLon
	obsY := (re*c + altKm) * cosLat * sinLon
	obsZ := (re*s + altKm) * sinLat

	// Calculate topocentric coordinates (station-centered)
	rx := sv.X - obsX
	ry := sv.Y - obsY
	rz := sv.Z - obsZ
	rvx := sv.VX
	rvy := sv.VY
	rvz := sv.VZ

	// Calculate look angles
	range_ := math.Sqrt(rx*rx + ry*ry + rz*rz)
	rangeRate := (rx*rvx + ry*rvy + rz*rvz) / range_

	// Local horizontal coordinates
	top_s := sinLat*cosLon*rx + sinLat*sinLon*ry - cosLat*rz
	top_e := -sinLon*rx + cosLon*ry
	top_z := cosLat*cosLon*rx + cosLat*sinLon*ry + sinLat*rz

	azimuth := math.Atan2(top_e, top_s) * rad2deg
	if azimuth < 0 {
		azimuth += 360
	}

	elevation := math.Asin(top_z/range_) * rad2deg

	return &TopocentricCoords{
		Azimuth:   azimuth,
		Elevation: elevation,
		Range:     range_,
		RangeRate: rangeRate,
	}
}
