package sgp4

import (
	"fmt"
	"math"
	"strings"
)

// SVG plot constants
const (
	svgWidth               = 600
	svgHeight              = 600
	plotMargin             = 50
	plotCenterX            = svgWidth / 2
	plotCenterY            = svgHeight / 2
	plotRadius             = (svgWidth / 2) - plotMargin
	labelFontSize          = 16 // For N,E,S,W
	elevationLabelFontSize = 10 // For 10,30,60 deg labels
	foregroundColor        = "black"
	secondaryColor         = "dimgray"
	gridLineStrokeWidth    = "1"
	pathStrokeWidth        = "3"
	pointRadius            = 5.0
	labelOffsetPoints      = 8.0
)

// polarToCartesian (unchanged)
func polarToCartesian(azimuth, elevation float64, currentPlotRadius float64) (x, y float64) {
	r := currentPlotRadius * (1.0 - elevation/90.0)
	if elevation < 0 {
		r = currentPlotRadius * (1.0 + math.Abs(elevation)/90.0)
	}
	azRad := azimuth * deg2rad
	x = plotCenterX + r*math.Sin(azRad)
	y = plotCenterY - r*math.Cos(azRad)
	return
}

// elevationToColor (unchanged)
func elevationToColor(elevation float64) string {
	if elevation < 0 {
		elevation = 0
	}
	if elevation > 90 {
		elevation = 90
	}
	t := elevation / 90.0
	r := int(255 * (1.0 - t))
	g := int(255 * t)
	b := 0
	return fmt.Sprintf("#%02x%02x%02x", r, g, b)
}

// GeneratePassPolarSVG generates an SVG string representing a satellite pass in a polar view.
func (p *PassDetails) GeneratePassPolarSVG() string {
	if len(p.DataPoints) < 2 {
		return fmt.Sprintf(`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" style="background-color:white;"><rect width="100%%" height="100%%" fill="white"/><text x="50" y="50" fill="black">Not enough data points for pass plot.</text></svg>`, svgWidth, svgHeight)
	}

	var svgBuilder strings.Builder
	svgBuilder.WriteString(fmt.Sprintf(`<svg width="%d" height="%d" xmlns="http://www.w3.org/2000/svg" style="background-color:white;">`, svgWidth, svgHeight))

	// Draw horizon circle (0 degrees elevation)
	svgBuilder.WriteString(fmt.Sprintf(`<circle cx="%f" cy="%f" r="%d" stroke="%s" stroke-width="%s" fill="none"/>`, float64(plotCenterX), float64(plotCenterY), plotRadius, foregroundColor, gridLineStrokeWidth))

	// Draw elevation circles and their labels
	elevationsToMark := []float64{10.0, 30.0, 60.0}
	for _, el := range elevationsToMark {
		radius := plotRadius * (1.0 - el/90.0)
		svgBuilder.WriteString(fmt.Sprintf(`<circle cx="%f" cy="%f" r="%f" stroke="%s" stroke-width="0.5" fill="none" stroke-dasharray="4,4"/>`, float64(plotCenterX), float64(plotCenterY), radius, secondaryColor))

		// Place elevation labels along the North line (0 deg Azimuth), just above their circle
		// X coordinate will be plotCenterX (for North line)
		// Y coordinate is plotCenterY - radius (for the circle edge at North) - small_offset
		labelX := float64(plotCenterX + 5)          // Slightly offset to the right of the North line for readability
		labelY := float64(plotCenterY) - radius - 3 // Place text slightly above its circle line

		svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="%s" font-size="%f" text-anchor="start" dominant-baseline="alphabetic">%d°</text>`, labelX, labelY, secondaryColor, float64(elevationLabelFontSize), int(el)))
	}
	// Zenith label (90 degrees at the center)
	svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="%s" font-size="%f" text-anchor="middle" dominant-baseline="middle">90°</text>`, float64(plotCenterX), float64(plotCenterY), secondaryColor, float64(elevationLabelFontSize)))

	// Draw cardinal direction lines and labels
	cardinalPoints := map[string]float64{
		"N": 0, "E": 90, "S": 180, "W": 270,
	}
	labelRadius := plotRadius + 15.0
	tickMarkLength := 8.0

	for dir, az := range cardinalPoints {
		x_horizon, y_horizon := polarToCartesian(az, 0, plotRadius)
		x_tick_end, y_tick_end := polarToCartesian(az, 0, plotRadius+tickMarkLength)
		svgBuilder.WriteString(fmt.Sprintf(`<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="%s" stroke-width="%s"/>`, x_horizon, y_horizon, x_tick_end, y_tick_end, foregroundColor, gridLineStrokeWidth))

		lx, ly := polarToCartesian(az, 0, labelRadius)

		textAnchor := "middle"
		dominantBaseline := "middle"
		labelAdjustFactor := 0.4 // Factor of font size for minor nudging

		if dir == "N" {
			dominantBaseline = "alphabetic"
			ly -= float64(labelFontSize) * labelAdjustFactor
		}
		if dir == "S" {
			dominantBaseline = "hanging"
			ly += float64(labelFontSize) * labelAdjustFactor
		}
		if dir == "E" {
			textAnchor = "start"
			lx += float64(labelFontSize) * labelAdjustFactor
		}
		if dir == "W" {
			textAnchor = "end"
			lx -= float64(labelFontSize) * labelAdjustFactor
		}

		svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="%s" font-size="%f" text-anchor="%s" dominant-baseline="%s">%s</text>`, lx, ly, foregroundColor, float64(labelFontSize), textAnchor, dominantBaseline, dir))
	}

	// Draw the pass path (unchanged)
	for i := 0; i < len(p.DataPoints)-1; i++ {
		p1 := p.DataPoints[i]
		p2 := p.DataPoints[i+1]
		x1, y1 := polarToCartesian(p1.Azimuth, p1.Elevation, plotRadius)
		x2, y2 := polarToCartesian(p2.Azimuth, p2.Elevation, plotRadius)
		avgElevation := (p1.Elevation + p2.Elevation) / 2.0
		color := elevationToColor(avgElevation)
		svgBuilder.WriteString(fmt.Sprintf(`<line x1="%f" y1="%f" x2="%f" y2="%f" stroke="%s" stroke-width="%s"/>`, x1, y1, x2, y2, color, pathStrokeWidth))
	}

	// Mark AOS, LOS, and Max Elevation points (unchanged)
	if len(p.DataPoints) > 0 {
		aosX, aosY := polarToCartesian(p.AOSObservation.LookAngles.Azimuth, p.AOSObservation.LookAngles.Elevation, plotRadius)
		svgBuilder.WriteString(fmt.Sprintf(`<circle cx="%f" cy="%f" r="%f" fill="darkblue" stroke="black" stroke-width="0.5"/>`, aosX, aosY, pointRadius))
		svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="darkblue" font-size="12" text-anchor="middle" dominant-baseline="text-after-edge">AOS</text>`, aosX, aosY-labelOffsetPoints))

		losX, losY := polarToCartesian(p.LOSObservation.LookAngles.Azimuth, p.LOSObservation.LookAngles.Elevation, plotRadius)
		svgBuilder.WriteString(fmt.Sprintf(`<circle cx="%f" cy="%f" r="%f" fill="darkred" stroke="black" stroke-width="0.5"/>`, losX, losY, pointRadius))
		svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="darkred" font-size="12" text-anchor="middle" dominant-baseline="text-before-edge">LOS</text>`, losX, losY+labelOffsetPoints))

		maxElX, maxElY := polarToCartesian(p.MaxElObservation.LookAngles.Azimuth, p.MaxElObservation.LookAngles.Elevation, plotRadius)
		svgBuilder.WriteString(fmt.Sprintf(`<circle cx="%f" cy="%f" r="%f" fill="lime" stroke="darkgreen" stroke-width="1.5"/>`, maxElX, maxElY, pointRadius+1))
		textAnchorMaxEl := "middle"
		dxMaxEl := 0.0
		dyMaxEl := -(labelOffsetPoints + 2)
		if p.MaxElObservation.LookAngles.Azimuth > 45 && p.MaxElObservation.LookAngles.Azimuth < 135 {
			textAnchorMaxEl = "start"
			dxMaxEl = labelOffsetPoints
			dyMaxEl = 0
		} else if p.MaxElObservation.LookAngles.Azimuth > 225 && p.MaxElObservation.LookAngles.Azimuth < 315 {
			textAnchorMaxEl = "end"
			dxMaxEl = -labelOffsetPoints
			dyMaxEl = 0
		} else if p.MaxElObservation.LookAngles.Azimuth >= 135 && p.MaxElObservation.LookAngles.Azimuth <= 225 {
			dyMaxEl = labelOffsetPoints + 2
		}
		svgBuilder.WriteString(fmt.Sprintf(`<text x="%f" y="%f" fill="darkgreen" font-size="12" text-anchor="%s" dominant-baseline="middle">%.0f°</text>`, maxElX+dxMaxEl, maxElY+dyMaxEl, textAnchorMaxEl, p.MaxElevation))
	}

	svgBuilder.WriteString(`</svg>`)
	return svgBuilder.String()
}
