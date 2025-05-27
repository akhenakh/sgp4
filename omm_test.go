package sgp4

import (
	"math"
	"testing"
	"time"
)

// ommJsonExample constant remains as previously defined (containing ISS ZARYA as the first entry)
const ommJsonExample = `
[{"OBJECT_NAME":"ISS (ZARYA)","OBJECT_ID":"1998-067A","EPOCH":"2025-05-26T13:06:57.824640","MEAN_MOTION":15.4975272,"ECCENTRICITY":0.0002241,"INCLINATION":51.6382,"RA_OF_ASC_NODE":54.2937,"ARG_OF_PERICENTER":147.4648,"MEAN_ANOMALY":271.6158,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":25544,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":51180,"BSTAR":0.00019155,"MEAN_MOTION_DOT":0.00010397,"MEAN_MOTION_DDOT":0},{"OBJECT_NAME":"CSS (TIANHE)","OBJECT_ID":"2021-035A","EPOCH":"2025-05-25T23:00:12.248640","MEAN_MOTION":15.62412324,"ECCENTRICITY":0.0005017,"INCLINATION":41.463,"RA_OF_ASC_NODE":155.4996,"ARG_OF_PERICENTER":337.345,"MEAN_ANOMALY":22.7167,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":48274,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":23268,"BSTAR":0.00015624,"MEAN_MOTION_DOT":0.00013949,"MEAN_MOTION_DDOT":0},{"OBJECT_NAME":"FREGAT DEB","OBJECT_ID":"2011-037PF","EPOCH":"2025-05-19T00:59:35.639808","MEAN_MOTION":12.28834273,"ECCENTRICITY":0.0869949,"INCLINATION":51.6315,"RA_OF_ASC_NODE":92.6347,"ARG_OF_PERICENTER":128.5677,"MEAN_ANOMALY":239.6424,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":49271,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":17632,"BSTAR":0.03654,"MEAN_MOTION_DOT":0.00014961,"MEAN_MOTION_DDOT":0}]
`

func TestParseOMMs(t *testing.T) {
	omms, err := ParseOMMs([]byte(ommJsonExample))
	if err != nil {
		t.Fatalf("ParseOMMs failed: %v", err)
	}

	if len(omms) != 3 { // Example has more, but test string is truncated to 3 for brevity.
		t.Errorf("Expected 3 OMM objects, got %d", len(omms))
	}

	// Test first OMM (ISS ZARYA)
	issOMM := omms[0]
	if issOMM.ObjectName != "ISS (ZARYA)" {
		t.Errorf("Expected ObjectName 'ISS (ZARYA)', got '%s'", issOMM.ObjectName)
	}
	if issOMM.NoradCatID != 25544 {
		t.Errorf("Expected NoradCatID 25544, got %d", issOMM.NoradCatID)
	}
	if issOMM.EpochStr != "2025-05-26T13:06:57.824640" {
		t.Errorf("Expected EpochStr '2025-05-26T13:06:57.824640', got '%s'", issOMM.EpochStr)
	}
	expectedMeanMotion := 15.4975272
	if math.Abs(issOMM.MeanMotion-expectedMeanMotion) > 1e-9 {
		t.Errorf("Expected MeanMotion %f, got %f", expectedMeanMotion, issOMM.MeanMotion)
	}
}

func TestOMMToTLE(t *testing.T) {
	omms, err := ParseOMMs([]byte(ommJsonExample))
	if err != nil {
		t.Fatalf("ParseOMMs failed for ToTLE test: %v", err)
	}
	if len(omms) == 0 {
		t.Fatal("No OMMs parsed for ToTLE test.")
	}

	issOMM := omms[0] // ISS (ZARYA)
	tle, err := issOMM.ToTLE()
	if err != nil {
		t.Fatalf("issOMM.ToTLE() failed: %v", err)
	}

	// Verify TLE fields
	if tle.Name != "ISS (ZARYA)" {
		t.Errorf("TLE Name: got '%s', want 'ISS (ZARYA)'", tle.Name)
	}
	if tle.SatelliteNumber != 25544 {
		t.Errorf("TLE SatelliteNumber: got %d, want 25544", tle.SatelliteNumber)
	}
	if tle.Classification != 'U' {
		t.Errorf("TLE Classification: got '%c', want 'U'", tle.Classification)
	}
	if tle.International != "98067A" { // From "1998-067A"
		t.Errorf("TLE International: got '%s', want '98067A'", tle.International)
	}

	// Epoch verification for ISS
	// OMM EPOCH: "2025-05-26T13:06:57.824640"
	// Expected TLE EpochYear: 2025
	// Expected TLE EpochDay: 146.546502599999987 (from float64(47217824640000) / float64(86400000000000))
	if tle.EpochYear != 2025 {
		t.Errorf("TLE EpochYear: got %d, want 2025", tle.EpochYear)
	}
	expectedEpochDayISS := 146.546502599999987              // Value from float64 conversion of OMM epoch string for ISS
	if math.Abs(tle.EpochDay-expectedEpochDayISS) > 1e-15 { // Tighter tolerance for direct check
		t.Errorf("ISS TLE EpochDay: got %.17f, want %.17f, diff %.17e", tle.EpochDay, expectedEpochDayISS, math.Abs(tle.EpochDay-expectedEpochDayISS))
	}

	// Verify EpochTime() reconstruction for ISS
	expectedEpochTimeISS := time.Date(2025, 5, 26, 13, 6, 57, 824640000, time.UTC)
	if !tle.EpochTime().Equal(expectedEpochTimeISS) {
		t.Errorf("ISS TLE EpochTime(): got %v, want %v. Diff (ns): %d",
			tle.EpochTime(), expectedEpochTimeISS, tle.EpochTime().Sub(expectedEpochTimeISS).Nanoseconds())
	}

	if math.Abs(tle.MeanMotionDot-0.00010397) > 1e-9 {
		t.Errorf("TLE MeanMotionDot: got %f, want 0.00010397", tle.MeanMotionDot)
	}
	if math.Abs(tle.MeanMotionDot2-0.0) > 1e-9 {
		t.Errorf("TLE MeanMotionDot2: got %f, want 0.0", tle.MeanMotionDot2)
	}
	if math.Abs(tle.Bstar-0.00019155) > 1e-9 {
		t.Errorf("TLE Bstar: got %f, want 0.00019155", tle.Bstar)
	}
	if tle.ElementNumber != 999 {
		t.Errorf("TLE ElementNumber: got %d, want 999", tle.ElementNumber)
	}

	// Line 2 fields
	if math.Abs(tle.Inclination-51.6382) > 1e-9 {
		t.Errorf("TLE Inclination: got %f, want 51.6382", tle.Inclination)
	}
	if math.Abs(tle.RightAscension-54.2937) > 1e-9 {
		t.Errorf("TLE RightAscension: got %f, want 54.2937", tle.RightAscension)
	}
	if math.Abs(tle.Eccentricity-0.0002241) > 1e-9 {
		t.Errorf("TLE Eccentricity: got %f, want 0.0002241", tle.Eccentricity)
	}
	if math.Abs(tle.ArgOfPerigee-147.4648) > 1e-9 {
		t.Errorf("TLE ArgOfPerigee: got %f, want 147.4648", tle.ArgOfPerigee)
	}
	if math.Abs(tle.MeanAnomaly-271.6158) > 1e-9 {
		t.Errorf("TLE MeanAnomaly: got %f, want 271.6158", tle.MeanAnomaly)
	}
	if math.Abs(tle.MeanMotion-15.4975272) > 1e-9 {
		t.Errorf("TLE MeanMotion: got %f, want 15.4975272", tle.MeanMotion)
	}
	if tle.RevolutionNumber != 51180 {
		t.Errorf("TLE RevolutionNumber: got %d, want 51180", tle.RevolutionNumber)
	}

	// FREGAT DEB - Test with different OBJECT_ID format and non-zero BSTAR exponent in TLE terms
	fregatOMM := omms[2]
	fregatTLE, err := fregatOMM.ToTLE()
	if err != nil {
		t.Fatalf("fregatOMM.ToTLE() failed: %v", err)
	}
	if fregatTLE.International != "11037PF" { // From "2011-037PF"
		t.Errorf("FREGAT TLE International: got '%s', want '11037PF'", fregatTLE.International)
	}
	if math.Abs(fregatTLE.Bstar-0.03654) > 1e-9 {
		t.Errorf("FREGAT TLE Bstar: got %f, want 0.03654", fregatTLE.Bstar)
	}

	// Test epoch parsing for CSS (TIANHE)
	// OMM EPOCH:"2025-05-25T23:00:12.248640"
	cssOMM := omms[1]
	cssTLE, err := cssOMM.ToTLE()
	if err != nil {
		t.Fatalf("cssOMM.ToTLE() failed: %v", err)
	}
	if cssTLE.EpochYear != 2025 {
		t.Errorf("CSS TLE EpochYear: got %d, want 2025", cssTLE.EpochYear)
	}
	// Expected EpochDay for CSS, as calculated by ommEpochToTleEpoch:
	// 145.0 + (float64(82812248640000) / float64(86400000000000))
	// = 145.0 + 0.95847509999999998692... (actual float64 result is closer to ...98692 not ...00002)
	// The value 145.95847509999998692 is what `ommEpochToTleEpoch` will produce for this input.
	expectedCSSEpochDay := 145.95847509999998692
	if math.Abs(cssTLE.EpochDay-expectedCSSEpochDay) > 1e-15 { // Use a very tight tolerance here
		t.Errorf("CSS TLE EpochDay: got %.20f, want %.20f, diff %.20e", cssTLE.EpochDay, expectedCSSEpochDay, math.Abs(cssTLE.EpochDay-expectedCSSEpochDay))
	}
}

func TestOmmEpochToTleEpochVariations(t *testing.T) {
	tests := []struct {
		name           string
		epochStr       string
		wantYear       int
		wantDay        float64 // This should be the value ommEpochToTleEpoch is expected to produce
		wantErr        bool
		wantHour       int // For checking the parsed time.Time object
		wantMinute     int
		wantSecond     int
		wantNanosecond int
	}{
		{
			name:     "Full precision with Z",
			epochStr: "2024-03-16T12:22:20.979840Z",
			wantYear: 2024,
			wantDay:  76.51552060000000168, // Updated to match direct float64 calculation from ns
			wantErr:  false,
			wantHour: 12, wantMinute: 22, wantSecond: 20, wantNanosecond: 979840000,
		},
		{
			name:     "No fractional seconds with Z",
			epochStr: "2025-01-01T00:00:00Z",
			wantYear: 2025,
			wantDay:  1.0,
			wantErr:  false,
			wantHour: 0, wantMinute: 0, wantSecond: 0, wantNanosecond: 0,
		},
		{
			name:     "Implied UTC (no Z, no offset) ISS",
			epochStr: "2025-05-26T13:06:57.824640", // From user example
			wantYear: 2025,
			wantDay:  146.54650259999998692, // Updated
			wantErr:  false,
			wantHour: 13, wantMinute: 6, wantSecond: 57, wantNanosecond: 824640000,
		},
		{
			name:     "Implied UTC (no Z, no offset) CSS",
			epochStr: "2025-05-25T23:00:12.248640",
			wantYear: 2025,
			wantDay:  145.95847509999998692, // Updated to match actual ommEpochToTleEpoch output
			wantErr:  false,
			wantHour: 23, wantMinute: 0, wantSecond: 12, wantNanosecond: 248640000,
		},
		{
			name:     "Implied UTC, 3 fractional digits",
			epochStr: "2023-12-31T23:59:59.123",
			wantYear: 2023,
			wantDay:  365.99998984953703110, // Updated
			wantErr:  false,
			wantHour: 23, wantMinute: 59, wantSecond: 59, wantNanosecond: 123000000,
		},
		{
			name:     "With offset",
			epochStr: "2024-03-16T14:22:20.979840+02:00", // Equivalent to 12:22:20.979840Z
			wantYear: 2024,
			wantDay:  76.51552060000000168, // Updated
			wantErr:  false,
			wantHour: 12, wantMinute: 22, wantSecond: 20, wantNanosecond: 979840000,
		},
		{
			name:     "Invalid epoch string",
			epochStr: "NOT_A_DATE",
			wantErr:  true,
		},
	}

	epochDayTolerance := 1e-15 // Tighten tolerance now that wantDay should match ommEpochToTleEpoch output exactly

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			year, day, parsedTime, err := ommEpochToTleEpoch(tt.epochStr)
			if (err != nil) != tt.wantErr {
				t.Errorf("ommEpochToTleEpoch() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if tt.wantErr {
				return
			}
			if year != tt.wantYear {
				t.Errorf("ommEpochToTleEpoch() year = %v, want %v", year, tt.wantYear)
			}
			if math.Abs(day-tt.wantDay) > epochDayTolerance {
				// Use %0.20f for high precision printing if it fails
				t.Errorf("ommEpochToTleEpoch() day = %.20f, want %.20f (diff %.20e, tol %.1e)", day, tt.wantDay, math.Abs(day-tt.wantDay), epochDayTolerance)
			}
			// Check components of the parsed time.Time object
			baseDateForExpected := time.Date(tt.wantYear, 1, 1, 0, 0, 0, 0, time.UTC)
			expectedDateFromDay := baseDateForExpected.AddDate(0, 0, int(tt.wantDay)-1)

			if parsedTime.Year() != tt.wantYear ||
				parsedTime.Month() != expectedDateFromDay.Month() ||
				parsedTime.Day() != expectedDateFromDay.Day() ||
				parsedTime.Hour() != tt.wantHour ||
				parsedTime.Minute() != tt.wantMinute ||
				parsedTime.Second() != tt.wantSecond ||
				parsedTime.Nanosecond() != tt.wantNanosecond {
				t.Errorf("parsedTime components mismatch: got %04d-%02d-%02dT%02d:%02d:%02d.%09dZ, want %04d-%02d-%02dT%02d:%02d:%02d.%09dZ",
					parsedTime.Year(), parsedTime.Month(), parsedTime.Day(), parsedTime.Hour(), parsedTime.Minute(), parsedTime.Second(), parsedTime.Nanosecond(),
					tt.wantYear, expectedDateFromDay.Month(), expectedDateFromDay.Day(), tt.wantHour, tt.wantMinute, tt.wantSecond, tt.wantNanosecond)
			}
		})
	}
}

func TestComparePosition_TLE_vs_OMM_ISS(t *testing.T) {
	// TLE string for ISS (ZARYA) corresponding to the first OMM entry
	// Epoch: 25146.54650260 -> 2025, Day 146.54650260
	// Day 146 of 2025 (non-leap) is May 26.
	// .54650260 day fraction * 86400 s/day = 47217.82464 seconds
	// 47217.82464 s = 13h * 3600 + 6m * 60 + 57.82464s
	// So, epoch is 2025-05-26 13:06:57.824640 UTC
	issTLEString := `ISS (ZARYA)             
1 25544U 98067A   25146.54650260  .00010397  00000+0  19155-3 0  9999
2 25544  51.6382  54.2937 0002241 147.4648 271.6158 15.49752720511807`

	// 1. Parse TLE from string
	tleFromString, err := ParseTLE(issTLEString)
	if err != nil {
		t.Fatalf("Failed to parse TLE string: %v", err)
	}

	// 2. Parse OMM JSON and convert the first entry (ISS ZARYA) to TLE
	omms, err := ParseOMMs([]byte(ommJsonExample))
	if err != nil {
		t.Fatalf("ParseOMMs failed: %v", err)
	}
	if len(omms) == 0 {
		t.Fatal("No OMMs parsed from example JSON.")
	}
	issOMM := omms[0] // First entry is ISS (ZARYA) with matching epoch
	if issOMM.NoradCatID != 25544 {
		t.Fatalf("Expected OMM for NORAD ID 25544, got %d (%s)", issOMM.NoradCatID, issOMM.ObjectName)
	}

	tleFromOMM, err := issOMM.ToTLE()
	if err != nil {
		t.Fatalf("Failed to convert OMM to TLE: %v", err)
	}

	// 3. Sanity check: Compare parsed TLE objects for basic field equality
	// Check epoch times - allow small difference due to float conversions
	epochTimeFromString := tleFromString.EpochTime()
	epochTimeFromOMM := tleFromOMM.EpochTime()
	epochTimeDiffNs := epochTimeFromString.Sub(epochTimeFromOMM).Nanoseconds()
	maxEpochTimeDiffNs := int64(10) // Allow up to 10ns difference

	if math.Abs(float64(epochTimeDiffNs)) > float64(maxEpochTimeDiffNs) {
		t.Errorf("EpochTime mismatch too large: TLE string parsed to %v, OMM parsed to %v. Diff (ns): %d. Max allowed (ns): %d",
			epochTimeFromString, epochTimeFromOMM, epochTimeDiffNs, maxEpochTimeDiffNs)
		t.Logf("TLE EpochDay: %.17f, OMM EpochDay: %.17f", tleFromString.EpochDay, tleFromOMM.EpochDay)
	}

	// Also check EpochDay directly.
	// TLE string epoch day is parsed as float64. OMM-derived EpochDay is also float64.
	// Their representations might differ at the very last bits.
	epochDayTolerance := 1e-15 // Allow for minimal float64 representation differences
	if math.Abs(tleFromString.EpochDay-tleFromOMM.EpochDay) > epochDayTolerance {
		t.Errorf("EpochDay mismatch too large: TLE string %.20f, OMM conv %.20f. Diff: %.20e",
			tleFromString.EpochDay, tleFromOMM.EpochDay, tleFromString.EpochDay-tleFromOMM.EpochDay)
	}

	// Check a few other key elements to ensure they match before propagation
	if math.Abs(tleFromString.MeanMotion-tleFromOMM.MeanMotion) > 1e-10 { // Tighter check
		t.Errorf("MeanMotion mismatch: TLE string %.15f, OMM %.15f", tleFromString.MeanMotion, tleFromOMM.MeanMotion)
	}
	if math.Abs(tleFromString.Eccentricity-tleFromOMM.Eccentricity) > 1e-10 {
		t.Errorf("Eccentricity mismatch: TLE string %.15f, OMM %.15f", tleFromString.Eccentricity, tleFromOMM.Eccentricity)
	}

	// 4. Propagate both TLEs to their epoch (tsince = 0.0)
	eciFromTLE, errTLE := tleFromString.FindPosition(0.0)
	if errTLE != nil {
		t.Fatalf("Failed to propagate TLE from string: %v", errTLE)
	}

	eciFromOMM, errOMM := tleFromOMM.FindPosition(0.0)
	if errOMM != nil {
		t.Fatalf("Failed to propagate TLE from OMM: %v", errOMM)
	}

	// 5. Compare ECI position vectors
	posToleranceKm := 1e-12 // Kilometers, should be extremely close
	if math.Abs(eciFromTLE.Position.X-eciFromOMM.Position.X) > posToleranceKm {
		t.Errorf("Position.X mismatch: TLE_str=%.15e km, OMM_conv=%.15e km, Diff=%.15e km",
			eciFromTLE.Position.X, eciFromOMM.Position.X, eciFromTLE.Position.X-eciFromOMM.Position.X)
	}
	if math.Abs(eciFromTLE.Position.Y-eciFromOMM.Position.Y) > posToleranceKm {
		t.Errorf("Position.Y mismatch: TLE_str=%.15e km, OMM_conv=%.15e km, Diff=%.15e km",
			eciFromTLE.Position.Y, eciFromOMM.Position.Y, eciFromTLE.Position.Y-eciFromOMM.Position.Y)
	}
	if math.Abs(eciFromTLE.Position.Z-eciFromOMM.Position.Z) > posToleranceKm {
		t.Errorf("Position.Z mismatch: TLE_str=%.15e km, OMM_conv=%.15e km, Diff=%.15e km",
			eciFromTLE.Position.Z, eciFromOMM.Position.Z, eciFromTLE.Position.Z-eciFromOMM.Position.Z)
	}

	// Compare ECI velocity vectors as well
	velToleranceKms := 1e-14 // km/s, extremely tight tolerance
	if math.Abs(eciFromTLE.Velocity.X-eciFromOMM.Velocity.X) > velToleranceKms {
		t.Errorf("Velocity.X mismatch: TLE_str=%.15e km/s, OMM_conv=%.15e km/s, Diff=%.15e km/s",
			eciFromTLE.Velocity.X, eciFromOMM.Velocity.X, eciFromTLE.Velocity.X-eciFromOMM.Velocity.X)
	}
	if math.Abs(eciFromTLE.Velocity.Y-eciFromOMM.Velocity.Y) > velToleranceKms {
		t.Errorf("Velocity.Y mismatch: TLE_str=%.15e km/s, OMM_conv=%.15e km/s, Diff=%.15e km/s",
			eciFromTLE.Velocity.Y, eciFromOMM.Velocity.Y, eciFromTLE.Velocity.Y-eciFromOMM.Velocity.Y)
	}
	if math.Abs(eciFromTLE.Velocity.Z-eciFromOMM.Velocity.Z) > velToleranceKms {
		t.Errorf("Velocity.Z mismatch: TLE_str=%.15e km/s, OMM_conv=%.15e km/s, Diff=%.15e km/s",
			eciFromTLE.Velocity.Z, eciFromOMM.Velocity.Z, eciFromTLE.Velocity.Z-eciFromOMM.Velocity.Z)
	}

	t.Logf("ISS ECI from TLE string: P(%.6f, %.6f, %.6f)km V(%.6f, %.6f, %.6f)km/s @ %s",
		eciFromTLE.Position.X, eciFromTLE.Position.Y, eciFromTLE.Position.Z,
		eciFromTLE.Velocity.X, eciFromTLE.Velocity.Y, eciFromTLE.Velocity.Z, eciFromTLE.DateTime.Format(time.RFC3339Nano))
	t.Logf("ISS ECI from OMM conv : P(%.6f, %.6f, %.6f)km V(%.6f, %.6f, %.6f)km/s @ %s",
		eciFromOMM.Position.X, eciFromOMM.Position.Y, eciFromOMM.Position.Z,
		eciFromOMM.Velocity.X, eciFromOMM.Velocity.Y, eciFromOMM.Velocity.Z, eciFromOMM.DateTime.Format(time.RFC3339Nano))

	// Ensure the DateTime of the ECI objects are also very close (reflecting the epoch time diff)
	eciDateTimeDiffNs := eciFromTLE.DateTime.Sub(eciFromOMM.DateTime).Nanoseconds()
	if math.Abs(float64(eciDateTimeDiffNs)) > float64(maxEpochTimeDiffNs) {
		t.Errorf("ECI.DateTime mismatch too large: TLE_str DT=%v, OMM_conv DT=%v. Diff (ns): %d. Max allowed (ns): %d",
			eciFromTLE.DateTime, eciFromOMM.DateTime, eciDateTimeDiffNs, maxEpochTimeDiffNs)
	}
}
