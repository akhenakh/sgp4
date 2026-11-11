package sgp4

import (
	"math"
	"strings"
	"testing"
)

const ommCSVExample = `OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
STARLINK-5195,2022-136A,2026-10-05T06:00:02.999808,15.33987296,.00014793,53.1614,120.4004,41.5220,36.3784,0,U,54051,999,576,.1629853E-2,.54371E-3,0
NEW-SAT,2026-001A,2026-10-05T06:00:02.999808,13.10000000,.00050000,97.5000,10.0000,20.0000,30.0000,0,U,100953,1,10,.1000000E-3,.10000E-4,0
`

func TestParseOMMsCSV(t *testing.T) {
	omms, err := ParseOMMsCSV(strings.NewReader(ommCSVExample))
	if err != nil {
		t.Fatalf("ParseOMMsCSV: %v", err)
	}
	if len(omms) != 2 {
		t.Fatalf("expected 2 OMMs, got %d", len(omms))
	}
	if omms[0].ObjectName != "STARLINK-5195" || omms[0].NoradCatID != 54051 {
		t.Errorf("first record wrong: %+v", omms[0])
	}
	if math.Abs(omms[0].MeanMotion-15.33987296) > 1e-9 {
		t.Errorf("mean motion: %v", omms[0].MeanMotion)
	}
	if omms[1].NoradCatID != 100953 {
		t.Errorf("6-digit catalog: %d", omms[1].NoradCatID)
	}
}

func TestParseOMMsCSVToTLESixDigit(t *testing.T) {
	omms, err := ParseOMMsCSV(strings.NewReader(ommCSVExample))
	if err != nil {
		t.Fatalf("ParseOMMsCSV: %v", err)
	}
	tle, err := omms[1].ToTLE()
	if err != nil {
		t.Fatalf("ToTLE: %v", err)
	}
	if tle.SatelliteNumber != 100953 {
		t.Fatalf("catalog number: %d", tle.SatelliteNumber)
	}
	if tle.International != "26001A" {
		t.Fatalf("international designator: %q", tle.International)
	}
	if tle.EpochYear != 2026 {
		t.Fatalf("epoch year: %d", tle.EpochYear)
	}
}

func TestParseOMMsReaderDetectsFormat(t *testing.T) {
	csv, err := ParseOMMsReader(strings.NewReader(ommCSVExample))
	if err != nil || len(csv) != 2 {
		t.Fatalf("CSV detection: n=%d err=%v", len(csv), err)
	}
	json, err := ParseOMMsReader(strings.NewReader(ommJsonExample))
	if err != nil || len(json) != 3 {
		t.Fatalf("JSON detection: n=%d err=%v", len(json), err)
	}
}

func TestParseOMMsReaderSingleObject(t *testing.T) {
	one := `{"OBJECT_NAME":"SINGLE","NORAD_CAT_ID":25544,"EPOCH":"2025-05-18T08:53:29.535936","MEAN_MOTION":15.49587957,"ECCENTRICITY":0.0002558,"INCLINATION":51.6369,"RA_OF_ASC_NODE":94.7823,"ARG_OF_PERICENTER":120.7586,"MEAN_ANOMALY":15.784}`
	omms, err := ParseOMMsReader(strings.NewReader(one))
	if err != nil {
		t.Fatalf("ParseOMMsReader: %v", err)
	}
	if len(omms) != 1 || omms[0].NoradCatID != 25544 || omms[0].ObjectName != "SINGLE" {
		t.Fatalf("single object not decoded: %+v", omms)
	}
}

// TestParseOMMsCSVRejectsMalformedValues guards the numeric fields: an
// unparsable value used to become 0 and reach the propagator as a valid set.
func TestParseOMMsCSVRejectsMalformedValues(t *testing.T) {
	tests := []struct {
		name string
		row  string
	}{
		{"bad mean motion", "BAD,2026-001A,2026-10-05T06:00:02.999808,abc,.00014793,53.1614,120.4004,41.5220,36.3784,0,U,100953,999,576,.1629853E-2,.54371E-3,0"},
		{"bad catalog number", "BAD,2026-001A,2026-10-05T06:00:02.999808,15.3,.00014793,53.1614,120.4004,41.5220,36.3784,0,U,not-a-number,999,576,.1629853E-2,.54371E-3,0"},
		{"bad bstar", "BAD,2026-001A,2026-10-05T06:00:02.999808,15.3,.00014793,53.1614,120.4004,41.5220,36.3784,0,U,100953,999,576,oops,.54371E-3,0"},
		{"empty required field", "BAD,2026-001A,2026-10-05T06:00:02.999808,,.00014793,53.1614,120.4004,41.5220,36.3784,0,U,100953,999,576,.1629853E-2,.54371E-3,0"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			data := strings.SplitN(ommCSVExample, "\n", 2)[0] + "\n" + tt.row + "\n"
			if _, err := ParseOMMsCSV(strings.NewReader(data)); err == nil {
				t.Fatal("ParseOMMsCSV accepted a malformed record")
			}
		})
	}
}

func TestParseOMMsCSVRequiresColumns(t *testing.T) {
	data := "OBJECT_NAME,NORAD_CAT_ID\nX,25544\n"
	_, err := ParseOMMsCSV(strings.NewReader(data))
	if err == nil || !strings.Contains(err.Error(), "EPOCH") {
		t.Fatalf("error = %v, want a missing EPOCH column error", err)
	}
}

func TestParseOMMsCSVOptionalFieldsMayBeAbsent(t *testing.T) {
	data := "EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,NORAD_CAT_ID\n" +
		"2026-10-05T06:00:02.999808,15.3,.00014793,53.1614,120.4004,41.5220,36.3784,25544\n"
	omms, err := ParseOMMsCSV(strings.NewReader(data))
	if err != nil {
		t.Fatalf("ParseOMMsCSV: %v", err)
	}
	if len(omms) != 1 || omms[0].BStar != 0 || omms[0].NoradCatID != 25544 {
		t.Fatalf("record = %+v", omms)
	}
}

func TestToTLEAllowsMissingObjectID(t *testing.T) {
	o := OMM{
		ObjectName:  "NO-ID",
		ObjectID:    "",
		EpochStr:    "2026-10-05T06:00:02.999808",
		MeanMotion:  15.0,
		Inclination: 51.6,
	}
	tle, err := o.ToTLE()
	if err != nil {
		t.Fatalf("ToTLE with missing ObjectID should succeed: %v", err)
	}
	if tle.SatelliteNumber != 0 {
		t.Fatalf("satellite number: %d", tle.SatelliteNumber)
	}
	if tle.International != "" {
		t.Fatalf("international should be blank, got %q", tle.International)
	}
}

func TestTLEToOMMRoundTrip(t *testing.T) {
	const tleStr = `ISS (ZARYA)
1 25544U 98067A   25247.10182809  .00011777  00000-0  21333-3 0  9997
2 25544  51.6327 275.9345 0004179 299.5263  60.5309 15.50088696528307`
	tle, err := ParseTLE(tleStr)
	if err != nil {
		t.Fatalf("parse TLE: %v", err)
	}
	omm := tle.ToOMM()
	if omm.ObjectID != "1998-067A" {
		t.Fatalf("OBJECT_ID: %q", omm.ObjectID)
	}
	back, err := omm.ToTLE()
	if err != nil {
		t.Fatalf("OMM.ToTLE: %v", err)
	}
	if back.SatelliteNumber != tle.SatelliteNumber ||
		back.Inclination != tle.Inclination ||
		back.MeanMotion != tle.MeanMotion ||
		back.RightAscension != tle.RightAscension {
		t.Fatalf("round trip mismatch:\n want %+v\n got  %+v", tle, back)
	}
	if back.International != "98067A" {
		t.Fatalf("international designator round trip: %q", back.International)
	}
}
