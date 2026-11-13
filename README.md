# 01-zz-gps-compare

Python GUI for extracting UBX NAV-PVT messages and comparing multiple GPS logs on an iTOW reference timeline.

## Run

Requires Python 3 with Tkinter.

```sh
python -m pip install -r requirements.txt
python gps_delay_analyzer.py
```

Select one to three GPS logs (`.log`, `.ubx`, or binary UBX files). Empty inputs are skipped. ULog input is optional.

The parser checks UBX checksums, extracts NAV-PVT fields, and filters records without valid 3D fixes. Extracted data is written to `analysis_output/*.pvt.csv` locally.

## Alignment

The first selected GPS provides the reference timeline when ULog is omitted. Other GPS records are matched by nearest iTOW within the configured tolerance in milliseconds. Set tolerance to zero for exact iTOW matches. A source record may match multiple reference timestamps. Original source iTOW values are retained in cursor readouts.

## Plot controls

- Tabs show individual PVT fields; legend names come from input filenames.
- Left-drag on empty plot space: rectangle zoom.
- Middle-drag: pan. Mouse wheel: zoom.
- Space or Home: restore the original view.
- Right-click: alternately place cursors A and B.
- Left-drag near a cursor: move it and update measurements.
- A/B keys: select a cursor; left/right arrows: move it.
- Readouts show interpolated curve values, delta time, delta Y, and nearest matched samples' original iTOW differences.
- Marker size adapts to visible density; all sample markers remain enabled.

## ULog limitations

The current `drone_fusion` support is preliminary. Its yaw units, reference north, coordinate conventions, and mapping from ULog timestamps to GPS iTOW require validation against the flight controller's message definitions. Do not treat its current automatic rotation/time mapping as verified ground truth or a measured GPS latency. GPS-only comparisons do not require ULog.

Raw flight logs and generated CSV files are intentionally excluded from version control.

## RTCM ephemeris inspection

Select the RTCM ephemeris project in the left panel and choose a binary RTCM3 log. `rtcm_ephemeris.py` validates CRC24Q and extracts GPS/GLO/GAL/BDS/QZSS MSM satellite masks and ephemeris satellite IDs.

The chart uses 10-second windows with three series:

- Red line: observed satellites whose ephemeris has not been received in this file by the window end, divided by the observed satellite count.
- Dashed line: observed satellites without an ephemeris message in that particular window, divided by the observed satellite count.
- Background bars on the right axis: distinct observed satellites in each window.

GPS MSM epochs anchor the timeline. Other constellation messages and ephemerides are associated with the preceding GPS epoch in file order. Logs without a usable GPS MSM anchor are rejected. This is an approximate file-order association, not an independent reception timestamp.

Previously received ephemerides count as present; validity, health and receiver caches from before the recording are not evaluated. The dashed series indicates message transmission coverage and does not imply an invalid ephemeris. Window statistics are exported to `analysis_output/<filename>.ephemeris.csv`.
