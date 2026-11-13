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
