# Local MacBook conversion test

Place `local_convert.py`, `run_local.sh`, and `requirements-local.txt` together in your repository folder. This runner uses the new bounded-memory converter. The existing Condor wrapper and converter are not used by this local command.

## One-time setup

Use Python 3.10 or newer. From the repository folder:

```bash
python3 --version
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-local.txt
```

## First test: five input events from one file

Replace the input path with the actual location on your Mac. Keep the quotes if the path contains spaces.

```bash
bash run_local.sh "/Users/sarahwaldych/Documents/pixel_clusters_d16401.out.gz" \
  --output-dir local_test_d16401 \
  --max-events 5
```

The default is five input events, including events excluded by the `z-entry == 100` filter. Fewer than five output rows can therefore be correct. The three files are written under `local_test_d16401/unflipped/`:

- `labels_d16401.parquet`
- `recon2D_d16401.parquet`
- `recon3D_d16401.parquet`

The gzip input is left intact and never decompressed to a separate disk file. The output directory must not already contain these filenames. An interrupted or failed conversion does not publish completed-looking final filenames; a forcibly killed process may leave a hidden temporary directory.

## Convert the entire one-file input

After inspecting the small test:

```bash
bash run_local.sh "/Users/sarahwaldych/Documents/pixel_clusters_d16401.out.gz" \
  --output-dir local_full_d16401 \
  --max-events 0
```

This still processes one event at a time. Full conversion can take substantial time and disk space. Memory does not accumulate all event payloads, although Parquet row-group metadata grows with the number of events. No MacBook peak-memory measurement is claimed.

## Read back one event

Run with `.venv/bin/python` (or in a notebook using that environment):

```python
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq

root = Path('local_test_d16401/unflipped')
labels = pq.read_table(root / 'labels_d16401.parquet')
print(labels.to_pylist())

pf = pq.ParquetFile(root / 'recon3D_d16401.parquet')
print('Saved events:', pf.metadata.num_rows)
if pf.metadata.num_rows:
    first = pf.read_row_group(0, columns=['charge'])
    charge = first.column('charge').chunk(0).values.to_numpy()
    charge = charge.reshape(400, 48, 192)
    print(charge.shape)
```

## Output compatibility

This local format stores one row per retained event. `recon2D` and `recon3D` have one `charge` list column (float64) instead of one column per charge value. The three files retain identical event ordering. The labels retain the original named physics and boundary fields; the implicit pandas index is not stored.

The flattened charge order, float64 precision, final-time-slice choice for recon2D, 100-micrometre sensor calculation, and `z-entry == 100` selection are preserved. The legacy boundary routine counts corners twice; this behavior is preserved for comparison.

The `(400, 48, 192)` readback convention matches the old Python converter. The C writer actually writes 192 rows of 48 numbers per slice. This local change preserves the existing converter's reshape convention rather than reinterpreting detector axes; physical axis/boundary interpretation needs a separate review before analysis.

Existing readers that expect `table.to_pandas().to_numpy().reshape(...)` on the charge files must use the list-column readback above. Limited runs stop after the requested complete events and do not validate the remaining gzip data. Full runs validate the gzip stream through end-of-file. The parser reports missing truth values, malformed numeric rows, incomplete slices/events, and nonfinite values.

## Changes made

- Added `local_convert.py`: incremental `.gz`/`.out` input; one preallocated float64 event buffer; event-wise Parquet writes; progress messages; explicit input validation; event limit; output collision checks; temporary outputs until successful conversion.
- Added `run_local.sh`: resolves its own directory and launches the local virtual environment without CERN software or network transfers. `PIXELAV_PYTHON` can override the Python executable with an absolute path.
- Added `requirements-local.txt`: NumPy and PyArrow dependencies. No pandas dependency is needed for this local path.
- Added this setup, execution, output-format, and readback guide.

The original Condor commands remain available separately. Do not use `datagen.sh` for this MacBook test.
