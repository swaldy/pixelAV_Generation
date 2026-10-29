#!/usr/bin/env python3
"""Validate PixelAV 48x192 parquet outputs and build the inputs for size-plots.
 
For every index in INPUT_DIR this reads the triplet written by
pixelAV_Generation/local_convert.py:
 
    labels_dNNNNN.parquet   recon2D_dNNNNN.parquet   recon3D_dNNNNN.parquet
 
and then
 
  1. validates it: file set, SUCCESS markers, schemas, row counts, label
     consistency, charge sanity, pixel layout vs. truth position,
     recon3D last time slice == recon2D, duplicate tracks across files;
  2. computes per-cluster x-size / y-size (the "recon0D" quantities that
     smart-pixels-ml/plotting/size-plots.ipynb reads from CSV), charge
     centroids and edge flags on the correct pixel layout;
  3. writes OUTPUT_DIR/labels_with_sizes.parquet and
     OUTPUT_DIR/validation_report.txt.
 
Pixel layout
------------
For each time slice the PixelAV driver prints 192 lines (y, 12.5 um pitch)
of 48 values (x, 50 um pitch), so each flat 'charge' list is C-ordered as
(y, x) = (192, 48), and (400, 192, 48) for recon3D. local_convert.py labels
the arrays "48,192" and computes the stored original_atEdge columns with that
transposed reshape; this script recomputes the edge flags on (192, 48) and
also checks the layout against the truth impact position.
 
Usage (lxplus or SWAN, where /eos/user is mounted):
 
    python3 prep_size_inputs.py \\
        /eos/user/s/swaldych/smart_pix/pixelav_48x192_10ps/parquet_streaming/<run_id>/unflipped \\
        --first 16401 --last 16440 --output-dir prepared_48x192
 
Each finished file is saved under OUTPUT_DIR/per_file/, so if a run is
interrupted (Ctrl-C), running the same command again skips the files already
done (use --no-cache to redo everything). --check-3d 0 skips the recon3D
spot check, which is the only step that touches the large recon3D files.
 
Exit status is 1 if any check FAILs, 0 otherwise.
"""
import argparse
import json
from pathlib import Path
import re
import sys
import time
 
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
 
# ---------------------------------------------------------------- geometry
NY, NX, NT = 192, 48, 400                  # rows = y (fine), cols = x (coarse)
PIXELS = NY * NX
PITCH_X, PITCH_Y, THICK = 50.0, 12.5, 100.0  # um
XOFF = (NX + 1) // 2 - 0.5                 # 23.5, as in detect() in ppixelav2.c
YOFF = (NY + 1) // 2 - 0.5                 # 95.5
# Impact points are sampled uniformly over the full matrix at the midplane:
#   x_mid = 0.5*xsize + NX*xsize*(r - 0.5)  ->  [-1175, 1225] um
X_MID_RANGE = (0.5 * PITCH_X - NX * PITCH_X / 2, 0.5 * PITCH_X + NX * PITCH_X / 2)
Y_MID_RANGE = (0.5 * PITCH_Y - NY * PITCH_Y / 2, 0.5 * PITCH_Y + NY * PITCH_Y / 2)
COTBETA_CUT = 2.0                          # filter_tracklist.py cut on |cotBeta|
COT_MAX = 10.0                             # PixelAV driver skips |cot| > 10
 
TRUTH = ['x-entry', 'y-entry', 'z-entry', 'n_x', 'n_y', 'n_z',
         'number_eh_pairs', 'y-local', 'pt']
DERIVED = ['cotAlpha', 'cotBeta', 'y-midplane', 'x-midplane']
STORED_EDGE = ['original_atEdge', 'chargeOriginal_atEdge', 'nPixOriginalAbove1e_atEdge']
LABEL_COLUMNS = TRUTH + DERIVED + STORED_EDGE
 
BORDER = np.zeros((NY, NX), dtype=bool)
BORDER[0, :] = BORDER[-1, :] = BORDER[:, 0] = BORDER[:, -1] = True
EDGE_THRESHOLD = 1.0                       # e-, same as the converter
NONZERO = 1e-6                             # charges are printed as %2.1f, so any
                                           # real nonzero profile bin is >= 0.1
 
 
class Report:
    def __init__(self):
        self.lines = []
        self.issues = []          # (index, level, message)
 
    def add(self, index, level, message):
        self.issues.append((index, level, message))
        print(f'  [{level}] d{index}: {message}' if index else f'  [{level}] {message}', flush=True)
 
    def n(self, level):
        return sum(1 for _, lvl, _ in self.issues if lvl == level)
 
 
# ---------------------------------------------------------------- helpers
def find_indices(input_dir, first, last):
    found = sorted(int(m[1]) for p in input_dir.glob('labels_d*.parquet')
                   if (m := re.fullmatch(r'labels_d(\d+)\.parquet', p.name)))
    if first is not None or last is not None:
        lo = first if first is not None else (found[0] if found else 0)
        hi = last if last is not None else (found[-1] if found else -1)
        return list(range(lo, hi + 1))
    if found:
        return list(range(found[0], found[-1] + 1))   # so gaps are reported
    return []
 
 
CACHE_VERSION = 3      # bump when the checks or output columns change
 
 
def cache_key(args, index):
    """Identifies one file's result: options + size/mtime of every input file."""
    sig = {}
    for name in (f'labels_d{index}.parquet', f'recon2D_d{index}.parquet',
                 f'recon3D_d{index}.parquet', f'SUCCESS_d{index}.txt'):
        p = args.input_dir / name
        sig[name] = [p.stat().st_size, int(p.stat().st_mtime)] if p.is_file() else None
    return {'version': CACHE_VERSION, 'threshold': args.threshold, 'check_3d': args.check_3d,
            'files': sig}
 
 
def read_success(input_dir, index):
    path = input_dir / f'SUCCESS_d{index}.txt'
    if not path.is_file():
        return None
    out = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out
 
 
def charge_batches(pf, expected_len, batch_size):
    """Yield (n, expected_len) float64 arrays; raise ValueError on bad rows."""
    for batch in pf.iter_batches(batch_size=batch_size, columns=['charge']):
        col = batch.column(0)
        if col.null_count:
            raise ValueError(f'{col.null_count} null charge lists')
        lengths = pc.list_value_length(col).to_numpy(zero_copy_only=False)
        bad = lengths != expected_len
        if bad.any():
            raise ValueError(f'{bad.sum()} rows with list length != {expected_len} '
                             f'(e.g. {lengths[bad][0]})')
        flat = col.flatten().to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        yield flat.reshape(len(lengths), expected_len)
 
 
def size_from_profile(profile):
    return (np.abs(profile) > NONZERO).sum(axis=1)
 
 
def analyse_recon2d(batch, threshold):
    """Per-event quantities from a (n, 9216) batch of last-time-slice charge."""
    n = batch.shape[0]
    img = batch.reshape(n, NY, NX)                   # (y, x): the correct layout
    imgt = np.where(img >= threshold, img, 0.0) if threshold > 0 else img
    out = {
        'x-size': size_from_profile(imgt.sum(axis=1)),   # project onto x (48 cols)
        'y-size': size_from_profile(imgt.sum(axis=2)),   # project onto y (192 rows)
        'charge_total': img.sum(axis=(1, 2)),
    }
    # Charge centroid (positive charge only: induced charge can be negative).
    w = np.clip(img, 0.0, None)
    wsum = w.sum(axis=(1, 2))
    safe = np.where(wsum > 0, wsum, np.nan)
    cx = w.sum(axis=1) @ (np.arange(NX) + 0.5) / safe
    cy = w.sum(axis=2) @ (np.arange(NY) + 0.5) / safe
    out['x-centroid'] = (cx - XOFF) * PITCH_X       # um, same frame as x-midplane
    out['y-centroid'] = (cy - YOFF) * PITCH_Y
    # Same centroid if the flat list were (wrongly) read as (48 x-rows, 192 y-cols).
    wt = np.clip(batch.reshape(n, NX, NY), 0.0, None)
    out['_x_centroid_T'] = (wt.sum(axis=2) @ (np.arange(NX) + 0.5) / safe - XOFF) * PITCH_X
    out['_y_centroid_T'] = (wt.sum(axis=1) @ (np.arange(NY) + 0.5) / safe - YOFF) * PITCH_Y
    # Edge flags on the correct layout (each border pixel counted once).
    border = img[:, BORDER]
    above = np.abs(border) > EDGE_THRESHOLD
    out['atEdge'] = above.any(axis=1)
    out['nPixAbove1e_atEdge'] = above.sum(axis=1).astype(np.int64)
    out['charge_atEdge'] = np.where(out['atEdge'], border.sum(axis=1), 0.0)
    # Reproduce the converter's stored flags (it duplicates corner pixels) for both
    # SHAPE settings: 'old' = (48, 192) transposed, 'fixed' = (192, 48).
    for tag, shape in (('old', (n, NX, NY)), ('fixed', (n, NY, NX))):
        m = batch.reshape(shape)
        edges = np.concatenate([m[:, 0, :], m[:, -1, :], m[:, :, 0], m[:, :, -1]], axis=1)
        a = np.abs(edges) > EDGE_THRESHOLD
        out[f'_{tag}_atEdge'] = a.any(axis=1)
        out[f'_{tag}_nPix'] = a.sum(axis=1).astype(np.int64)
        out[f'_{tag}_charge'] = np.where(a.any(axis=1), edges.sum(axis=1), 0.0)
    out['_finite'] = np.isfinite(batch).all(axis=1)
    return out
 
 
def check_labels(index, df, rep):
    missing = [c for c in LABEL_COLUMNS if c not in df.columns]
    if missing:
        rep.add(index, 'FAIL', f'labels missing columns {missing}')
        return None
    extra = [c for c in df.columns if c not in LABEL_COLUMNS]
    if extra:
        rep.add(index, 'INFO', f'labels has extra columns {extra}')
    num = df[TRUTH + DERIVED].to_numpy(dtype=np.float64)
    if not np.isfinite(num).all():
        rep.add(index, 'FAIL', f'{(~np.isfinite(num)).any(axis=1).sum()} label rows with NaN/inf')
    if not (df['z-entry'] == 100).all():
        rep.add(index, 'FAIL', f'{(df["z-entry"] != 100).sum()} rows with z-entry != 100 '
                               '(converter should keep only unflipped tracks)')
    if not (df['n_z'] < 0).all():
        rep.add(index, 'WARN', f'{(df["n_z"] >= 0).sum()} rows with n_z >= 0 '
                               '(unflipped tracks should travel toward -z)')
    with np.errstate(divide='ignore', invalid='ignore'):
        ok_a = np.allclose(df['cotAlpha'], df['n_x'] / df['n_z'], rtol=1e-9, atol=0)
        ok_b = np.allclose(df['cotBeta'], df['n_y'] / df['n_z'], rtol=1e-9, atol=0)
    if not (ok_a and ok_b):
        rep.add(index, 'FAIL', 'cotAlpha/cotBeta != n_x/n_z, n_y/n_z')
    xm = df['x-entry'] + df['cotAlpha'] * (THICK / 2 - df['z-entry'])
    ym = df['y-entry'] + df['cotBeta'] * (THICK / 2 - df['z-entry'])
    if not (np.allclose(xm, df['x-midplane'], atol=1e-6) and np.allclose(ym, df['y-midplane'], atol=1e-6)):
        rep.add(index, 'FAIL', 'x/y-midplane inconsistent with entry point and slopes')
    tol = 0.1
    out_x = ((df['x-midplane'] < X_MID_RANGE[0] - tol) | (df['x-midplane'] > X_MID_RANGE[1] + tol)).sum()
    out_y = ((df['y-midplane'] < Y_MID_RANGE[0] - tol) | (df['y-midplane'] > Y_MID_RANGE[1] + tol)).sum()
    if out_x or out_y:
        rep.add(index, 'WARN', f'{out_x} x-midplane / {out_y} y-midplane values outside the '
                               f'sampled matrix {X_MID_RANGE} / {Y_MID_RANGE} um')
    nb = (df['cotBeta'].abs() >= COTBETA_CUT).sum()
    if nb:
        rep.add(index, 'WARN', f'{nb} rows with |cotBeta| >= {COTBETA_CUT} '
                               '(track list was filtered to |cotBeta| < 2: check alpha/beta naming)')
    nc = ((df['cotAlpha'].abs() > COT_MAX) | (df['cotBeta'].abs() > COT_MAX)).sum()
    if nc:
        rep.add(index, 'FAIL', f'{nc} rows with |cot| > {COT_MAX} (driver should have skipped them)')
    extra = {}
    neh = df['number_eh_pairs'].to_numpy(np.float64)
    zero = int((neh <= 0).sum())
    if zero:
        rep.add(index, 'FAIL' if zero > 0.01 * len(df) else 'WARN',
                f'{zero} of {len(df)} tracks ({zero / len(df):.0%}) made no e-h pairs. A track crossing '
                f'{THICK:g} um of Si should make thousands. Most likely cause: PixelAV could not find '
                'SIRUTH.SPR in its working directory (look for "cant find SIRUTH.SPR" in the '
                'pixelav_block_NN.log); it then runs with a flat placeholder cross-section')
    path = THICK * np.sqrt(1 + df['cotAlpha'].to_numpy() ** 2 + df['cotBeta'].to_numpy() ** 2)
    if (neh > 0).any():
        extra['med_eh_per_um'] = float(np.median(neh[neh > 0] / path[neh > 0]))
        if not 40 <= extra['med_eh_per_um'] <= 200:
            rep.add(index, 'WARN', f'median e-h pairs per um of track = {extra["med_eh_per_um"]:.0f} '
                                   '(silicon MIP: roughly 70-110)')
    extra['frac_no_eh'] = zero / len(df)
    if (df['pt'] == 0).any():
        rep.add(index, 'WARN', f'{(df["pt"] == 0).sum()} rows with pt == 0 (charge sign lost)')
    return extra
 
 
MAX_RG_EVENTS = 50     # never read a recon3D row group holding more events than this
 
 
def check_recon3d(index, path, samples, rep):
    """samples: {row: recon2D row array}. Reads only the row group holding each
    sampled event (local_convert.py writes one event per row group)."""
    pf = pq.ParquetFile(path)
    md = pf.metadata
    starts = np.cumsum([0] + [md.row_group(g).num_rows for g in range(md.num_row_groups)])
    results = []
    for r in sorted(samples):
        g = int(np.searchsorted(starts, r, side='right') - 1)
        size = md.row_group(g).num_rows
        if size > MAX_RG_EVENTS:
            rep.add(index, 'WARN', f'recon3D row {r} is stored in a row group of {size} events '
                                   f'(~{size * NT * PIXELS * 8 / 1e9:.0f} GB in memory); spot check skipped. '
                                   'This file was not written one event per row group like local_convert.py does')
            continue
        print(f'    recon3D: reading event {r} (row group {g}/{md.num_row_groups}) ...', flush=True)
        col = pf.read_row_group(g, columns=['charge']).column(0)
        if isinstance(col, pa.ChunkedArray):
            col = pa.concat_arrays(col.chunks)
        col = col.slice(int(r - starts[g]), 1)
        n = pc.list_value_length(col).to_numpy(zero_copy_only=False)
        if col.null_count or n.size != 1 or n[0] != NT * PIXELS:
            rep.add(index, 'FAIL', f'recon3D row {r}: list length {n.tolist()} != {NT * PIXELS}')
            continue
        arr = col.flatten().to_numpy(zero_copy_only=False).reshape(NT, PIXELS)
        if not np.isfinite(arr).all():
            rep.add(index, 'FAIL', f'recon3D row {r}: NaN/inf values')
            continue
        if not np.array_equal(arr[-1], samples[r]):
            rep.add(index, 'FAIL', f'recon3D row {r}: last time slice differs from recon2D')
        qt = arr.sum(axis=1)
        final = qt[-1] if qt[-1] != 0 else np.nan
        results.append({'row': r, 'q_1ns_over_final': qt[99] / final,
                        'q_first_over_final': qt[0] / final})
    return results
 
 
# ---------------------------------------------------------------- per file
def process_index(index, input_dir, args, rep):
    paths = {k: input_dir / f'{k}_d{index}.parquet' for k in ('labels', 'recon2D', 'recon3D')}
    missing = [k for k, p in paths.items() if not p.is_file()]
    if missing:
        rep.add(index, 'FAIL', f'missing files: {", ".join(missing)}')
        return None, None
    success = read_success(input_dir, index)
    if success is None:
        rep.add(index, 'WARN', 'no SUCCESS marker (uploads are not atomic; job may not have finished)')
    elif success.get('max_events', '0') != '0':
        rep.add(index, 'FAIL', f'SUCCESS marker says max_events={success["max_events"]}: '
                               'this is an event-limited test conversion, not the full file')
 
    pfs = {k: pq.ParquetFile(p) for k, p in paths.items()}
    nrows = {k: pf.metadata.num_rows for k, pf in pfs.items()}
    if len(set(nrows.values())) != 1:
        rep.add(index, 'FAIL', f'row counts differ: {nrows}')
        return None, None
    n = nrows['labels']
    if n == 0:
        rep.add(index, 'FAIL', 'files are empty')
        return None, None
    for k in ('recon2D', 'recon3D'):
        schema = pfs[k].schema_arrow
        if schema.names != ['charge'] or str(schema.field('charge').type) not in (
                'list<item: double>', 'list<element: double>'):
            rep.add(index, 'FAIL', f'{k} schema is {schema} (expected one list<double> column "charge")')
            return None, None
    meta = pfs['recon2D'].schema_arrow.metadata or {}
    stats = {'index': index, 'rows': n,
             'shape_metadata': meta.get(b'shape', b'?').decode(),
             'max_events': None if success is None else success.get('max_events')}
 
    t0 = time.time()
    progress = lambda msg: print(f'    {msg} [{time.time() - t0:.0f} s]', flush=True)
    labels = pd.read_parquet(paths['labels'])
    extra = check_labels(index, labels, rep)
    if extra is None:
        return None, None
    stats.update(extra)
    progress(f'labels: {n} rows checked')
 
    # recon2D pass
    k3d = min(max(args.check_3d, 0), n)
    want3d = sorted(set(np.linspace(0, n - 1, k3d).round().astype(int).tolist())) if k3d else []
    samples, parts, offset = {}, [], 0
    try:
        for batch in charge_batches(pfs['recon2D'], PIXELS, args.batch_size):
            for r in want3d:
                if offset <= r < offset + len(batch):
                    samples[r] = batch[r - offset].copy()
            parts.append(analyse_recon2d(batch, args.threshold))
            offset += len(batch)
    except ValueError as err:
        rep.add(index, 'FAIL', f'recon2D: {err}')
        return None, None
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    progress('recon2D: sizes computed')
 
    bad = (~d['_finite']).sum()
    if bad:
        rep.add(index, 'FAIL', f'{bad} recon2D rows with NaN/inf')
    empty = (d['charge_total'] <= 0).sum()
    if empty:
        rep.add(index, 'WARN', f'{empty} clusters with no positive charge in the matrix')
 
    # Stored edge flags: which converter SHAPE produced them?
    st_at = labels['original_atEdge'].to_numpy(bool)
    st_n = labels['nPixOriginalAbove1e_atEdge'].to_numpy(np.int64)
    st_q = labels['chargeOriginal_atEdge'].to_numpy(np.float64)
    layout = 'neither'
    for tag in ('fixed', 'old'):
        if (np.array_equal(st_at, d[f'_{tag}_atEdge']) and np.array_equal(st_n, d[f'_{tag}_nPix'])
                and np.allclose(st_q, d[f'_{tag}_charge'], atol=1e-6)):
            layout = tag
            break
    disagree = (st_at != d['atEdge']).sum()
    stats.update(stored_edge_layout=layout,
                 atEdge_frac=d['atEdge'].mean(), stored_atEdge_frac=st_at.mean(),
                 edge_flag_disagree=int(disagree))
    if layout == 'neither':
        rep.add(index, 'WARN', 'stored original_atEdge columns cannot be reproduced from recon2D with '
                               "the converter's logic: labels and recon2D rows may not belong together")
 
    # Layout check against truth (contained clusters with charge).
    sel = (~d['atEdge']) & (d['charge_total'] > 0)
    if sel.sum() >= 5:
        dx = np.abs(d['x-centroid'][sel] - labels['x-midplane'].to_numpy()[sel]) / PITCH_X
        dy = np.abs(d['y-centroid'][sel] - labels['y-midplane'].to_numpy()[sel]) / PITCH_Y
        dxT = np.abs(d['_x_centroid_T'][sel] - labels['x-midplane'].to_numpy()[sel]) / PITCH_X
        dyT = np.abs(d['_y_centroid_T'][sel] - labels['y-midplane'].to_numpy()[sel]) / PITCH_Y
        stats.update(med_dx_pix=np.median(dx), med_dy_pix=np.median(dy),
                     med_dx_pix_T=np.median(dxT), med_dy_pix_T=np.median(dyT))
        res = stats['med_dx_pix'] + stats['med_dy_pix']
        res_T = stats['med_dx_pix_T'] + stats['med_dy_pix_T']
        on_track = stats['med_dx_pix'] < 1.0 and stats['med_dy_pix'] < 3.0
        if res_T <= 2 * res:
            rep.add(index, 'FAIL',
                    f'pixel layout: reading the data as (192 y, 48 x) does not clearly beat (48, 192) '
                    f'(median |d| {res:.2f} vs {res_T:.2f} px): x/y-size may be swapped')
        elif not on_track:
            rep.add(index, 'FAIL',
                    f'charge is not centred where the track crosses the sensor midplane (median '
                    f'|dx|={stats["med_dx_pix"]:.2f} px, |dy|={stats["med_dy_pix"]:.2f} px; expect '
                    f'< 1 px and < 3 px). The pixel layout is fine ((48,192) reading is worse: '
                    f'{res_T:.1f} px); the deposited charge itself looks wrong')
        # Size vs angle sanity: y-size ~ 8|cotBeta| + O(1), x-size ~ 2|cotAlpha| + O(1).
        ry = d['y-size'][sel] - THICK / PITCH_Y * labels['cotBeta'].abs().to_numpy()[sel]
        rx = d['x-size'][sel] - THICK / PITCH_X * labels['cotAlpha'].abs().to_numpy()[sel]
        stats.update(med_ysize_minus_geom=np.median(ry), med_xsize_minus_geom=np.median(rx),
                     contained=int(sel.sum()))
        with np.errstate(invalid='ignore', divide='ignore'):
            frac = d['charge_total'][sel] / labels['number_eh_pairs'].to_numpy()[sel]
        stats['med_charge_over_neh'] = np.median(frac)
    else:
        rep.add(index, 'WARN', f'only {sel.sum()} contained clusters: layout/size checks skipped')
 
    if want3d:
        try:
            res = check_recon3d(index, paths['recon3D'], samples, rep)
            if res:
                stats['q_1ns_over_final'] = np.median([r['q_1ns_over_final'] for r in res])
        except Exception as err:          # report, don't abort the whole run
            rep.add(index, 'FAIL', f'recon3D spot check failed: {type(err).__name__}: {err}')
        progress('recon3D: spot check done')
 
    out = labels.copy()
    out.insert(0, 'row', np.arange(n, dtype=np.int64))
    out.insert(0, 'source_index', np.int64(index))
    for k in ('x-size', 'y-size', 'charge_total', 'x-centroid', 'y-centroid',
              'atEdge', 'nPixAbove1e_atEdge', 'charge_atEdge'):
        out[k] = d[k]
    out['x-size'] = out['x-size'].astype(np.int64)
    out['y-size'] = out['y-size'].astype(np.int64)
    return out, stats
 
 
# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input_dir', type=Path, help='directory with labels/recon2D/recon3D_dNNNNN.parquet')
    ap.add_argument('--output-dir', type=Path, default=Path('prepared_48x192'))
    ap.add_argument('--first', type=int, help='first index expected (e.g. 16401)')
    ap.add_argument('--last', type=int, help='last index expected (e.g. 16440)')
    ap.add_argument('--threshold', type=float, default=0.0,
                    help='pixel charge threshold [e-] before counting sizes '
                         '(default 0: count rows/cols with nonzero net charge, as in the paper)')
    ap.add_argument('--check-3d', type=int, default=3,
                    help='recon3D events to spot-check per file, evenly spaced (default 3; 0 to skip)')
    ap.add_argument('--batch-size', type=int, default=500, help='recon2D rows per read batch')
    ap.add_argument('--no-cache', action='store_true',
                    help='redo every file; by default finished files are reused from '
                         'OUTPUT_DIR/per_file/ if their inputs and options are unchanged')
    args = ap.parse_args(argv)
 
    if not args.input_dir.is_dir():
        ap.error(f'not a directory: {args.input_dir}')
    indices = find_indices(args.input_dir, args.first, args.last)
    if not indices:
        ap.error(f'no labels_d*.parquet files in {args.input_dir}')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output_dir / 'per_file'
    cache_dir.mkdir(exist_ok=True)
 
    rep = Report()
    frames, stats = [], []
    print(f'Checking {len(indices)} indices d{indices[0]}..d{indices[-1]} in {args.input_dir}', flush=True)
    for index in indices:
        key = cache_key(args, index)
        cj, cp = cache_dir / f'd{index}.json', cache_dir / f'd{index}.parquet'
        if not args.no_cache and cj.is_file():
            cached = json.loads(cj.read_text())
            if cached.get('key') == key and (cached['stats'] is None or cp.is_file()):
                print(f'd{index} (reused from {cache_dir})', flush=True)
                for level, message in cached['issues']:
                    rep.add(index, level, message)
                if cached['stats'] is not None:
                    frames.append(pd.read_parquet(cp))
                    stats.append(cached['stats'])
                continue
        print(f'd{index}', flush=True)
        before = len(rep.issues)
        df, st = process_index(index, args.input_dir, args, rep)
        if df is not None:
            frames.append(df)
            stats.append(st)
            df.to_parquet(cp, index=False)
        cj.write_text(json.dumps({'key': key, 'stats': st,     # written last: marks the file done
                                  'issues': [[lvl, msg] for _, lvl, msg in rep.issues[before:]]},
                                 default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
 
    if not frames:
        print('No usable files.', flush=True)
        return 1
    combined = pd.concat(frames, ignore_index=True)
    dup = combined.duplicated(subset=TRUTH, keep=False)
    if dup.any():
        rep.add(None, 'WARN', f'{dup.sum()} rows share identical truth values across/within files '
                              '(same track simulated twice?)')
    st = pd.DataFrame(stats).set_index('index')
    if len(st) > 3:
        med = st['rows'].median()
        low = st.index[st['rows'] < 0.8 * med].tolist()
        if low:
            rep.add(None, 'WARN', f'files with < 80% of the median row count ({med:.0f}): {low}')
    if 'med_dx_pix' in st and ((st['med_dx_pix_T'] + st['med_dy_pix_T'])
                               > 2 * (st['med_dx_pix'] + st['med_dy_pix'])).all():
        rep.add(None, 'INFO', 'pixel layout confirmed: charge centroids are far closer to the truth '
                              'position when the flat list is read as (192 y-rows, 48 x-cols)')
    if 'frac_no_eh' in st and (st['frac_no_eh'] > 0.01).any():
        bad = st.index[st['frac_no_eh'] > 0.01].tolist()
        rep.add(None, 'INFO', f'{len(bad)} of {len(st)} files have tracks with zero e-h pairs: '
                              f'{bad[0]}..{bad[-1]}' if len(bad) > 1 else f'1 file has tracks with zero e-h pairs: {bad}')
    old = st['stored_edge_layout'] == 'old'
    if old.any():
        n_changed = int(st.loc[old, 'edge_flag_disagree'].sum())
        rep.add(None, 'INFO', f'{old.sum()} files were converted with SHAPE = (400, 48, 192): their stored '
                              f'original_atEdge is computed on a transposed grid ({n_changed} events change '
                              'edge status). Use the new atEdge column; set SHAPE = (400, 192, 48) in '
                              'local_convert.py for future conversions.')
 
    out_parquet = args.output_dir / 'labels_with_sizes.parquet'
    combined.to_parquet(out_parquet, index=False)
 
    with pd.option_context('display.width', 250, 'display.max_columns', 30,
                           'display.float_format', '{:.3g}'.format):
        table = st.to_string()
    summary = combined[['x-size', 'y-size']].describe().to_string()
    contained = combined[~combined['atEdge']]
    verdict = 'FAIL' if rep.n('FAIL') else ('PASS with warnings' if rep.n('WARN') else 'PASS')
    lines = [
        f'Input:      {args.input_dir}',
        f'Indices:    d{indices[0]}..d{indices[-1]} ({len(indices)} expected, {len(frames)} usable)',
        f'Events:     {len(combined)} total, {len(contained)} contained (no border pixel > {EDGE_THRESHOLD:g} e-)',
        f'Size rule:  rows/cols whose summed charge is nonzero'
        + (f' after a {args.threshold:g} e- pixel threshold' if args.threshold > 0 else ''),
        f'Verdict:    {verdict}  ({rep.n("FAIL")} FAIL, {rep.n("WARN")} WARN, {rep.n("INFO")} INFO)',
        '', 'Issues:',
        *([f'  [{lvl}] ' + (f'd{i}: ' if i else '') + msg for i, lvl, msg in rep.issues] or ['  none']),
        '', 'Per-file statistics:',
        '  med_dx_pix / med_dy_pix: median |charge centroid - truth midplane| in pixels, (192,48) layout',
        '  *_T: same for the transposed (48,192) reading; should be much larger',
        '  med_ysize_minus_geom: median(y-size - 8|cotBeta|); med_xsize_minus_geom: median(x-size - 2|cotAlpha|)',
        '  q_1ns_over_final: recon3D total charge at 1 ns / at 4 ns (spot-checked events)',
        '  frac_no_eh: fraction of tracks with number_eh_pairs == 0 (should be 0)',
        '  med_eh_per_um: median e-h pairs per um of track path (silicon MIP: roughly 70-110)',
        table, '', 'Cluster sizes (all events):', summary, '',
        f'Wrote {out_parquet}',
    ]
    report_path = args.output_dir / 'validation_report.txt'
    report_path.write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:5]), flush=True)
    print(f'Report: {report_path}\nPlot input: {out_parquet}', flush=True)
    return 1 if rep.n('FAIL') else 0
 
 
if __name__ == '__main__':
    sys.exit(main())
