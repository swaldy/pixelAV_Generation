"""Bounded-memory local PixelAV conversion; charge arrays use list columns."""
import argparse
from contextlib import ExitStack
import gzip
from itertools import islice
from pathlib import Path
import re
import tempfile

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

SHAPE = (400, 48, 192)  # Preserve the existing converter's reshape convention.
PIXELS = SHAPE[1] * SHAPE[2]
TRUTH = ['x-entry', 'y-entry', 'z-entry', 'n_x', 'n_y', 'n_z',
         'number_eh_pairs', 'y-local', 'pt']
DERIVED = ['cotAlpha', 'cotBeta', 'y-midplane', 'x-midplane']
LABEL_SCHEMA = pa.schema(
    [(name, pa.float64()) for name in TRUTH + DERIVED]
    + [('original_atEdge', pa.bool_()), ('chargeOriginal_atEdge', pa.float64()),
       ('nPixOriginalAbove1e_atEdge', pa.int64())])


def iter_events(path):
    """Yield one validated event; never retain previous events or full input text."""
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as source:
        header, stats = source.readline(), source.readline()
        if not header or not stats:
            raise ValueError('Missing input header or pixel metadata')
        print('Header:', header.strip(), flush=True)
        print('Pixelstats:', stats.strip(), flush=True)
        event = None
        truth = None
        need_truth = False
        time_index = -1
        offset = 0
        event_number = 0
        for line_number, line in enumerate(source, 3):
            text = line.strip()
            if not text:
                continue
            if text == '<cluster>':
                if event is not None or need_truth:
                    raise ValueError(f'Incomplete event before line {line_number}')
                need_truth = True
                time_index, offset = -1, 0
                event_number += 1
                continue
            if need_truth:
                truth = np.array([float(v) for v in text.split()], dtype=np.float64)
                if truth.size != len(TRUTH) or not np.isfinite(truth).all():
                    raise ValueError(f'Invalid truth row at line {line_number}')
                event = np.empty((SHAPE[0], PIXELS), dtype=np.float64)
                need_truth = False
                continue
            if 'time slice' in text:
                if event is None or (time_index >= 0 and offset != PIXELS):
                    raise ValueError(f'Unexpected/incomplete time slice at line {line_number}')
                time_index += 1
                if time_index >= SHAPE[0]:
                    raise ValueError(f'Too many time slices at line {line_number}')
                offset = 0
                continue
            if event is None or time_index < 0:
                raise ValueError(f'Unexpected data at line {line_number}')
            values = np.array([float(v) for v in text.split()], dtype=np.float64)
            end = offset + values.size
            if end > PIXELS or not np.isfinite(values).all():
                raise ValueError(f'Invalid pixel row at line {line_number}')
            event[time_index, offset:end] = values
            offset = end
            if time_index == SHAPE[0] - 1 and offset == PIXELS:
                yield event_number, truth, event
                event = None
                truth = None
        if event is not None or need_truth:
            raise ValueError('Input ended inside an incomplete event')


def label_row(truth, event):
    row = dict(zip(TRUTH, map(float, truth)))
    if row['n_z'] == 0:
        raise ValueError('n_z is zero; cannot calculate track slopes')
    row['cotAlpha'] = row['n_x'] / row['n_z']
    row['cotBeta'] = row['n_y'] / row['n_z']
    row['y-midplane'] = row['y-entry'] + row['cotBeta'] * (50 - row['z-entry'])
    row['x-midplane'] = row['x-entry'] + row['cotAlpha'] * (50 - row['z-entry'])
    matrix = event[-1].reshape(SHAPE[1:])
    # Preserve legacy boundary accounting, including duplicated corner entries.
    edges = np.concatenate([matrix[0, :], matrix[-1, :], matrix[:, 0], matrix[:, -1]])
    above = np.abs(edges) > 1
    row['original_atEdge'] = bool(above.any())
    row['chargeOriginal_atEdge'] = float(edges.sum()) if above.any() else 0.0
    row['nPixOriginalAbove1e_atEdge'] = int(above.sum()) if above.any() else 0
    return row


def charge_schema(shape):
    return pa.schema([('charge', pa.list_(pa.float64()))], metadata={
        b'layout': b'pixelav-local-list-v1',
        b'shape': ','.join(map(str, shape)).encode(),
        b'order': b'C; same flattened order as legacy converter'})


def charge_table(values, schema):
    flat = values.reshape(-1)
    column = pa.ListArray.from_arrays(pa.array([0, flat.size], type=pa.int32()),
                                    pa.array(flat, type=pa.float64()))
    return pa.Table.from_arrays([column], schema=schema)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Local pixel_clusters_dNNNNN.out.gz or .out')
    parser.add_argument('--output-dir', type=Path, default=Path('local_parquet_test'))
    parser.add_argument('--max-events', type=int, default=5,
                        help='Input events to examine (default: 5; 0: entire file)')
    args = parser.parse_args(argv)
    if args.max_events < 0:
        parser.error('--max-events must be nonnegative')
    match = re.fullmatch(r'pixel_clusters_d(\d+)\.out(?:\.gz)?', args.input.name)
    if not match:
        parser.error('Input filename must be pixel_clusters_dNNNNN.out[.gz]')
    if not args.input.is_file():
        parser.error(f'Input file does not exist: {args.input}')
    target = args.output_dir / 'unflipped'
    target.mkdir(parents=True, exist_ok=True)
    names = [f'{kind}_d{match[1]}.parquet' for kind in ('labels', 'recon2D', 'recon3D')]
    for name in names:
        if (target / name).exists():
            parser.error(f'Output exists: {target / name}; choose a new --output-dir')
    schemas = [LABEL_SCHEMA, charge_schema(SHAPE[1:]), charge_schema(SHAPE)]
    processed = written = 0
    # Only promote complete outputs after all requested events have been converted.
    with tempfile.TemporaryDirectory(prefix='.conversion-', dir=target) as temporary:
        paths = [Path(temporary) / name for name in names]
        with ExitStack() as stack:
            writers = [stack.enter_context(pq.ParquetWriter(
                path, schema, compression='zstd', use_dictionary=False,
                write_batch_size=65536)) for path, schema in zip(paths, schemas)]
            events = iter_events(args.input)
            stack.callback(events.close)
            selected = islice(events, args.max_events) if args.max_events else events
            for processed, truth, event in selected:
                if truth[2] == 100:  # Match the legacy unflipped selection.
                    writers[0].write_table(pa.Table.from_pylist([label_row(truth, event)], schema=schemas[0]))
                    writers[1].write_table(charge_table(event[-1], schemas[1]))
                    writers[2].write_table(charge_table(event, schemas[2]))
                    written += 1
                print(f'Processed {processed}; saved {written} unflipped events', flush=True)
        if processed == 0:
            raise ValueError('No complete events found')
        for path, name in zip(paths, names):
            path.rename(target / name)
    print(f'Done: {processed} input events, {written} output rows per file.', flush=True)
    print(f'Output directory: {target.resolve()}', flush=True)
    if args.max_events:
        print('Event-limited test: the remainder of the input was not validated.', flush=True)


if __name__ == '__main__':
    main()
