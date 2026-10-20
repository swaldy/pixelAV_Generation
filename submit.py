#!/usr/bin/env python3
"""Submit one gzip conversion job; --dry-run only writes the submit file."""
import argparse
from datetime import datetime
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=int, action='append', default=None,
                        help='One file index; repeat for several (default: 16401)')
    parser.add_argument('--range', nargs=2, type=int, metavar=('FIRST', 'LAST'),
                        help='Inclusive index range, e.g. --range 16401 16440')
    parser.add_argument('--max-events', type=int, default=5,
                        help='Input events to process (default: 5; 0: entire file)')
    parser.add_argument('--memory-gb', type=int, default=8)
    parser.add_argument('--disk-gb', type=int, default=100,
                        help='Scratch space including input and all outputs')
    parser.add_argument('--flavour',
                        choices=['espresso', 'microcentury', 'longlunch', 'workday',
                                 'tomorrow', 'testmatch', 'nextweek'],
                        default='workday')
    parser.add_argument('--run-id', default=None,
                        help='Reuse an existing run directory instead of making a new one')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.index and args.range:
        parser.error('Use --index or --range, not both')
    if args.range:
        first, last = args.range
        if first < 1 or last < first:
            parser.error('--range must be two ascending positive integers')
        indices = list(range(first, last + 1))
    else:
        indices = sorted(set(args.index or [16401]))
    if any(i < 1 for i in indices):
        parser.error('Indices must be positive')
    if args.max_events < 0 or args.memory_gb < 1 or args.disk_gb < 1:
        parser.error('Resource requests must be positive; max-events must be nonnegative')
    repo = Path(__file__).resolve().parent
    for name in ('datagen.sh', 'local_convert.py'):
        if not (repo / name).is_file():
            parser.error(f'Missing {repo / name}')
    run_id = args.run_id or datetime.now().strftime('run_%Y%m%d_%H%M%S_%f')
    if not run_id.startswith('run_'):
        parser.error('--run-id must start with "run_" (datagen.sh validates it)')
    logdir = repo / 'logs' / run_id
    logdir.mkdir(parents=True, exist_ok=True)
    relative = logdir.relative_to(repo)
    tag = f'd{indices[0]}' if len(indices) == 1 else f'd{indices[0]}_d{indices[-1]}'
    submit_file = logdir / f'{tag}.sub'
    queue = ' '.join(str(i) for i in indices)
    logdir = repo / 'logs' / run_id
    logdir.mkdir(parents=True)
    relative = logdir.relative_to(repo)
    submit_file = logdir / f'd{args.index}.sub'
    submit_file.write_text(f'''universe                = vanilla
executable              = datagen.sh
arguments               = $(index) {args.max_events} {run_id}
request_cpus            = 1
request_memory          = {args.memory_gb} GB
request_disk            = {args.disk_gb} GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = local_convert.py
transfer_output_files   = ""
output                  = {relative}/d$(index).out
error                   = {relative}/d$(index).err
log                     = {relative}/d$(index).log
+JobFlavour             = "{args.flavour}"
queue index in ({queue})
''')
    print(f'Submit file: {submit_file}', flush=True)
    print(f'Indices    : {len(indices)} job(s): {indices[0]}..{indices[-1]}', flush=True)
    print('EOS outputs: /eos/user/s/swaldych/smart_pix/pixelav_48x192_10ps/'
          f'parquet_streaming/{run_id}/unflipped', flush=True)
    if args.dry_run:
        print('Dry run: no job submitted.', flush=True)
        print(submit_file.read_text())
    else:
        subprocess.run(['condor_submit', str(submit_file.relative_to(repo))], cwd=repo, check=True)


if __name__ == '__main__':
    main()
