#!/usr/bin/env python3
"""Submit one gzip conversion job; --dry-run only writes the submit file."""
import argparse
from datetime import datetime
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=int, default=16401)
    parser.add_argument('--max-events', type=int, default=5,
                        help='Input events to process (default: 5; 0: entire file)')
    parser.add_argument('--memory-gb', type=int, default=8)
    parser.add_argument('--disk-gb', type=int, default=100,
                        help='Scratch space including input and all outputs')
    parser.add_argument('--flavour', choices=['longlunch', 'workday', 'tomorrow'], default='workday')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.index < 1 or args.max_events < 0 or args.memory_gb < 1 or args.disk_gb < 1:
        parser.error('Index and resource requests must be positive; max-events must be nonnegative')
    repo = Path(__file__).resolve().parent
    for name in ('datagen.sh', 'local_convert.py'):
        if not (repo / name).is_file():
            parser.error(f'Missing {repo / name}')
    run_id = datetime.now().strftime('run_%Y%m%d_%H%M%S_%f')
    logdir = repo / 'logs' / run_id
    logdir.mkdir(parents=True)
    relative = logdir.relative_to(repo)
    submit_file = logdir / f'd{args.index}.sub'
    submit_file.write_text(f'''universe                = vanilla
executable              = datagen.sh
arguments               = {args.index} {args.max_events} {run_id}
request_cpus            = 1
request_memory          = {args.memory_gb} GB
request_disk            = {args.disk_gb} GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = local_convert.py
transfer_output_files   = ""
output                  = {relative}/d{args.index}.out
error                   = {relative}/d{args.index}.err
log                     = {relative}/d{args.index}.log
+JobFlavour             = "{args.flavour}"
queue 1
''')
    print(f'Submit file: {submit_file}', flush=True)
    print(f'EOS outputs: /eos/user/s/swaldych/smart_pix/pixelav_48x192_10ps/parquet_streaming/{run_id}/unflipped', flush=True)
    if args.dry_run:
        print('Dry run: no job submitted.', flush=True)
        print(submit_file.read_text())
    else:
        subprocess.run(['condor_submit', str(submit_file.relative_to(repo))], cwd=repo, check=True)


if __name__ == '__main__':
    main()
