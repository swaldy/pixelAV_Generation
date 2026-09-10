#!/usr/bin/env bash
# Uses a local Python environment; no Condor, EOS, CVMFS, or decompression step.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PIXELAV_PYTHON:-$script_dir/.venv/bin/python}"
if [[ ! -x "$python_bin" ]]; then
    echo "Missing Python environment. In $script_dir run:" >&2
    echo 'python3 -m venv .venv' >&2
    echo '.venv/bin/python -m pip install -r requirements-local.txt' >&2
    exit 1
fi
exec "$python_bin" -u "$script_dir/local_convert.py" "$@"
