#!/usr/bin/env bash
# Condor wrapper for the same streaming converter used by run_local.sh.
set -eo pipefail
if [[ $# -ne 3 || ! "$1" =~ ^[0-9]+$ || ! "$2" =~ ^[0-9]+$ || ! "$3" =~ ^run_[0-9_]+$ ]]; then
    echo 'Usage: datagen.sh INDEX MAX_EVENTS RUN_ID' >&2
    exit 2
fi
i=$1
max_events=$2
run_id=$3
# LCG setup can reference unset variables: enable nounset only afterwards.
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh || exit $?
set -u
trap 'status=$?; echo "ERROR: line $LINENO exited with status $status" >&2; exit "$status"' ERR
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export ARROW_NUM_THREADS=1
ulimit -c 0

endpoint=root://eosuser.cern.ch
input_dir=/eos/user/s/swaldych/smart_pix/pixelav_48x192_10ps
# Unique run directory separates list-column outputs from legacy Parquet files.
remote_dir="$input_dir/parquet_streaming/$run_id/unflipped"
input_file="pixel_clusters_d${i}.out.gz"
echo "Worker: $(hostname); directory: $PWD; index: $i; max-events: $max_events"
python3 -c 'import sys, numpy, pyarrow; print(sys.version); print("NumPy", numpy.__version__, "PyArrow", pyarrow.__version__)'
command -v xrdcp
command -v xrdfs
xrdfs "$endpoint" mkdir -p "$remote_dir"
xrdcp "$endpoint/$input_dir/$input_file" "$input_file"
# No gunzip: local_convert.py streams directly from the compressed input.
/usr/bin/time -v python3 -u -X faulthandler local_convert.py "$input_file" \
    --output-dir converted --max-events "$max_events"

for kind in labels recon2D recon3D; do
    file="${kind}_d${i}.parquet"
    test -s "converted/unflipped/$file"
    xrdcp "converted/unflipped/$file" "$endpoint/$remote_dir/$file"
done
# Consumers should require this marker: the three uploads are not atomic.
printf 'index=%s\nmax_events=%s\n' "$i" "$max_events" > "SUCCESS_d${i}.txt"
xrdcp "SUCCESS_d${i}.txt" "$endpoint/$remote_dir/SUCCESS_d${i}.txt"
echo "SUCCESS: $endpoint/$remote_dir (d$i)"
