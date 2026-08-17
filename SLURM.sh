#!/bin/bash
#SBATCH --job-name=pixelav_200ps_48x192
#SBATCH --account=pi-badea
#SBATCH --partition=caslake
#SBATCH --array=1-40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=36:00:00
#SBATCH --output=/home/cmspixelsim/Documents/sarah/pixelAV_Generation/binary_output_48x192_200ps/%A_%a.txt

DATADIR=/home/cmspixelsim/Documents/sarah/pixelAV_Generation

FILEIND=$(printf "%05d" $((16400 + SLURM_ARRAY_TASK_ID)))

cd ${DATADIR}

/home/cmspixelsim/Documents/sarah/pixelAV_Generation/ppixelav2_list_trkpy_n_2f_200ps_48x192 \
    ${SLURM_ARRAY_TASK_ID} 5000

gzip pixel_clusters_d${FILEIND}.out
