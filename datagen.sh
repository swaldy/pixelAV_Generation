#!/bin/bash
# datagen.sh
# Fetches a raw pixel_clusters .out.gz file from EOS, runs datagen.py to add the
# contained-cluster columns, and pushes the updated labels/recon2D/recon3D parquet
# files back to EOS, replacing the existing (column-less) versions.

i=$1

source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

INPUT_EOS_DIR=/eos/project/s/smartpix-box/pixelAV_datasets/unshuffled_DO_NOT_DELETE/temporary/dataset_3srb_16x16_50x12P5_centeredIncidence_10ps_300k
outdir=/eos/project/s/smartpix-box/pixelAV_datasets/unshuffled_DO_NOT_DELETE/temporary/dataset_3srb_16x16_50x12P5_centeredIncidence_10ps_300k/dataset_3srb_16x16_50x12P5_centeredIncidence_10ps_300k_parquets/

mkdir -p unflipped

xrdcp root://eosproject.cern.ch/$INPUT_EOS_DIR/pixel_clusters_d${i}.out.gz pixel_clusters_d${i}.out.gz
pwd

gunzip pixel_clusters_d${i}.out.gz

python3 datagen.py $i

xrdcp -f unflipped/labels_d${i}.parquet root://eosproject.cern.ch/$outdir/unflipped/labels_d${i}.parquet
xrdcp -f unflipped/recon2D_d${i}.parquet root://eosproject.cern.ch/$outdir/unflipped/recon2D_d${i}.parquet
xrdcp -f unflipped/recon3D_d${i}.parquet root://eosproject.cern.ch/$outdir/unflipped/recon3D_d${i}.parquet
