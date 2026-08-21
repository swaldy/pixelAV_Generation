"""
datagen_10ps_16x16.py
Parquet conversion for 10ps 16x16 pixelAV datasets.

Changes vs upstream datagen.py (github.com/smart-pix/filter/semiprocessing_datagen/datagen.py):
  - len(cur_slice) == 16*16  (was 13*21)
  - len(cur_cluster) == 400  (was 20)
  - reads .out.gz directly (no manual gunzip needed)

Usage (run from dataset directory):
    python datagen_10ps_16x16.py <file_index>

Example:
    cd /project/badea/smartpix/harshul/dataset_3sr_16x16_50x12P5_centeredIncidence_10ps_243k
    python /home/harshul/smartpixels/datagen_10ps_16x16.py 16401
"""

import sys
import gzip
import os
import numpy as np
import pandas as pd


def split(index, df1, df2, df3):
    df1.columns = df1.columns.astype(str)
    df2.columns = df2.columns.astype(str)
    df3.columns = df3.columns.astype(str)

    os.makedirs("unflipped", exist_ok=True)

    mask = df1['z-entry'] == 100
    tag  = "d" + str(index)

    # NOTE: labels parquet does NOT include contained/peripheral cluster columns
    # (e.g. original_atEdge, chargeOriginal_atEdge per Danush's PDF definition).
    # Containment is currently computed only in validation scripts (plot_truth_individual.py),
    # not persisted to the output parquets.
    df1[mask].to_parquet(f"unflipped/labels_{tag}.parquet")
    df2[mask].to_parquet(f"unflipped/recon2D_{tag}.parquet")
    df3[mask].to_parquet(f"unflipped/recon3D_{tag}.parquet")


def parseFile(filein):
    gz_path  = filein + ".gz"
    out_path = filein

    if os.path.exists(gz_path):
        f = gzip.open(gz_path, 'rt')
    elif os.path.exists(out_path):
        f = open(out_path, 'r')
    else:
        raise FileNotFoundError(f"Neither {gz_path} nor {out_path} found")

    lines = f.readlines()
    f.close()

    header     = lines[0].strip()
    pixelstats = lines[1].strip()
    print("Header:", header)
    print("Pixelstats:", pixelstats)

    readyToGetTruth     = False
    readyToGetTimeSlice = False

    clusterctr    = 0
    cluster_truth = []
    timeslice     = 0
    cur_slice     = []
    cur_cluster   = []
    events        = []

    for line in lines:
        if "<cluster>" in line:
            readyToGetTruth     = True
            readyToGetTimeSlice = False
            clusterctr += 1
            cur_cluster = []
            timeslice   = 0
            continue

        if readyToGetTruth:
            cluster_truth.append(line.strip().split())
            readyToGetTruth = False
            continue

        if "time slice" in line:
            readyToGetTimeSlice = True
            cur_slice = []
            timeslice += 1
            continue

        if readyToGetTimeSlice:
            cur_row    = line.strip().split()
            cur_slice += [float(item) for item in cur_row]

            if len(cur_slice) == 16 * 16:
                cur_cluster.append(cur_slice)

            if len(cur_cluster) == 400:
                events.append(cur_cluster)
                readyToGetTimeSlice = False

    print("Number of clusters =", len(cluster_truth))
    print("Number of events =",   len(events))
    print("Number of time slices in cluster =", len(events[0]))

    arr_truth  = np.array(cluster_truth)
    arr_events = np.array(events)

    return arr_events, arr_truth


def main():
    if len(sys.argv) != 2:
        print("Usage: python datagen_10ps_16x16.py <file_index>")
        sys.exit(1)

    index  = int(sys.argv[1])
    tag    = "d" + str(index)
    arr_events, arr_truth = parseFile(filein=f"./pixel_clusters_d{index}.out")

    df = pd.DataFrame(arr_truth,
                      columns=['x-entry', 'y-entry', 'z-entry',
                               'n_x', 'n_y', 'n_z',
                               'number_eh_pairs', 'y-local', 'pt'])
    for col in df.columns:
        df[col] = df[col].astype(float)

    sensor_thickness = 100  # um
    df['cotAlpha']   = df['n_x'] / df['n_z']
    df['cotBeta']    = df['n_y'] / df['n_z']
    df['y-midplane'] = df['y-entry'] + df['cotBeta']  * (sensor_thickness / 2 - df['z-entry'])
    df['x-midplane'] = df['x-entry'] + df['cotAlpha'] * (sensor_thickness / 2 - df['z-entry'])

    df2 = pd.DataFrame([np.array(e[-1]).flatten() for e in arr_events])
    df3 = pd.DataFrame([np.array(e).flatten()     for e in arr_events])

    split(index, df, df2, df3)


if __name__ == "__main__":
    main()
