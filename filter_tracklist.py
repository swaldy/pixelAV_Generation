"""
filter_tracklist.py
Filter a pixelAV track list by |cotBeta| (col1) and save a subset.

This is how track_list_400k_cotBeta_lt2.txt was created from
track_list_2MflatPt_regression.txt (Danush's 2M flat-pT list).

Usage:
    python filter_tracklist.py <input_track_list> <max_cotbeta> <n_rows> <output_file>

Example:
    python filter_tracklist.py track_list_2MflatPt_regression.txt 2.0 400000 track_list_400k_cotBeta_lt2.txt
"""

import sys
import numpy as np
import random as random


def main():
    if len(sys.argv) != 5:
        print("Usage: python filter_tracklist.py <input> <max_cotbeta> <n_rows> <output>")
        sys.exit(1)

    input_file  = sys.argv[1]
    max_cotbeta = float(sys.argv[2])
    n_rows      = int(sys.argv[3])
    output_file = sys.argv[4]

    print(f"Loading {input_file} ...")
    data = np.loadtxt(input_file)
    print(f"Total rows: {len(data):,}")

    # col1 -> labels cotBeta (see column mapping in README)
    mask = np.abs(data[:, 0]) < max_cotbeta
    print(data[mask])
    rng=np.random.default_rng()
    #filtered=rng.choice(data[mask], size=n_rows,replace=False)
    filtered=data[mask]
    print(f"Rows with |cotBeta| < {max_cotbeta}: {mask.sum():,}")
    print(f"Saving first {len(filtered):,} rows to {output_file} ...")

    np.savetxt(output_file, filtered, fmt='%.6f %.6f %.4f %d %.6f %.6f %.4f')
    print("Done.")


if __name__ == "__main__":
    main()
