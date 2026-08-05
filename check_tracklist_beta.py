"""
check_tracklist_beta.py
Print cotBeta (col1) and cotAlpha (col2) statistics for a pixelAV track list.

Column mapping (matches datagen convention):
  col1 -> labels cotBeta  (cotalpha_track)
  col2 -> labels cotAlpha (cotbeta_track)

Usage:
    python check_tracklist_beta.py <track_list.txt>
"""

import sys
import numpy as np

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "track_list.txt"
    print(f"Loading {path} ...")

    # Read in chunks to handle large files without hitting memory issues
    col1, col2 = [], []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if len(vals) < 2:
                continue
            col1.append(float(vals[0]))
            col2.append(float(vals[1]))

    col1 = np.array(col1)
    col2 = np.array(col2)

    print(f"\nTotal rows: {len(col1):,}")

    for name, arr in [("cotBeta (col1)", col1), ("cotAlpha (col2)", col2)]:
        # Filter obvious outliers (|val| > 50 is unphysical for cot angles in CMS)
        clean = arr[np.abs(arr) <= 50]
        frac_out = 1 - len(clean) / len(arr)
        print(f"\n{name}:")
        print(f"  rows used (|val|<=50): {len(clean):,}  ({frac_out:.3%} outliers excluded)")
        print(f"  mean={clean.mean():.4f}  std={clean.std():.4f}")
        print(f"  min={clean.min():.4f}  max={clean.max():.4f}")
        pcts = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        vals = np.percentile(clean, pcts)
        print("  percentiles: " + "  ".join(f"p{p}={v:.3f}" for p, v in zip(pcts, vals)))

if __name__ == "__main__":
    main()
