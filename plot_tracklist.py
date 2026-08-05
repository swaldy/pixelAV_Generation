import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(
    description="Plot PixelAV track-list distributions."
)
parser.add_argument("tracklist", help="PixelAV track-list text file")
parser.add_argument(
    "--output",
    default="tracklist_distributions.png",
    help="Output PNG filename",
)
args = parser.parse_args()


# Expected columns:
# 0: cotBeta according to final label convention
# 1: cotAlpha according to final label convention
# 2: pion momentum
# 3: flipped
# 4: x-local
# 5: y-local
# 6: signed pT
data = np.loadtxt(args.tracklist)

if data.ndim == 1:
    data = data.reshape(1, -1)

if data.shape[1] != 7:
    raise ValueError(
        f"Expected 7 columns, but found {data.shape[1]}"
    )

cot_beta = data[:, 0]
cot_alpha = data[:, 1]
pt = data[:, 6]
abs_pt = np.abs(pt)

print(f"File: {args.tracklist}")
print(f"Number of tracks: {len(data):,}")

print(
    "cotBeta:",
    f"min={cot_beta.min():.5f},",
    f"max={cot_beta.max():.5f},",
    f"mean={cot_beta.mean():.5f},",
    f"std={cot_beta.std():.5f}",
)

print(
    "cotAlpha:",
    f"min={cot_alpha.min():.5f},",
    f"max={cot_alpha.max():.5f},",
    f"mean={cot_alpha.mean():.5f},",
    f"std={cot_alpha.std():.5f}",
)

print(
    "signed pT:",
    f"min={pt.min():.5f},",
    f"max={pt.max():.5f},",
    f"mean={pt.mean():.5f},",
    f"std={pt.std():.5f}",
)

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

axes[0, 0].hist(
    cot_beta,
    bins=100,
    histtype="step",
    linewidth=1.5,
    color="tab:blue",
)
axes[0, 0].axvline(-2, color="black", linestyle="--", alpha=0.6)
axes[0, 0].axvline(2, color="black", linestyle="--", alpha=0.6)
axes[0, 0].set_xlabel(r"$\mathrm{cot}\beta$")
axes[0, 0].set_ylabel("Tracks")
axes[0, 0].set_title(r"$\mathrm{cot}\beta$ distribution")
axes[0, 0].grid(alpha=0.25)

axes[0, 1].hist(
    cot_alpha,
    bins=100,
    histtype="step",
    linewidth=1.5,
    color="tab:orange",
)
axes[0, 1].set_xlabel(r"$\mathrm{cot}\alpha$")
axes[0, 1].set_ylabel("Tracks")
axes[0, 1].set_title(r"$\mathrm{cot}\alpha$ distribution")
axes[0, 1].grid(alpha=0.25)

axes[1, 0].hist(
    pt,
    bins=100,
    histtype="step",
    linewidth=1.5,
    color="tab:green",
)
axes[1, 0].set_xlabel(r"Signed $p_T$ [GeV]")
axes[1, 0].set_ylabel("Tracks")
axes[1, 0].set_title(r"Signed $p_T$ distribution")
axes[1, 0].grid(alpha=0.25)

axes[1, 1].hist(
    abs_pt,
    bins=100,
    histtype="step",
    linewidth=1.5,
    color="tab:red",
)
axes[1, 1].set_xlabel(r"$|p_T|$ [GeV]")
axes[1, 1].set_ylabel("Tracks")
axes[1, 1].set_title(r"$|p_T|$ distribution")
axes[1, 1].grid(alpha=0.25)

fig.suptitle(
    f"PixelAV track list: {Path(args.tracklist).name}\n"
    f"{len(data):,} tracks",
    fontsize=13,
)

fig.tight_layout()
fig.savefig(args.output, dpi=200, bbox_inches="tight")

print(f"Saved plot: {args.output}")
