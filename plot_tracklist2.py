import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(
    description="Plot PixelAV track-list distributions and cotBeta passing fractions."
)

parser.add_argument(
    "tracklist",
    help="PixelAV track-list text file"
)

parser.add_argument(
    "--output",
    default="tracklist_distributions.png",
    help="Output PNG filename",
)

parser.add_argument(
    "--cotbeta-max",
    type=float,
    default=2.0,
    help="Maximum |cotBeta| used for passing fraction (default: 2.0)",
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


# -------------------------------------------------------------------
# Extract variables
# -------------------------------------------------------------------

cot_beta = data[:, 0]
cot_alpha = data[:, 1]
pt = data[:, 6]
abs_pt = np.abs(pt)


# -------------------------------------------------------------------
# Define cotBeta passing requirement
# -------------------------------------------------------------------

cotbeta_limit = args.cotbeta_max

pass_mask = np.abs(cot_beta) < cotbeta_limit

n_pass = np.count_nonzero(pass_mask)
overall_fraction = n_pass / len(data)


# -------------------------------------------------------------------
# Print basic information
# -------------------------------------------------------------------

print(f"File: {args.tracklist}")
print(f"Number of tracks: {len(data):,}")

print(
    "cotBeta:",
    f"min={cot_beta.min():.5f},",
    f"max={cot_beta.max():.5f},",
    f"mean={cot_beta.mean():.5f},",
    f"std={cot_beta.std():.5f}"
)

print(
    "cotAlpha:",
    f"min={cot_alpha.min():.5f},",
    f"max={cot_alpha.max():.5f},",
    f"mean={cot_alpha.mean():.5f},",
    f"std={cot_alpha.std():.5f}"
)

print(
    "signed pT:",
    f"min={pt.min():.5f},",
    f"max={pt.max():.5f},",
    f"mean={pt.mean():.5f},",
    f"std={pt.std():.5f}"
)

print()
print(
    f"Tracks passing |cotBeta| < {cotbeta_limit}: "
    f"{n_pass:,} / {len(data):,}"
)

print(
    f"Overall passing fraction: {overall_fraction:.4f} "
    f"({100 * overall_fraction:.2f}%)"
)


# -------------------------------------------------------------------
# Calculate passing fraction versus signed pT
# -------------------------------------------------------------------

pt_bins = np.linspace(
    pt.min(),
    pt.max(),
    501
)

all_signed_counts, pt_edges = np.histogram(
    pt,
    bins=pt_bins
)

pass_signed_counts, _ = np.histogram(
    pt[pass_mask],
    bins=pt_bins
)

pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])

signed_passing_fraction = np.divide(
    pass_signed_counts,
    all_signed_counts,
    out=np.full(
        all_signed_counts.shape,
        np.nan,
        dtype=float
    ),
    where=all_signed_counts > 0
)


# -------------------------------------------------------------------
# Calculate passing fraction versus |pT|
# -------------------------------------------------------------------

abs_pt_bins = np.linspace(
    abs_pt.min(),
    abs_pt.max(),
    501
)

all_abs_counts, abs_pt_edges = np.histogram(
    abs_pt,
    bins=abs_pt_bins
)

pass_abs_counts, _ = np.histogram(
    abs_pt[pass_mask],
    bins=abs_pt_bins
)

abs_pt_centers = 0.5 * (
    abs_pt_edges[:-1] + abs_pt_edges[1:]
)

abs_passing_fraction = np.divide(
    pass_abs_counts,
    all_abs_counts,
    out=np.full(
        all_abs_counts.shape,
        np.nan,
        dtype=float
    ),
    where=all_abs_counts > 0
)


# -------------------------------------------------------------------
# Make plots
# -------------------------------------------------------------------

fig, axes = plt.subplots(
    3,
    2,
    figsize=(12, 13)
)


# -------------------------------------------------------------------
# Calculate histogram bins for distributions
# -------------------------------------------------------------------

n_bins = 500


# cotBeta
cot_beta_counts, cot_beta_edges = np.histogram(
    cot_beta,
    bins=n_bins
)

cot_beta_centers = 0.5 * (
    cot_beta_edges[:-1] + cot_beta_edges[1:]
)

cot_beta_widths = np.diff(cot_beta_edges)


# cotAlpha
cot_alpha_counts, cot_alpha_edges = np.histogram(
    cot_alpha,
    bins=n_bins
)

cot_alpha_centers = 0.5 * (
    cot_alpha_edges[:-1] + cot_alpha_edges[1:]
)

cot_alpha_widths = np.diff(cot_alpha_edges)


# Signed pT
pt_counts, pt_hist_edges = np.histogram(
    pt,
    bins=n_bins
)

pt_hist_centers = 0.5 * (
    pt_hist_edges[:-1] + pt_hist_edges[1:]
)

pt_hist_widths = np.diff(pt_hist_edges)


# |pT|
abs_pt_counts, abs_pt_hist_edges = np.histogram(
    abs_pt,
    bins=n_bins
)

abs_pt_hist_centers = 0.5 * (
    abs_pt_hist_edges[:-1] + abs_pt_hist_edges[1:]
)

abs_pt_hist_widths = np.diff(abs_pt_hist_edges)


# -------------------------
# cotBeta distribution
# -------------------------

axes[0, 0].bar(
    cot_beta_centers,
    cot_beta_counts,
    width=cot_beta_widths,
    align="center",
    color="tab:blue",
    linewidth=0,
)

axes[0, 0].axvline(
    -cotbeta_limit,
    color="black",
    linestyle="--",
    alpha=0.6
)

axes[0, 0].axvline(
    cotbeta_limit,
    color="black",
    linestyle="--",
    alpha=0.6
)

axes[0, 0].set_xlabel(r"$\mathrm{cot}\beta$")
axes[0, 0].set_ylabel("Tracks")

axes[0, 0].set_title(
    r"$\mathrm{cot}\beta$ distribution"
)

axes[0, 0].grid(alpha=0.25)


# -------------------------
# cotAlpha distribution
# -------------------------

axes[0, 1].bar(
    cot_alpha_centers,
    cot_alpha_counts,
    width=cot_alpha_widths,
    align="center",
    color="tab:orange",
    linewidth=0,
)

axes[0, 1].set_xlabel(r"$\mathrm{cot}\alpha$")
axes[0, 1].set_ylabel("Tracks")

axes[0, 1].set_title(
    r"$\mathrm{cot}\alpha$ distribution"
)

axes[0, 1].grid(alpha=0.25)


# -------------------------
# Signed pT distribution
# -------------------------

axes[1, 0].bar(
    pt_hist_centers,
    pt_counts,
    width=pt_hist_widths,
    align="center",
    color="tab:green",
    linewidth=0,
)

axes[1, 0].set_xlabel(r"Signed $p_T$ [GeV]")
axes[1, 0].set_ylabel("Tracks")

axes[1, 0].set_title(
    r"Signed $p_T$ distribution"
)

axes[1, 0].grid(alpha=0.25)


# -------------------------
# |pT| distribution
# -------------------------

axes[1, 1].bar(
    abs_pt_hist_centers,
    abs_pt_counts,
    width=abs_pt_hist_widths,
    align="center",
    color="tab:red",
    linewidth=0,
)

axes[1, 1].set_xlabel(r"$|p_T|$ [GeV]")
axes[1, 1].set_ylabel("Tracks")

axes[1, 1].set_title(
    r"$|p_T|$ distribution"
)

axes[1, 1].grid(alpha=0.25)


# -------------------------
# Passing fraction vs signed pT
# -------------------------

pt_bin_widths = np.diff(pt_edges)

axes[2, 0].bar(
    pt_centers,
    signed_passing_fraction,
    width=pt_bin_widths,
    align="center",
    color="tab:blue",
    linewidth=0,
)

axes[2, 0].set_xlabel(r"Signed $p_T$ [GeV]")
axes[2, 0].set_ylabel("Passing fraction")

axes[2, 0].set_title(
    rf"Fraction passing $|\mathrm{{cot}}\beta| < {cotbeta_limit}$"
)

axes[2, 0].set_ylim(0, 1.05)
axes[2, 0].grid(alpha=0.25)


# -------------------------
# Passing fraction vs |pT|
# -------------------------

abs_pt_bin_widths = np.diff(abs_pt_edges)

axes[2, 1].bar(
    abs_pt_centers,
    abs_passing_fraction,
    width=abs_pt_bin_widths,
    align="center",
    color="tab:blue",
    linewidth=0,
)

axes[2, 1].set_xlabel(r"$|p_T|$ [GeV]")
axes[2, 1].set_ylabel("Passing fraction")

axes[2, 1].set_title(
    rf"Fraction passing $|\mathrm{{cot}}\beta| < {cotbeta_limit}$"
)

axes[2, 1].set_ylim(0, 1.05)
axes[2, 1].grid(alpha=0.25)


# -------------------------------------------------------------------
# Overall title
# -------------------------------------------------------------------

fig.suptitle(
    f"PixelAV track list: {Path(args.tracklist).name}\n"
    f"{len(data):,} tracks",
    fontsize=13,
)

fig.tight_layout()

fig.savefig(
    args.output,
    dpi=200,
    bbox_inches="tight"
)

# -------------------------------------------------------------------
# Histogram sanity checks
# -------------------------------------------------------------------

counts, edges = np.histogram(
    cot_beta,
    bins=500
)
print("cotBeta entries:", counts.sum())

counts, edges = np.histogram(
    cot_alpha,
    bins=500
)
print("cotAlpha entries:", counts.sum())

counts, edges = np.histogram(
    pt,
    bins=500
)
print("Signed pT entries:", counts.sum())

counts, edges = np.histogram(
    abs_pt,
    bins=500
)
print("|pT| entries:", counts.sum())

print(f"Saved plot: {args.output}")
