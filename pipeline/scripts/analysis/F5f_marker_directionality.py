#!/usr/bin/env python3
"""
F5f — Pseudotime directionality diagnostic via canonical M1/M2 markers.

Question: does the M1↔M2 trajectory show monotonic transition in canonical
markers (e.g., M1 markers decrease and M2 markers increase along Slingshot
pseudotime), or are the dynamics chaotic / bidirectional?

Markers used (panel-available):
  M1 (pro-inflammatory):  CD86, CXCL10, CXCL9, CD80
  M2 (immunosuppressive): CD163, MS4A4A, TREM2
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import results_root

BASE = results_root()
INDIR = BASE / "03_phase3_spatial/F5b_macrophage_pseudotime"
FIG = BASE / "figures/phase5_pseudotime_v2"

print("Loading ...")
counts = pd.read_csv(INDIR / "macrophage_counts_for_R.csv", index_col=0)
sl = pd.read_csv(INDIR / "macrophage_slingshot_pseudotime.csv", index_col=0)
assert len(counts) == len(sl)

M1 = ["CD86", "CXCL10", "CXCL9", "CD80"]
M2 = ["CD163", "MS4A4A", "TREM2"]
both = M1 + M2 + ["ADAM17"]

# Build long-form table
df = pd.DataFrame({"slingshot": sl["slingPseudotime_1"].values,
                   "dpt": sl["dpt_pseudotime"].values,
                   "celltype": sl["cell_type_L2"].values})
for g in both:
    if g in counts.columns:
        df[g] = counts[g].values

# Compute Spearman of marker vs Slingshot pseudotime
print("\nMarker correlation vs SLINGSHOT pseudotime (M1 forced as root):")
print(f"{'Marker':<10} {'Class':<3} {'Spearman r':>11} {'p':>9} {'Median M1':>11} {'Median M2':>11}")
results = []
for g in both:
    if g not in counts.columns:
        continue
    cls = "M1" if g in M1 else ("M2" if g in M2 else "?")
    r, p = spearmanr(df["slingshot"], df[g])
    mM1 = df.loc[df.celltype == "Macrophage_M1", g].median()
    mM2 = df.loc[df.celltype == "Macrophage_M2", g].median()
    print(f"{g:<10} {cls:<3} {r:>11.3f} {p:>9.1e} {mM1:>11.2f} {mM2:>11.2f}")
    results.append({"gene": g, "class": cls, "r_slingshot": r,
                     "p_slingshot": p, "median_M1": mM1, "median_M2": mM2})

# Visualization: smoothed mean expression per pseudotime bin
n_bins = 25
df["bin"] = pd.cut(df["slingshot"], bins=n_bins, labels=False)
binned = df.groupby("bin")[both].mean().reset_index()
binned_dpt = df.groupby("bin")["dpt"].mean().reset_index()

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor="white")

# Panel A — M1 markers
ax = axes[0]
for g in M1:
    if g in binned.columns:
        ax.plot(binned["bin"], binned[g], "-o", label=g, alpha=0.85, ms=4)
for g in M2:
    if g in binned.columns:
        ax.plot(binned["bin"], binned[g], "--s", label=g, alpha=0.85, ms=4)
if "ADAM17" in binned.columns:
    ax.plot(binned["bin"], binned["ADAM17"], "-^", label="ADAM17", alpha=0.85,
            ms=5, color="black", lw=1.6)
ax.set_xlabel("Slingshot pseudotime bin (M1 root → outward)", fontsize=11)
ax.set_ylabel("Mean expression (log-norm)", fontsize=11)
ax.set_title("A  Canonical M1/M2 markers along Slingshot trajectory\n"
             "Solid = M1 markers · Dashed = M2 markers · Black = ADAM17",
             loc="left", fontsize=11, fontweight="bold")
ax.legend(loc="best", fontsize=8, ncol=2)

# Panel B — markers vs DPT (independent direction)
ax2 = axes[1]
df["dpt_bin"] = pd.cut(df["dpt"], bins=n_bins, labels=False)
binned_dpt = df.groupby("dpt_bin")[both].mean().reset_index()
for g in M1:
    if g in binned_dpt.columns:
        ax2.plot(binned_dpt["dpt_bin"], binned_dpt[g], "-o", label=g, alpha=0.85, ms=4)
for g in M2:
    if g in binned_dpt.columns:
        ax2.plot(binned_dpt["dpt_bin"], binned_dpt[g], "--s", label=g, alpha=0.85, ms=4)
if "ADAM17" in binned_dpt.columns:
    ax2.plot(binned_dpt["dpt_bin"], binned_dpt["ADAM17"], "-^",
             label="ADAM17", alpha=0.85, ms=5, color="black", lw=1.6)
ax2.set_xlabel("DPT pseudotime bin (M1 cluster centroid root)", fontsize=11)
ax2.set_ylabel("Mean expression (log-norm)", fontsize=11)
ax2.set_title("B  Same markers along independent DPT pseudotime\n"
              "(used to test directionality robustness)",
              loc="left", fontsize=11, fontweight="bold")
ax2.legend(loc="best", fontsize=8, ncol=2)

plt.tight_layout()
out = FIG / "F5f_marker_directionality.png"
plt.savefig(out, dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved: {out}")

# Quantitative directionality summary
res_df = pd.DataFrame(results)
res_df.to_csv(INDIR / "F5f_marker_directionality.csv", index=False)
print("\n" + "=" * 60)
print("DIRECTIONALITY VERDICT")
print("=" * 60)
m1_rs = res_df.loc[res_df["class"] == "M1", "r_slingshot"]
m2_rs = res_df.loc[res_df["class"] == "M2", "r_slingshot"]
print(f"M1 markers — mean Spearman vs Slingshot pseudotime: {m1_rs.mean():.3f}")
print(f"M2 markers — mean Spearman vs Slingshot pseudotime: {m2_rs.mean():.3f}")
print(f"  Expected if M1→M2 directionality is real: M1 markers r<0, M2 markers r>0")
print(f"  Magnitude > 0.1 needed for clean directionality")
