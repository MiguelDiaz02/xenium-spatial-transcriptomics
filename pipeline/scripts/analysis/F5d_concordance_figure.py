#!/usr/bin/env python3
"""
F5d — Generate the Slingshot×DPT concordance figure (replaces aborted R block).

Run: conda run -n xenium_pipeline python F5d_concordance_figure.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import results_root

BASE = results_root()
INFILE = BASE / "03_phase3_spatial/F5b_macrophage_pseudotime/macrophage_slingshot_pseudotime.csv"
FIG = BASE / "figures/phase5_pseudotime_v2"
FIG.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE, index_col=0)
sl = df["slingPseudotime_1"]
dpt = df["dpt_pseudotime"]
mask = sl.notna() & dpt.notna()
sp_r, sp_p = spearmanr(sl[mask], dpt[mask])
pe_r, pe_p = pearsonr(sl[mask], dpt[mask])

fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor="white")

# Panel A: scatter with cell type
ax = axes[0]
colors = {"Macrophage_M1": "#1f77b4", "Macrophage_M2": "#d62728"}
for ct, c in colors.items():
    sub = df[df["cell_type_L2"] == ct]
    ax.scatter(sub["slingPseudotime_1"], sub["dpt_pseudotime"], s=2, alpha=0.25,
               c=c, label=f"{ct} (n={len(sub):,})")
ax.set_xlabel("Slingshot pseudotime (Lineage 1)", fontsize=11)
ax.set_ylabel("Diffusion pseudotime (DPT)", fontsize=11)
ax.set_title(f"Pseudotime concordance — Spearman r = {sp_r:.3f}",
             fontsize=12, fontweight="bold")
ax.legend(loc="upper left", fontsize=9, markerscale=4)

# Panel B: density heatmap
ax2 = axes[1]
ax2.hexbin(sl[mask], dpt[mask], gridsize=60, cmap="viridis", mincnt=1)
ax2.set_xlabel("Slingshot pseudotime (Lineage 1)", fontsize=11)
ax2.set_ylabel("Diffusion pseudotime (DPT)", fontsize=11)
ax2.set_title(f"Density (hexbin) — Spearman r = {sp_r:.3f}, Pearson r = {pe_r:.3f}",
              fontsize=12, fontweight="bold")

plt.tight_layout()
out = FIG / "F5c_concordance_slingshot_vs_dpt.png"
plt.savefig(out, dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")

# Also: pseudotime distribution per cell type
fig2, ax3 = plt.subplots(1, 2, figsize=(13, 4), facecolor="white")
for i, (col, name) in enumerate([("slingPseudotime_1", "Slingshot Lineage 1"),
                                  ("dpt_pseudotime", "Diffusion pseudotime")]):
    for ct, c in colors.items():
        sub = df[df["cell_type_L2"] == ct]
        ax3[i].hist(sub[col].dropna(), bins=50, alpha=0.55, color=c, label=ct)
    ax3[i].set_xlabel(name)
    ax3[i].set_ylabel("Cell count")
    ax3[i].set_title(name)
    ax3[i].legend()
plt.tight_layout()
out2 = FIG / "F5b_pseudotime_distributions.png"
plt.savefig(out2, dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved: {out2}")

print(f"\nConcordance summary:")
print(f"  Spearman r = {sp_r:.4f}, p = {sp_p:.2e}")
print(f"  Pearson  r = {pe_r:.4f}, p = {pe_p:.2e}")
print(f"  N cells (both pseudotimes non-NA): {mask.sum():,}")
