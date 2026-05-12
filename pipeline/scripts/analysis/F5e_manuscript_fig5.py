#!/usr/bin/env python3
"""
F5e — Generate the manuscript Figure 5 (M1↔M2 macrophage continuum).
Replaces the OLD `Fig5B_pseudotime_proportions.png` with the audited M1↔M2 narrative.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
INFILE = BASE / "human_lung_cancer/results/03_phase3_spatial/F5b_macrophage_pseudotime/macrophage_slingshot_pseudotime.csv"
OUT_MS = BASE / "manuscript/figures/Fig5B_pseudotime_proportions.png"
OUT_RES = BASE / "human_lung_cancer/results/figures/ms_main_fig5B_pseudotime_proportions.png"

df = pd.read_csv(INFILE, index_col=0)
sl = df["slingPseudotime_1"]
dpt = df["dpt_pseudotime"]
mask = sl.notna() & dpt.notna()
sp_r, _ = spearmanr(sl[mask], dpt[mask])

colors = {"Macrophage_M1": "#1f77b4", "Macrophage_M2": "#d62728"}

fig = plt.figure(figsize=(13, 5.5), facecolor="white")
gs = GridSpec(1, 3, figure=fig, wspace=0.32)

# Panel A: UMAP colored by cell type + DPT
ax = fig.add_subplot(gs[0, 0])
sc = ax.scatter(df["umap_1"], df["umap_2"], c=df["dpt_pseudotime"],
                s=2, cmap="viridis", alpha=0.7)
plt.colorbar(sc, ax=ax, label="DPT pseudotime", shrink=0.7)
# Overlay M1/M2 centroids
for ct, c in colors.items():
    sub = df[df["cell_type_L2"] == ct]
    cx, cy = sub["umap_1"].mean(), sub["umap_2"].mean()
    ax.text(cx, cy, ct.replace("Macrophage_", ""),
            fontsize=14, fontweight="bold", color=c,
            ha="center", va="center",
            bbox=dict(boxstyle="round", fc="white", ec=c, lw=2, alpha=0.85))
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("A  Macrophage manifold (PAGA-init UMAP)\ncolored by diffusion pseudotime",
             loc="left", fontsize=11, fontweight="bold")
ax.set_xticks([])
ax.set_yticks([])

# Panel B: Slingshot vs DPT concordance (hexbin)
ax2 = fig.add_subplot(gs[0, 1])
hb = ax2.hexbin(sl[mask], dpt[mask], gridsize=50, cmap="viridis", mincnt=1)
plt.colorbar(hb, ax=ax2, label="Cells / hex", shrink=0.7)
ax2.set_xlabel("Slingshot pseudotime (Lineage 1)")
ax2.set_ylabel("Diffusion pseudotime (DPT)")
ax2.set_title(f"B  Two-method concordance\nSpearman r = {sp_r:.3f}, n = {mask.sum():,}",
              loc="left", fontsize=11, fontweight="bold")

# Panel C: distribution of Slingshot pseudotime by cell_type
ax3 = fig.add_subplot(gs[0, 2])
for ct, c in colors.items():
    sub = df[df["cell_type_L2"] == ct]
    ax3.hist(sub["slingPseudotime_1"].dropna(), bins=50, alpha=0.55,
             color=c, label=f"{ct.replace('Macrophage_', '')}  (n={len(sub):,})")
ax3.set_xlabel("Slingshot pseudotime (Lineage 1)")
ax3.set_ylabel("Cell count")
ax3.set_title("C  M1-rooted Slingshot ordering\nseparates M1 and M2 along the trajectory",
              loc="left", fontsize=11, fontweight="bold")
ax3.legend(loc="upper right", fontsize=10)

plt.tight_layout()
for outpath in (OUT_MS, OUT_RES):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"Saved: {outpath}")
plt.close()
