#!/usr/bin/env python3
"""
F5g — PAGA macrofágico anotado biológicamente.

Mejora la visualización original de PAGA (que sólo mostraba clusters numerados
0..10 sin contexto). Esta versión:
  - Tamaño de nodo = # de células en el sub-cluster
  - Color de nodo = % de células M2 (azul=M1, rojo=M2, púrpura=mixto)
  - Grosor de arista = conectividad PAGA
  - Layout reproducible
  - Marcadores M1 (CD86, CD80, CXCL9, CXCL10) y M2 (CD163, MS4A4A, TREM2)
    overlay como dot-plot por sub-cluster
"""
from __future__ import annotations
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

import spatialdata as sd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

warnings.filterwarnings("ignore", category=FutureWarning)

BASE = project_root()
ZARR = BASE / "human_lung_cancer/results/sdata.zarr"
FIG = BASE / "human_lung_cancer/results/figures/phase5_pseudotime_audited"

print("Loading sdata.zarr ...")
sdata = sd.read_zarr(str(ZARR))
adata = list(sdata.tables.values())[0]
mask = adata.obs["cell_type_L2"].isin(["Macrophage_M1", "Macrophage_M2"])
adata_m = adata[mask].copy()
adata_m.obs["cell_type_L2"] = adata_m.obs["cell_type_L2"].cat.remove_unused_categories()
print(f"Macrophages: {adata_m.n_obs}")

# Reproduce F5b pipeline (deterministic via seed=42)
sc.tl.pca(adata_m, n_comps=30, random_state=42)
sc.pp.neighbors(adata_m, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.leiden(adata_m, resolution=0.5, random_state=42, key_added="macro_leiden",
             flavor="igraph", directed=False, n_iterations=2)
sc.tl.paga(adata_m, groups="macro_leiden")
# Compute PAGA positions (required before UMAP init_pos="paga")
sc.pl.paga(adata_m, plot=False)
sc.tl.umap(adata_m, init_pos="paga", random_state=42)

# Compute per-cluster stats
clusters = sorted(adata_m.obs["macro_leiden"].cat.categories,
                  key=lambda x: int(x))
cluster_stats = []
for c in clusters:
    sub = adata_m.obs[adata_m.obs["macro_leiden"] == c]
    n = len(sub)
    n_m1 = (sub["cell_type_L2"] == "Macrophage_M1").sum()
    n_m2 = (sub["cell_type_L2"] == "Macrophage_M2").sum()
    pct_m2 = 100 * n_m2 / n
    cluster_stats.append({"cluster": c, "n_cells": n,
                          "n_M1": n_m1, "n_M2": n_m2, "pct_M2": pct_m2})
stats_df = pd.DataFrame(cluster_stats)
print(stats_df.to_string(index=False))

# Color map: blue (M1) → light gray (mixed) → red (M2)
m1m2_cmap = LinearSegmentedColormap.from_list(
    "m1m2", ["#1f77b4", "#dddddd", "#d62728"], N=256
)

# Add a per-cell `pct_M2_in_cluster` column so sc.pl.paga can color by mean
adata_m.obs["is_M2_int"] = (adata_m.obs["cell_type_L2"] == "Macrophage_M2").astype(int)
cluster_to_pct = dict(zip(stats_df.cluster, stats_df.pct_M2 / 100.0))
adata_m.obs["pct_M2_cluster"] = (
    adata_m.obs["macro_leiden"].astype(str).map(cluster_to_pct).astype(float)
)

# ── Figure: 2 panels (annotated PAGA + marker dotplot) ───────────────────────
M1 = ["CD86", "CD80", "CXCL9", "CXCL10"]
M2 = ["CD163", "MS4A4A", "TREM2"]
markers_in = [g for g in (M1 + M2) if g in adata_m.var_names]

fig = plt.figure(figsize=(14, 6.5), facecolor="white")
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=0.3)

# ---- Panel A: Annotated PAGA ----
ax = fig.add_subplot(gs[0, 0])
node_colors = [m1m2_cmap(stats_df.loc[stats_df.cluster == c, "pct_M2"].iloc[0] / 100.0)
               for c in clusters]
node_sizes = [stats_df.loc[stats_df.cluster == c, "n_cells"].iloc[0] for c in clusters]
max_size = max(node_sizes)
node_sizes_scaled = [400 + 2000 * (s / max_size) for s in node_sizes]

sc.pl.paga(
    adata_m, ax=ax, show=False,
    color="pct_M2_cluster",
    cmap=m1m2_cmap,
    node_size_scale=2.5,
    node_size_power=0.5,
    fontsize=11,
    edge_width_scale=0.5,
    threshold=0.05,
    colorbar=False,  # we add a custom one below
)
ax.set_title("A  PAGA — sub-clusters macrofágicos\n"
             "Color: fracción M2 (azul=M1 puro, rojo=M2 puro) · tamaño: # células",
             loc="left", fontsize=11, fontweight="bold")

# Custom colorbar for M1/M2 composition
cax = fig.add_axes([0.395, 0.18, 0.012, 0.45])
sm = plt.cm.ScalarMappable(cmap=m1m2_cmap, norm=plt.Normalize(vmin=0, vmax=100))
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("% M2 cells in sub-cluster", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# ---- Panel B: Dotplot of M1/M2 markers per sub-cluster ----
ax2 = fig.add_subplot(gs[0, 1])
gene_idx = {g: i for i, g in enumerate(adata_m.var_names)}
mean_expr = np.zeros((len(clusters), len(markers_in)))
pct_expr = np.zeros((len(clusters), len(markers_in)))
for i, c in enumerate(clusters):
    cells = adata_m.obs["macro_leiden"] == c
    X_sub = adata_m[cells].X
    if hasattr(X_sub, "toarray"):
        X_sub = X_sub.toarray()
    for j, g in enumerate(markers_in):
        col = X_sub[:, gene_idx[g]]
        mean_expr[i, j] = col[col > 0].mean() if (col > 0).any() else 0.0
        pct_expr[i, j] = 100 * (col > 0).mean()

# Normalize mean expression per gene (column) to 0-1
col_max = mean_expr.max(axis=0, keepdims=True)
col_max[col_max == 0] = 1
mean_norm = mean_expr / col_max

# Reorder clusters by pct_M2 (low → high) so dotplot reads left-bottom→top-right
order = stats_df.sort_values("pct_M2").cluster.tolist()
order_idx = [clusters.index(c) for c in order]
mean_norm = mean_norm[order_idx]
pct_expr = pct_expr[order_idx]

for i in range(mean_norm.shape[0]):
    for j in range(mean_norm.shape[1]):
        size = pct_expr[i, j] * 12  # max ~1200 if 100%
        is_m1 = markers_in[j] in M1
        color_val = mean_norm[i, j]
        color = (plt.cm.Blues(color_val) if is_m1 else plt.cm.Reds(color_val))
        ax2.scatter(j, i, s=max(size, 4), c=[color], edgecolors="black",
                    linewidths=0.4)

ax2.set_yticks(range(len(order)))
m2_pcts = [f"{stats_df.loc[stats_df.cluster == c, 'pct_M2'].iloc[0]:.0f}"
           for c in order]
ax2.set_yticklabels([f"L{c} ({n} cells, {pct}% M2)"
                     for c, n, pct in zip(order,
                                           [stats_df.loc[stats_df.cluster == c, "n_cells"].iloc[0]
                                            for c in order],
                                           m2_pcts)],
                    fontsize=9)
ax2.set_xticks(range(len(markers_in)))
ax2.set_xticklabels(markers_in, fontsize=9, rotation=45, ha="right")
ax2.set_xlabel("Markers (M1: blue scale | M2: red scale)", fontsize=10)
ax2.set_title("B  Marker expression per sub-cluster\n"
              "Tamaño = % células expresando · color = expresión media (normalizada por gen)",
              loc="left", fontsize=11, fontweight="bold")

# Vertical separator between M1 and M2 markers
ax2.axvline(len(M1) - 0.5, color="gray", linestyle="--", lw=0.7)
ax2.text((len(M1) - 1) / 2, len(order) + 0.4, "M1 markers",
         ha="center", fontsize=9, color="#1f77b4", fontweight="bold")
ax2.text(len(M1) + (len(M2) - 1) / 2, len(order) + 0.4, "M2 markers",
         ha="center", fontsize=9, color="#d62728", fontweight="bold")
ax2.set_xlim(-0.5, len(markers_in) - 0.5)
ax2.set_ylim(-0.5, len(order) - 0.5 + 0.5)

plt.tight_layout()
out = FIG / "F5g_paga_annotated.png"
plt.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print(f"\nSaved: {out}")

# Also print readable cluster summary
print("\nCluster composition summary:")
print(stats_df.sort_values("pct_M2").to_string(index=False))
