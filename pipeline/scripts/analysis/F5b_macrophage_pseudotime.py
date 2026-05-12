#!/usr/bin/env python3
"""
F5b — Macrophage M1↔M2 pseudotime (audit-fixed sub-trajectory)

Replaces the F5 tumor→multi-end Slingshot abuse with a biologically valid
pseudotime restricted to a single, plausible developmental axis: macrophage
M1↔M2 polarization continuum (or transitional state if no clear directionality).

Methods (≥2 per family, per consensus framework):
  - PAGA + DPT (this script, scanpy-native)
  - Slingshot (F5c R script, separate)

Input:  sdata.zarr (cell_type_L2 column with Macrophage_M1, Macrophage_M2)
Output: F5b_macrophage_pseudotime/
  - macrophage_pseudotime.csv (per-cell pseudotime + metadata)
  - macrophage_counts_for_R.csv (for downstream Slingshot)
  - macrophage_summary.json
  - figures (UMAP, PAGA, DPT)
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import spatialdata as sd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
ZARR = BASE / "human_lung_cancer/results/sdata.zarr"
OUT = BASE / "human_lung_cancer/results/03_phase3_spatial/F5b_macrophage_pseudotime"
FIG = BASE / "human_lung_cancer/results/figures/phase5_pseudotime_v2"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG)
sc.settings.set_figure_params(dpi=200, figsize=(6, 5), fontsize=10)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    log("=" * 60)
    log("F5b — Macrophage M1↔M2 pseudotime (audit-fixed)")
    log("=" * 60)

    log("Loading sdata.zarr ...")
    sdata = sd.read_zarr(str(ZARR))
    table_key = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_key]
    log(f"  Total: {adata.n_obs:,} cells × {adata.n_vars} genes")

    # Subset M1 + M2 macrophages
    mask = adata.obs["cell_type_L2"].isin(["Macrophage_M1", "Macrophage_M2"])
    adata_m = adata[mask].copy()
    adata_m.obs["cell_type_L2"] = adata_m.obs["cell_type_L2"].cat.remove_unused_categories()
    log(f"  M1+M2 subset: {adata_m.n_obs:,} cells")
    log(f"  Type counts: {dict(adata_m.obs['cell_type_L2'].value_counts())}")

    # PCA on log-normalized expression (use real counts, not embeddings)
    log("\n[1/5] PCA + neighbors + UMAP ...")
    if "X_pca" in adata_m.obsm:
        del adata_m.obsm["X_pca"]
    sc.tl.pca(adata_m, n_comps=30, random_state=42)
    sc.pp.neighbors(adata_m, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.umap(adata_m, random_state=42)

    # Higher-resolution Leiden within macrophages
    log("[2/5] Sub-clustering (Leiden) ...")
    sc.tl.leiden(adata_m, resolution=0.5, random_state=42, key_added="macro_leiden")
    log(f"  Sub-clusters: {adata_m.obs['macro_leiden'].nunique()}")

    # PAGA on Leiden sub-clusters (graph between cluster connectivities)
    log("[3/5] PAGA ...")
    sc.tl.paga(adata_m, groups="macro_leiden")
    sc.pl.paga(adata_m, plot=False)  # required for downstream UMAP-init
    sc.tl.umap(adata_m, init_pos="paga", random_state=42)

    # DPT — diffusion pseudotime
    # Root: cell with strongest M1 marker score (CD86 + IL1B + TNF if in panel,
    # otherwise pick an M1-labeled cell with highest within-cluster marker)
    log("[4/5] Diffusion pseudotime (DPT) — root from M1 cluster ...")

    # Pick a root cell: M1-labeled, near M1 cluster centroid in PCA
    m1_cells = adata_m.obs.index[adata_m.obs["cell_type_L2"] == "Macrophage_M1"]
    m1_pca = adata_m.obsm["X_pca"][adata_m.obs["cell_type_L2"] == "Macrophage_M1"]
    m1_centroid = m1_pca.mean(axis=0)
    dists = np.linalg.norm(m1_pca - m1_centroid, axis=1)
    root_idx_local = int(np.argmin(dists))
    root_cell = m1_cells[root_idx_local]
    adata_m.uns["iroot"] = int(np.where(adata_m.obs.index == root_cell)[0][0])
    log(f"  Root cell: {root_cell} (M1, PCA-centroid)")

    sc.tl.diffmap(adata_m, n_comps=15, random_state=42)
    sc.tl.dpt(adata_m, n_branchings=0)

    # Save figures
    log("[5/5] Figures ...")
    sc.pl.umap(adata_m, color=["cell_type_L2", "macro_leiden", "dpt_pseudotime"],
               save="_F5b_macrophage_overview.png", show=False, ncols=3, wspace=0.3)
    sc.pl.paga(adata_m, color="macro_leiden", save="_F5b_macrophage_paga.png", show=False)

    # Pseudotime distribution by cell type
    fig, ax = plt.subplots(figsize=(6, 4))
    for ct, color in zip(["Macrophage_M1", "Macrophage_M2"], ["#1f77b4", "#d62728"]):
        ax.hist(adata_m.obs.loc[adata_m.obs["cell_type_L2"] == ct, "dpt_pseudotime"],
                bins=40, alpha=0.6, label=ct, color=color)
    ax.set_xlabel("DPT pseudotime")
    ax.set_ylabel("Cell count")
    ax.set_title("DPT pseudotime distribution by macrophage state")
    ax.legend()
    plt.tight_layout()
    plt.savefig(FIG / "F5b_dpt_distribution.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Save per-cell pseudotime
    out_df = adata_m.obs[["cell_type_L2", "macro_leiden", "dpt_pseudotime"]].copy()
    out_df["umap_1"] = adata_m.obsm["X_umap"][:, 0]
    out_df["umap_2"] = adata_m.obsm["X_umap"][:, 1]
    if "spatial" in adata_m.obsm:
        out_df["x"] = adata_m.obsm["spatial"][:, 0]
        out_df["y"] = adata_m.obsm["spatial"][:, 1]
    out_df.to_csv(OUT / "macrophage_pseudotime.csv")
    log(f"  Saved: {OUT / 'macrophage_pseudotime.csv'}")

    # Export counts + metadata for R Slingshot
    log("\nExporting counts for R Slingshot ...")
    if hasattr(adata_m.X, "toarray"):
        counts = adata_m.X.toarray()
    else:
        counts = adata_m.X
    # Use raw counts if available
    if "counts" in adata_m.layers:
        counts_raw = adata_m.layers["counts"]
        if hasattr(counts_raw, "toarray"):
            counts_raw = counts_raw.toarray()
    else:
        # Fall back to current X (log-normalized — flag in summary)
        counts_raw = counts
        log("  WARNING: no 'counts' layer; using adata.X (likely log-norm) — Slingshot+tradeSeq need raw counts")

    counts_df = pd.DataFrame(counts_raw, index=adata_m.obs_names, columns=adata_m.var_names)
    counts_df.to_csv(OUT / "macrophage_counts_for_R.csv")
    log(f"  Saved counts: {counts_df.shape}")

    meta_df = adata_m.obs[["cell_type_L2", "macro_leiden", "dpt_pseudotime"]].copy()
    meta_df["pca_1"] = adata_m.obsm["X_pca"][:, 0]
    meta_df["pca_2"] = adata_m.obsm["X_pca"][:, 1]
    meta_df["umap_1"] = adata_m.obsm["X_umap"][:, 0]
    meta_df["umap_2"] = adata_m.obsm["X_umap"][:, 1]
    meta_df.to_csv(OUT / "macrophage_metadata_for_R.csv")

    # Save full PCA for R
    pca_df = pd.DataFrame(adata_m.obsm["X_pca"], index=adata_m.obs_names,
                          columns=[f"PC{i+1}" for i in range(adata_m.obsm["X_pca"].shape[1])])
    pca_df.to_csv(OUT / "macrophage_pca_for_R.csv")

    summary = {
        "n_cells_total": int(adata.n_obs),
        "n_cells_macrophage": int(adata_m.n_obs),
        "n_M1": int((adata_m.obs["cell_type_L2"] == "Macrophage_M1").sum()),
        "n_M2": int((adata_m.obs["cell_type_L2"] == "Macrophage_M2").sum()),
        "n_subclusters_leiden": int(adata_m.obs["macro_leiden"].nunique()),
        "root_cell": str(root_cell),
        "root_cell_type": "Macrophage_M1",
        "method_1": "PAGA+DPT (scanpy)",
        "method_2_pending": "Slingshot (F5c R script)",
        "raw_counts_used": "counts" in adata_m.layers,
    }
    with open(OUT / "macrophage_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log("\n=== F5b COMPLETE ===")
    log(f"Outputs: {OUT}")
    log(f"Figures: {FIG}")


if __name__ == "__main__":
    main()
