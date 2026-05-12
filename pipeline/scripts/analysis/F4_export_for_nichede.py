#!/usr/bin/env python3
"""
F4 — Export data for NicheDE analysis.

Exports from SpatialData:
  - counts_matrix.mtx  : genes × cells raw counts (Matrix Market)
  - gene_names.csv     : gene names
  - cell_names.csv     : cell IDs
  - spatial_coords.csv : x, y per cell
  - cell_types.csv     : cell_id, cell_type (L1 broad), cell_type_fine (L2)
  - banksy_domains.csv : cell_id, banksy_domain (14 Banksy spatial domains)

NicheDE requires raw counts + spatial coordinates + cell type labels.
Domains (from Phase B Banksy) serve as the niche context.

Citation: Mason JC et al. Genome Biol 25, 107 (2024).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import spatialdata as sd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata",   required=True, type=Path)
    p.add_argument("--domains", required=True, type=Path,
                   help="banksy_domains.csv from Phase B (F1_spatial_domains/)")
    p.add_argument("--outdir",  required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 70)
    log.info("F4 — Export for NicheDE")
    log.info("=" * 70)

    log.info(f"Loading SpatialData from {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  {adata.n_obs:,} cells × {adata.n_vars} genes")

    # ── Counts matrix ────────────────────────────────────────────────────────
    log.info("\n[1/4] Exporting counts matrix...")
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    X_gc = X.T.tocsc() if sp.issparse(X) else sp.csc_matrix(X.T)  # genes × cells

    scipy.io.mmwrite(str(args.outdir / "counts_matrix.mtx"), X_gc)
    pd.Series(adata.var_names.tolist()).to_csv(
        args.outdir / "gene_names.csv", index=False, header=False)
    pd.Series(adata.obs_names.tolist()).to_csv(
        args.outdir / "cell_names.csv", index=False, header=False)
    log.info(f"  Counts shape: {X_gc.shape} (genes × cells)")

    # ── Spatial coordinates ──────────────────────────────────────────────────
    log.info("\n[2/4] Exporting spatial coordinates...")
    coords = adata.obsm["spatial"]
    coords_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "x": coords[:, 0],
        "y": coords[:, 1],
    })
    coords_df.to_csv(args.outdir / "spatial_coords.csv", index=False)
    log.info(f"  Coords: {coords_df.shape[0]:,} cells")

    # ── Cell type labels ─────────────────────────────────────────────────────
    log.info("\n[3/4] Exporting cell type labels...")
    ct_df = adata.obs[["cell_type", "cell_type_fine"]].copy()
    ct_df.index.name = "cell_id"
    ct_df = ct_df.reset_index()

    # Summary
    log.info(f"  cell_type (L1) unique: {ct_df['cell_type'].nunique()}")
    for ct, n in ct_df["cell_type"].value_counts().items():
        log.info(f"    {ct:25s} {n:7,}")
    log.info(f"  cell_type_fine (L2) unique: {ct_df['cell_type_fine'].nunique()}")

    ct_df.to_csv(args.outdir / "cell_types.csv", index=False)

    # ── Banksy domain labels ─────────────────────────────────────────────────
    log.info("\n[4/4] Merging Banksy domain labels...")
    domains_df = pd.read_csv(args.domains)
    # Align cell_id types (adata obs_names may be strings "0","1",... while CSV may have int)
    domains_df["cell_id"] = domains_df["cell_id"].astype(str)
    # Align to adata cell order
    domain_merged = ct_df[["cell_id"]].merge(domains_df, on="cell_id", how="left")
    domain_merged["banksy_domain"] = domain_merged["banksy_domain"].fillna("Unknown")

    n_with_domain = (domain_merged["banksy_domain"] != "Unknown").sum()
    log.info(f"  Cells with Banksy domain: {n_with_domain:,}/{len(domain_merged):,}")
    log.info(f"  Domain breakdown:")
    for d, n in domain_merged["banksy_domain"].value_counts().items():
        log.info(f"    {d:8s} {n:7,}")

    domain_merged.to_csv(args.outdir / "cell_domains.csv", index=False)

    log.info(f"\n✓ Export complete → {args.outdir}")


if __name__ == "__main__":
    main()
