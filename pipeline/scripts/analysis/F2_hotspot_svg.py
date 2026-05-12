#!/usr/bin/env python3
"""
F2 — Hotspot: Automatic Expression Histology Analysis

Identifies Spatially Variable Genes (SVGs) using information-theoretic approach.
Hotspot is top-ranked in Salas et al. 2025 benchmark (99.2% TPR, 0.8% FPR).

Citation: DeTomaso D, Yosef N. Cell Syst 12, 446-456 (2021).
Benchmarking: Salas et al. Nat Methods 22, 813-823 (2025).
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

import spatialdata as sd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--n_neighbors", type=int, default=30)
    p.add_argument("--model", default="danb",
                   help="Hotspot model: danb | bernoulli | normal | none")
    p.add_argument("--fdr_threshold", type=float, default=0.05)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log.info("=" * 70)
    log.info("F2 — Hotspot SVG Detection")
    log.info("=" * 70)

    log.info(f"Loading SpatialData from {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  AnnData: {adata.shape[0]:,} cells × {adata.shape[1]} genes")

    if "spatial" not in adata.obsm:
        raise RuntimeError("Missing 'spatial' in obsm")
    log.info(f"  Spatial coords present: {adata.obsm['spatial'].shape}")

    # Hotspot needs raw counts in layers["counts"] and spatial coords in obsm
    # Use counts layer if available, else fall back to X
    if "counts" in adata.layers:
        adata_hs = adata.copy()
        adata_hs.X = adata.layers["counts"].copy()
        log.info("  Using layers['counts'] as expression input")
    else:
        adata_hs = adata.copy()
        log.info("  Using X as expression input (no counts layer found)")

    try:
        import hotspot
    except ImportError:
        log.error("hotspot not installed. Run: pip install hotspot-sc")
        sys.exit(1)

    log.info(f"\n[1/3] Initializing Hotspot (n_neighbors={args.n_neighbors}, model={args.model})...")
    t = time.time()

    # Hotspot API: Hotspot(adata, latent_obsm_key=..., model=...)
    hs = hotspot.Hotspot(
        adata_hs,
        latent_obsm_key="spatial",
        model=args.model,
    )
    log.info(f"  Initialized ({time.time() - t:.1f}s)")

    log.info(f"\n[2/3] Creating KNN graph...")
    t = time.time()
    hs.create_knn_graph(weighted_graph=False, n_neighbors=args.n_neighbors)
    log.info(f"  KNN graph ready ({time.time() - t:.1f}s)")

    log.info(f"\n[3/3] Computing local autocorrelations (this is the slow step)...")
    t = time.time()
    hs_results = hs.compute_autocorrelations(jobs=1)
    log.info(f"  Done ({time.time() - t:.1f}s)")

    # hs_results is a DataFrame indexed by Gene with columns: C, Z, Pval, FDR
    df = hs_results.copy()
    df.index.name = "gene"
    df = df.reset_index()   # gene, C, Z, Pval, FDR
    df.rename(columns={"Gene": "gene"}, inplace=True)  # handle capital-G index name

    # Hotspot already returns BH-corrected FDR; use it directly
    df["significant"] = df["FDR"] < args.fdr_threshold
    df = df.sort_values("C", ascending=False).reset_index(drop=True)

    out_path = args.outdir / "hotspot_svg_scores.csv"
    df.to_csv(out_path, index=False)
    log.info(f"  Saved: {out_path}")

    n_sig = int(df["significant"].sum())
    log.info(f"\n  Hotspot SVGs (FDR<{args.fdr_threshold}): {n_sig}/{len(df)}")
    log.info("  Top 10 SVGs:")
    gene_col = "gene" if "gene" in df.columns else df.columns[0]
    for _, row in df.head(10).iterrows():
        log.info(f"    {row[gene_col]:15s} C={row['C']:7.4f} Z={row['Z']:7.2f} FDR={row['FDR']:.2e}")

    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_neighbors": int(args.n_neighbors),
        "model": args.model,
        "fdr_threshold": float(args.fdr_threshold),
        "n_svg_significant": n_sig,
        "execution_seconds": round(time.time() - t0, 1),
    }
    (args.outdir / "hotspot_summary.json").write_text(json.dumps(summary, indent=2))
    log.info(f"\n✓ Hotspot complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
