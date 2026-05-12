#!/usr/bin/env python3
"""
F1 — Novae Spatial Domain Assignment.

Thin wrapper that consumes F8 outputs (embeddings + multi-level domains) and
emits a clean CSV with the user-selected primary domain level for downstream F4 (Niche-DE)
and F5 (Slingshot start cluster).

If F8 has not been run, this script can run Novae inference inline (slower path).

Citation: Blampey Q et al. Nat Methods 22, 2539-2550 (2025).
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--pretrained", default="MICS-Lab/novae-human-0")
    p.add_argument("--n-domains-levels", type=int, nargs="+", default=[5, 10, 20])
    p.add_argument("--primary-level", type=int, default=10,
                   help="Which level to report as primary `novae_domain` column")
    p.add_argument("--radius", type=float, default=200.0)
    p.add_argument("--f8-dir", type=Path, default=None,
                   help="Optional F8 output dir; if present, reuse embeddings/domains")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log.info("=" * 70)
    log.info("F1 — Novae Spatial Domains")
    log.info("=" * 70)

    # Try to reuse F8 output first (fast path)
    f8_domains_csv = None
    if args.f8_dir:
        f8_domains_csv = args.f8_dir / "novae_domains_multilevel.csv"
    else:
        # Auto-detect F8 output relative to this rule's expected location
        candidate = args.outdir.parent / "F8_novae" / "novae_domains_multilevel.csv"
        if candidate.exists():
            f8_domains_csv = candidate

    if f8_domains_csv and f8_domains_csv.exists():
        log.info(f"[fast path] Reusing F8 output: {f8_domains_csv}")
        domains_df = pd.read_csv(f8_domains_csv, index_col=0)
        primary_col = f"novae_domains_L{args.primary_level}"
        if primary_col not in domains_df.columns:
            available = list(domains_df.columns)
            raise RuntimeError(
                f"F8 output missing column {primary_col}. Available: {available}"
            )
    else:
        log.info("[slow path] F8 output not found; running Novae inference inline")
        import spatialdata as sd
        import novae
        sdata = sd.read_zarr(str(args.sdata))
        adata = sdata.tables[list(sdata.tables.keys())[0]]
        log.info(f"  AnnData: {adata.shape}")
        novae.spatial_neighbors(adata, technology="xenium", radius=args.radius)
        model = novae.Novae.from_pretrained(args.pretrained)
        model.compute_representations(adata, zero_shot=True, accelerator="gpu", num_workers=4)
        domain_cols = []
        for level in args.n_domains_levels:
            col = f"novae_domains_L{level}"
            model.assign_domains(adata, level=level, key_added=col)
            domain_cols.append(col)
        domains_df = adata.obs[domain_cols].copy()
        domains_df.index.name = "cell_id"
        primary_col = f"novae_domains_L{args.primary_level}"

    # Add primary column alias
    domains_df["novae_domain"] = domains_df[primary_col]
    out_path = args.outdir / "novae_domains.csv"
    domains_df.to_csv(out_path)
    log.info(f"  Saved: {out_path} (primary level={args.primary_level})")

    # Quick stats
    n_domains = domains_df["novae_domain"].nunique()
    log.info(f"\nDomain counts (primary level L{args.primary_level}: {n_domains} unique):")
    log.info(domains_df["novae_domain"].value_counts().head(15).to_string())

    summary = {
        "primary_level": args.primary_level,
        "n_unique_domains_primary": int(n_domains),
        "n_cells": len(domains_df),
        "execution_seconds": round(time.time() - t0, 1),
    }
    (args.outdir / "F1_novae_summary.json").write_text(json.dumps(summary, indent=2))

    log.info(f"\n✓ F1 (Novae) complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
