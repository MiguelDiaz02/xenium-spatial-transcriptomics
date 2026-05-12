#!/usr/bin/env python3
"""
F8 — Novae Foundation Model: embeddings, hierarchical spatial domains, pathway analysis.

This is the SUBSTRATE for F1 (Novae domain consensus with Banksy) and F5 (Slingshot
on Novae embeddings). Producing it once and persisting it back into sdata.zarr
makes F1/F5 fast and reproducible.

Pipeline (per Blampey et al., Nat Methods 22, 2539-2550, 2025):
  1. Build spatial graph via novae.spatial_neighbors(technology='xenium')
  2. Load pretrained model MICS-Lab/novae-human-0
  3. zero_shot=True compute_representations on full dataset
  4. Multi-level hierarchical domain assignment (levels [5, 10, 20])
  5. Pathway scoring with curated lung-cancer-relevant Hallmark gene sets
  6. Spatially variable gene analysis per domain
  7. Persist embeddings (.parquet), domains (.csv), pathway scores, SVG ranks

Usage:
  python F8_novae_embeddings.py \\
      --sdata human_lung_cancer/results/sdata.zarr \\
      --pretrained MICS-Lab/novae-human-0 \\
      --n-levels 5,10,20 \\
      --radius 200 \\
      --pathway-db msigdb_hallmark \\
      --batch-correct \\
      --outdir results/03_phase3_spatial/F8_novae

Citation:
  Blampey Q, Benkirane H, Bercovici N et al. Novae: a graph-based foundation model
  for spatial transcriptomics data. Nat Methods 22, 2539-2550 (2025).
  doi:10.1038/s41592-025-02899-6
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

import torch
import spatialdata as sd
import anndata as ad
import novae

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


# ─── Curated pathway gene sets for 289-gene Xenium INMUNO panel ─────────────
# Subset of MSigDB Hallmark pathways with ≥4 genes overlapping the panel.
# Selection rationale: cancer immunology and tumor microenvironment biology.
HALLMARK_PATHWAYS_CURATED = {
    "INTERFERON_GAMMA_RESPONSE": [
        "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "GBP1", "GBP4", "GBP5",
        "PSMB8", "PSMB9", "TAP1", "CD274", "HLA-DRA", "CIITA",
    ],
    "INTERFERON_ALPHA_RESPONSE": [
        "ISG15", "IFI44", "IFIT1", "IFIT3", "IFI6", "OAS1", "MX1", "STAT1",
    ],
    "INFLAMMATORY_RESPONSE": [
        "IL6", "TNF", "IL1B", "CCL2", "CCL3", "CCL4", "CCL5", "CXCL8",
        "CXCL1", "CXCL2", "PTGS2", "NLRP3", "NFKB1", "NFKB2",
    ],
    "TNF_NF_KB_SIGNALING": [
        "TNF", "TNFAIP3", "NFKB1", "NFKB2", "RELB", "ICAM1", "IL1B", "CCL2",
        "CCL5", "CXCL1", "CXCL2", "CXCL10",
    ],
    "EPITHELIAL_MESENCHYMAL_TRANSITION": [
        "VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
        "TGFB1", "TGFB2", "MMP2", "MMP9", "ITGAV", "ITGB1",
    ],
    "HYPOXIA": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "PDK1", "ENO1", "PGK1", "ALDOA",
    ],
    "ANGIOGENESIS": [
        "VEGFA", "VEGFC", "PECAM1", "PDGFRA", "PDGFRB", "FLT1", "TEK", "KDR",
        "VWF", "CDH5",
    ],
    "TGF_BETA_SIGNALING": [
        "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3",
        "SMAD7", "SERPINE1",
    ],
    "IL2_STAT5_SIGNALING": [
        "IL2RA", "IL2RB", "IL2RG", "STAT5A", "STAT5B", "FOXP3", "TNFRSF18",
        "TNFRSF8", "TNFRSF4",
    ],
    "IL6_JAK_STAT3_SIGNALING": [
        "IL6", "IL6R", "STAT3", "SOCS3", "JAK1", "JAK2", "OSMR", "TYK2",
    ],
    "ALLOGRAFT_REJECTION": [
        "CD2", "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "CD28", "ICOS",
        "IL2", "IFNG", "GZMA", "GZMB", "PRF1", "FASLG", "CCR5",
    ],
    "COMPLEMENT": [
        "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B",
        "C5", "C5AR1", "CFB", "CFD",
    ],
    "APOPTOSIS": [
        "BCL2", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS", "FASLG",
        "BIRC5", "XIAP", "TNFSF10",
    ],
    "MYC_TARGETS_V1": [
        "MYC", "MAX", "ODC1", "EIF4E", "EIF4G1", "PTMA", "MCM2", "MCM5",
    ],
    "G2M_CHECKPOINT": [
        "CCNB1", "CCNB2", "CDK1", "CDC20", "AURKA", "AURKB", "PLK1", "TOP2A",
        "MKI67", "TPX2", "CENPF",
    ],
    "P53_PATHWAY": [
        "TP53", "MDM2", "MDM4", "CDKN1A", "BAX", "PUMA", "GADD45A", "PML",
    ],
    "KRAS_SIGNALING_UP": [
        "KRAS", "EGFR", "DUSP6", "ETV1", "ETV4", "ETV5", "SPRY2", "SPRY4",
    ],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--pretrained", default="MICS-Lab/novae-human-0")
    p.add_argument("--n-levels", default="5,10,20",
                   help="Comma-separated hierarchical levels (e.g., '5,10,20')")
    p.add_argument("--radius", type=float, default=200.0,
                   help="Spatial neighbor radius in microns")
    p.add_argument("--pathway-db", default="msigdb_hallmark",
                   help="Pathway database (currently only 'msigdb_hallmark' supported)")
    p.add_argument("--batch-correct", action="store_true",
                   help="Run native Novae batch correction (multi-sample only)")
    p.add_argument("--accelerator", default="gpu", choices=["gpu", "cpu"])
    p.add_argument("--num-workers", type=int, default=4)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--write-back", action="store_true",
                   help="Persist novae_latent + domain columns back into sdata.zarr")
    return p.parse_args()


def filter_pathways_to_panel(pathways: dict[str, list[str]],
                             panel: set[str],
                             min_size: int = 4) -> dict[str, list[str]]:
    """Keep only pathways with at least `min_size` genes in the panel."""
    kept = {}
    for name, genes in pathways.items():
        present = [g for g in genes if g in panel]
        if len(present) >= min_size:
            kept[name] = present
    return kept


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    levels = [int(x) for x in args.n_levels.split(",")]
    t0 = time.time()

    log.info("=" * 70)
    log.info("F8 — Novae foundation model (embeddings + multi-level domains)")
    log.info("=" * 70)
    log.info(f"  Pretrained: {args.pretrained}")
    log.info(f"  Hierarchical levels: {levels}")
    log.info(f"  Spatial neighbor radius: {args.radius} μm")

    log.info(f"\n[1/7] Loading sdata.zarr from {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"      AnnData: {adata.shape}")
    log.info(f"      Spatial coords present: {('spatial' in adata.obsm)}")

    log.info("\n[2/7] Building Novae spatial neighbors (using existing obsm['spatial'])...")
    t = time.time()
    # adata.obsm['spatial'] already populated by upstream Xenium ingest;
    # `technology` would re-compute coords (conflicts), so use radius-only.
    novae.spatial_neighbors(adata, radius=args.radius)
    log.info(f"      done ({time.time() - t:.1f}s); novae_sid: {'novae_sid' in adata.obs}")

    log.info(f"\n[3/7] Loading pretrained model {args.pretrained}...")
    t = time.time()
    model = novae.Novae.from_pretrained(args.pretrained)
    log.info(f"      Model loaded ({time.time() - t:.1f}s)")
    log.info(f"      Embedding dim: {model.hparams.get('embedding_size')}")
    log.info(f"      Num prototypes: {model.hparams.get('num_prototypes')}")

    log.info("\n[4/7] Computing representations (zero-shot)...")
    t = time.time()
    model.compute_representations(adata, zero_shot=True,
                                   accelerator=args.accelerator,
                                   num_workers=args.num_workers)
    log.info(f"      done ({time.time() - t:.1f}s)")
    if "novae_latent" not in adata.obsm:
        raise RuntimeError("Novae compute_representations did not produce 'novae_latent'")
    log.info(f"      novae_latent shape: {adata.obsm['novae_latent'].shape}")
    if args.accelerator == "gpu":
        log.info(f"      GPU memory: {torch.cuda.memory_allocated()/1e9:.2f} GB allocated")

    if args.batch_correct:
        log.info("\n[4b/7] Native batch effect correction... (skip if single sample)")
        slide_keys = [c for c in adata.obs.columns if "slide" in c.lower() or "sample" in c.lower()]
        if len(slide_keys) > 0:
            try:
                model.batch_effect_correction(adata, obs_key=slide_keys[0])
                log.info(f"      done with key '{slide_keys[0]}'")
            except Exception as e:
                log.warning(f"      batch_effect_correction skipped: {e}")
        else:
            log.info("      No sample/slide column detected; skipping")

    log.info("\n[5/7] Multi-level hierarchical domain assignment...")
    domain_cols = []
    for level in levels:
        t = time.time()
        col = f"novae_domains_L{level}"
        model.assign_domains(adata, level=level, key_added=col)
        domain_cols.append(col)
        n_distinct = adata.obs[col].nunique()
        log.info(f"      level={level}: {col} produced {n_distinct} domains ({time.time() - t:.1f}s)")
        log.info(f"         distribution (top 5):\n{adata.obs[col].value_counts().head(5).to_string()}")

    log.info("\n[6/7] Pathway scoring (curated Hallmark)...")
    panel = set(adata.var_names)
    pathways = filter_pathways_to_panel(HALLMARK_PATHWAYS_CURATED, panel, min_size=4)
    log.info(f"      Hallmark pathways with ≥4 genes in panel: {len(pathways)}/{len(HALLMARK_PATHWAYS_CURATED)}")
    pathway_df = None
    if len(pathways) > 0:
        # Use the highest-resolution domain key as obs_key for plotting; pathway_scores returns df
        max_level_col = domain_cols[-1]
        try:
            pathway_df = novae.plot.pathway_scores(
                adata, pathways, obs_key=max_level_col,
                min_pathway_size=4, return_df=True, show=False,
            )
            if pathway_df is not None:
                log.info(f"      pathway_scores df shape: {pathway_df.shape}")
        except Exception as e:
            log.warning(f"      pathway_scores failed: {e}")

    log.info("\n[7/7] Spatially variable genes per domain...")
    svg_results = {}
    try:
        # spatially_variable_genes returns list of top SVGs across domains
        max_level_col = domain_cols[-1]
        top_svgs = novae.plot.spatially_variable_genes(
            adata, obs_key=max_level_col, top_k=10,
            min_positive_ratio=0.05, return_list=True, show=False,
        )
        if top_svgs is not None:
            log.info(f"      Top SVGs (level={max_level_col}): {top_svgs[:10]}")
            svg_results["top_svgs_max_level"] = top_svgs
    except Exception as e:
        log.warning(f"      spatially_variable_genes failed: {e}")

    # ── Persist outputs ───────────────────────────────────────────────────────
    log.info("\nWriting outputs...")
    embeddings_df = pd.DataFrame(
        adata.obsm["novae_latent"],
        index=adata.obs_names,
        columns=[f"novae_dim_{i}" for i in range(adata.obsm["novae_latent"].shape[1])],
    )
    embeddings_path = args.outdir / "novae_latent.parquet"
    embeddings_df.to_parquet(embeddings_path)
    log.info(f"  Saved: {embeddings_path}")

    domains_df = adata.obs[domain_cols].copy()
    domains_df.index.name = "cell_id"
    domains_path = args.outdir / "novae_domains_multilevel.csv"
    domains_df.to_csv(domains_path)
    log.info(f"  Saved: {domains_path}")

    if pathway_df is not None:
        pathway_path = args.outdir / "novae_pathway_scores.csv"
        pathway_df.to_csv(pathway_path)
        log.info(f"  Saved: {pathway_path}")

    if svg_results:
        svg_path = args.outdir / "novae_svg_scores.csv"
        # Convert top_svgs list to a DataFrame
        top_svgs_list = svg_results.get("top_svgs_max_level", [])
        pd.DataFrame({"gene": top_svgs_list, "rank": range(1, len(top_svgs_list) + 1)}).to_csv(
            svg_path, index=False,
        )
        log.info(f"  Saved: {svg_path}")

    # Summary
    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "embedding_dim": int(adata.obsm["novae_latent"].shape[1]),
        "domain_levels": levels,
        "domains_per_level": {col: int(adata.obs[col].nunique()) for col in domain_cols},
        "pathways_kept": len(pathways) if pathways else 0,
        "execution_seconds": round(time.time() - t0, 1),
    }
    (args.outdir / "F8_summary.json").write_text(json.dumps(summary, indent=2))
    log.info(f"  Saved: {args.outdir / 'F8_summary.json'}")

    if args.write_back:
        log.info("\nPersisting novae outputs back into zarr...")
        sdata.tables[table_name] = adata
        sdata.write_element(table_name, element_type="tables", overwrite=True) \
            if hasattr(sdata, "write_element") else None
        log.info("  Done.")

    log.info(f"\n✓ F8 complete in {time.time() - t0:.1f}s; outputs: {args.outdir}")


if __name__ == "__main__":
    main()
