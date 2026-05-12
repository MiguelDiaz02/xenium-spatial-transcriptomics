#!/usr/bin/env python3
"""
F0 reannotation v3 — Anotación híbrida anclada en Leiden + immune_granular

Estrategia:
  1. Células INMUNES → cell_type_immune_granular (subclustering computacional real)
     Granularidad adicional (CD8_exhausted vs cytotoxic, B_naive vs memory,
     Plasma vs Plasmablast, cDC1/cDC2/pDC/DC_mature) resuelta por score_genes.
  2. Células NO-INMUNES → Leiden cluster + marker score per-célula
     Clusters puros (≥70%): asignación directa.
     Clusters mixtos: argmax de score_genes dentro del contexto del cluster.
  3. Todo mapeado al esquema L1/L2/L3 existente.
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import spatialdata as sd
import scanpy as sc
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
from utils.paths import sdata_path, results_root, project_root
log = get_logger(__name__)

DEFAULT_MARKERS_YAML = project_root() / "pipeline" / "config" / "markers" / "lung_pilot.yaml"


def load_markers_config(path: Path) -> dict:
    """Load the markers + hierarchy + annotation parameters from YAML."""
    with open(path) as f:
        cfg = yaml.safe_load(f)
    required = ["markers", "hierarchy_l1", "hierarchy_l3",
                "immune_granular_to_l2", "nonimmune_l2", "annotation_params"]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"{path}: missing required keys {missing}")
    return cfg


def pick_best_score(obs: pd.DataFrame, candidates: list[str]) -> pd.Series:
    """Return the candidate L2 type with highest z-score per cell."""
    score_cols = [f"score_{c}" for c in candidates if f"score_{c}" in obs.columns]
    if not score_cols:
        return pd.Series(candidates[0], index=obs.index)
    return obs[score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)


def main():
    parser = argparse.ArgumentParser(description="F0 reannotation v3 — hybrid annotation pipeline")
    parser.add_argument("--markers", type=Path, default=DEFAULT_MARKERS_YAML,
                        help=f"Path to markers YAML config (default: {DEFAULT_MARKERS_YAML})")
    parser.add_argument("--sdata", type=Path, default=sdata_path(),
                        help="Path to sdata.zarr (default: results/sdata.zarr)")
    parser.add_argument("--outdir", type=Path,
                        default=results_root() / "02_biology/reannotation_v3",
                        help="Output directory")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log.info("=" * 70)
    log.info("F0 reannotation v3 — Leiden-anchor + immune_granular hybrid")
    log.info("=" * 70)
    log.info(f"  markers config: {args.markers}")
    log.info(f"  sdata path:     {args.sdata}")
    log.info(f"  output dir:     {args.outdir}")

    # ── Load markers config ───────────────────────────────────────────────
    cfg = load_markers_config(args.markers)
    markers = cfg["markers"]
    L2_TO_L1 = cfg["hierarchy_l1"]
    L2_TO_L3 = cfg["hierarchy_l3"]
    immune_to_l2 = cfg["immune_granular_to_l2"]
    nonimmune_l2 = cfg["nonimmune_l2"]
    leiden_direct = cfg.get("leiden_direct", {}) or {}
    params = cfg["annotation_params"]
    immune_col = params.get("immune_granular_col", "cell_type_immune_granular")
    leiden_col = params.get("leiden_col", "leiden")
    score_threshold = float(params.get("score_threshold", 0.3))
    use_leiden_direct = bool(params.get("use_leiden_direct", True))
    unassigned = params.get("unassigned_label", "Unassigned")
    min_markers = int(params.get("min_markers_per_signature", 2))
    log.info(f"  loaded {len(markers)} L2 marker signatures from YAML")

    # ── Load ─────────────────────────────────────────────────────────────
    log.info("\n[1/6] Cargando sdata.zarr...")
    sdata = sd.read_zarr(str(args.sdata))
    table_key = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_key]
    obs = adata.obs.copy()
    n = len(obs)
    log.info(f"  {n:,} células × {adata.n_vars} genes")

    # ── Compute missing score columns if needed ───────────────────────────
    log.info("\n[2/6] Verificando score columns...")
    adata_work = adata.copy()
    missing = [ct for ct in markers if f"score_{ct}" not in obs.columns]
    if missing:
        log.info(f"  Calculando scores para: {missing}")
        for ct in missing:
            genes = [g for g in markers[ct] if g in adata_work.var_names]
            if len(genes) >= min_markers:
                sc.tl.score_genes(adata_work, gene_list=genes, score_name=f"score_{ct}", use_raw=False)
                obs[f"score_{ct}"] = adata_work.obs[f"score_{ct}"].values
                log.info(f"    ✓ {ct} ({len(genes)} genes)")
            else:
                log.warning(f"    ✗ {ct} skipped: only {len(genes)} markers in panel (min {min_markers})")
    else:
        log.info(f"  Todos los score columns presentes ({len(markers)} tipos)")

    # ── Build v3 annotation ───────────────────────────────────────────────
    log.info("\n[3/6] Construyendo anotación v3...")
    l2 = pd.Series(unassigned, index=obs.index, name="cell_type_L2_v3")
    source = pd.Series("score_only", index=obs.index, name="annotation_source")
    leiden = obs[leiden_col].astype(str)
    ig = obs[immune_col].fillna("Non-immune").astype(str)

    # Step A: immune_granular → fine L2 (sub-typing by score)
    n_immune = 0
    for ig_type, l2_candidates in immune_to_l2.items():
        mask = ig == ig_type
        if mask.sum() == 0:
            continue
        if len(l2_candidates) == 1:
            l2[mask] = l2_candidates[0]
        else:
            l2[mask] = pick_best_score(obs[mask], l2_candidates)
        source[mask] = "immune_granular"
        n_immune += mask.sum()
        log.info(f"  [immune] {ig_type:25} → {l2[mask].value_counts().to_dict()}  (n={mask.sum():,})")

    # Step B: pure Leiden clusters (non-immune)
    n_leiden_direct = 0
    if use_leiden_direct:
        for leiden_id, l2_type in leiden_direct.items():
            mask = (leiden == str(leiden_id)) & (ig == "Non-immune") & (l2 == unassigned)
            l2[mask] = l2_type
            source[mask] = f"leiden_{leiden_id}_direct"
            n_leiden_direct += mask.sum()
            log.info(f"  [leiden {leiden_id}] direct → {l2_type}  (n={mask.sum():,})")
    else:
        log.info("  [leiden_direct] disabled by config")

    # Step C: remaining Non-immune → per-cell score argmax (non-immune types)
    mask_remaining = (ig == "Non-immune") & (l2 == unassigned)
    score_cols_ni = [f"score_{ct}" for ct in nonimmune_l2 if f"score_{ct}" in obs.columns]
    if score_cols_ni and mask_remaining.sum() > 0:
        best = obs.loc[mask_remaining, score_cols_ni].idxmax(axis=1).str.replace("score_", "", regex=False)
        best_score = obs.loc[mask_remaining, score_cols_ni].max(axis=1)
        l2[mask_remaining] = best.where(best_score > score_threshold, other=unassigned)
        source[mask_remaining] = "score_nonimmune"
        log.info(f"  [score_nonimmune] {mask_remaining.sum():,} células → {l2[mask_remaining].value_counts().to_dict()}")

    # ── Map to L1 / L3 ───────────────────────────────────────────────────
    log.info("\n[4/6] Mapeando L1 y L3...")
    l1 = l2.map(L2_TO_L1).fillna(unassigned)
    l3 = l2.map(L2_TO_L3).fillna("Steady_state")

    # ── Stats ─────────────────────────────────────────────────────────────
    log.info("\n[5/6] Distribución v3:")
    log.info(f"\n  cell_type_L2_v3:")
    for ct, cnt in l2.value_counts().items():
        log.info(f"    {ct:<30} {cnt:>8,}  ({100*cnt/n:.1f}%)")
    log.info(f"\n  Fuentes de anotación:")
    for src, cnt in source.value_counts().items():
        log.info(f"    {src:<30} {cnt:>8,}  ({100*cnt/n:.1f}%)")

    # ── Write to sdata.zarr ───────────────────────────────────────────────
    log.info("\n[6/6] Escribiendo en sdata.zarr...")
    new_cols = {
        "cell_type_L1_v3": l1.astype("category"),
        "cell_type_L2_v3": l2.astype("category"),
        "cell_type_L3_v3": l3.astype("category"),
        "annotation_source_v3": source.astype("category"),
    }
    for col, series in new_cols.items():
        adata.obs[col] = series.values

    # Overwrite table in zarr
    sdata.delete_element_from_disk(table_key)
    sdata.write_element(table_key)
    log.info(f"  ✓ sdata.zarr actualizado con columnas v3")

    # Save summary CSV
    summary_df = pd.DataFrame({
        "cell_type_L2_v3": l2.values,
        "cell_type_L1_v3": l1.values,
        "annotation_source_v3": source.values,
    }, index=obs.index)
    summary_df.to_csv(args.outdir / "annotation_v3_per_cell.csv")

    # Distribution table
    dist = pd.DataFrame({
        "L2": l2.value_counts().index,
        "L1": [L2_TO_L1.get(x, unassigned) for x in l2.value_counts().index],
        "n_cells": l2.value_counts().values,
        "pct": (l2.value_counts().values / n * 100).round(1),
    })
    dist.to_csv(args.outdir / "distribution_v3.csv", index=False)

    # Source breakdown
    source.value_counts().to_frame("n_cells").to_csv(args.outdir / "annotation_sources_v3.csv")

    elapsed = time.time() - t0
    summary = {
        "version": "v3",
        "strategy": "immune_granular_anchor + leiden_direct + score_nonimmune",
        "markers_config": str(args.markers),
        "n_cells": int(n),
        "n_immune_annotated": int(n_immune),
        "n_leiden_direct": int(n_leiden_direct),
        "n_score_nonimmune": int(mask_remaining.sum()),
        "n_unassigned": int((l2 == unassigned).sum()),
        "n_l2_types": int(l2.nunique()),
        "execution_seconds": round(elapsed, 1),
        "reference_models": [
            "cell_type_immune_granular (subclustering Leiden inmune, computacional)",
            "CellTypist Human_Lung_Atlas.pkl (referencia HCA Lung)",
            "CellTypist Human_PF_Lung.pkl",
            "CellTypist Human_IPF_Lung.pkl",
            "CellTypist Immune_All_High.pkl",
            "CellTypist Immune_All_Low.pkl",
            "Marker scoring sc.tl.score_genes (panel INP-Pulmón 2025, 380+29 genes)",
        ],
        "annotation_sources": source.value_counts().to_dict(),
        "l2_distribution": l2.value_counts().to_dict(),
    }
    with open(args.outdir / "annotation_v3_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    log.info(f"\n✓ v3 completo en {elapsed:.1f}s")
    log.info(f"  {l2.nunique()} tipos L2 | {(l2==unassigned).sum():,} Unassigned ({100*(l2==unassigned).sum()/n:.1f}%)")
    log.info(f"  immune_granular: {n_immune:,} | leiden_direct: {n_leiden_direct:,} | score: {mask_remaining.sum():,}")


if __name__ == "__main__":
    main()
