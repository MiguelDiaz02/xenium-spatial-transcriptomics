#!/usr/bin/env python3
"""
F0 — Hierarchical Cell Type Re-annotation (L1 / L2 / L3)
=========================================================

Reads `pipeline/config/cell_markers_lung.yaml`, scores each signature via
`scanpy.tl.score_genes`, and assigns cells at three resolution levels:

  L1 — Broad (existing 7 types, retained for backwards compat)
  L2 — Granular (~17 expected: AT1, Ciliated, CD8_effector, CD8_exhausted, ...)
  L3 — Functional states (Proliferating, Cytotoxic_effector, Exhaustion, TLS_signature, ...)

Outputs:
  - Adds columns `cell_type_L2`, `cell_type_L3_states` to adata.obs
  - Writes diagnostic CSV: results/02_biology/reannotation/{level}_scores.csv
  - Writes summary report: results/02_biology/reannotation/reannotation_summary.md
  - Persists updated table back into sdata.zarr

Usage:
  python F0_reannotation.py \\
      --sdata human_lung_cancer/results/sdata.zarr \\
      --config pipeline/config/cell_markers_lung.yaml \\
      --outdir human_lung_cancer/results/02_biology/reannotation

Citation context:
  - Marker-based scoring: Tirosh et al. Science 352, 189-196 (2016)
  - scanpy.tl.score_genes: Wolf et al. Genome Biol 19, 15 (2018)
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import matplotlib
matplotlib.use("Agg")

import spatialdata as sd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument(
        "--write-back",
        action="store_true",
        help="Write updated table back into sdata.zarr (in-place)",
    )
    return p.parse_args()


def load_markers(config_path: Path, available_genes: set[str]) -> dict:
    """Load marker hierarchy and filter to genes present in panel."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    def _filter(markers: list[str]) -> list[str]:
        return [g for g in markers if g in available_genes]

    filtered = {"level_1": {}, "level_2": {}, "level_3": {}}
    for level in ["level_1", "level_2", "level_3"]:
        for ctype, info in cfg.get(level, {}).items():
            present = _filter(info.get("markers", []))
            absent = [g for g in info.get("markers", []) if g not in available_genes]
            neg_present = _filter(info.get("negative_markers", []) or [])
            if len(present) >= 1:  # require ≥1 marker present
                filtered[level][ctype] = {
                    "markers": present,
                    "negative_markers": neg_present,
                    "absent_markers": absent,
                    "parent": info.get("parent"),
                    "description": info.get("description", ""),
                }
            else:
                log.warning(f"Skipping {ctype}: 0 of {len(info.get('markers', []))} markers in panel")

    params = cfg.get("annotation_params", {})
    return filtered, params


def score_signatures(adata, sig_dict: dict, prefix: str, ctrl_size: int = 50) -> pd.DataFrame:
    """Run scanpy.tl.score_genes for each signature; return scores DataFrame."""
    score_cols = []
    for ctype, info in sig_dict.items():
        gene_list = info["markers"]
        if len(gene_list) == 0:
            continue
        score_name = f"{prefix}__{ctype}"
        try:
            sc.tl.score_genes(
                adata,
                gene_list=gene_list,
                ctrl_size=min(ctrl_size, max(20, len(gene_list) * 5)),
                score_name=score_name,
                use_raw=False,
            )
            score_cols.append(score_name)
        except Exception as e:
            log.warning(f"Failed scoring {ctype}: {e}")
    return adata.obs[score_cols]


def assign_argmax(scores_df: pd.DataFrame, prefix: str, z_thr: float = 1.0,
                  ambig_thr: float = 0.3) -> pd.Series:
    """Per-cell argmax assignment with z-score threshold + ambiguity flag."""
    types = [c.replace(f"{prefix}__", "") for c in scores_df.columns]
    arr = scores_df.values
    # z-score across cells per signature
    arr_z = (arr - arr.mean(axis=0, keepdims=True)) / (arr.std(axis=0, keepdims=True) + 1e-9)

    # top1 / top2 ranking
    top1_idx = arr_z.argmax(axis=1)
    sorted_arr = np.sort(arr_z, axis=1)
    top1_z = sorted_arr[:, -1]
    top2_z = sorted_arr[:, -2] if arr_z.shape[1] > 1 else np.full_like(top1_z, -np.inf)

    assignments = np.array([types[i] for i in top1_idx], dtype=object)
    assignments[top1_z < z_thr] = "Unassigned"
    margin = top1_z - top2_z
    ambiguous = (top1_z >= z_thr) & (margin < ambig_thr)
    second_idx = np.argsort(-arr_z, axis=1)[:, 1] if arr_z.shape[1] > 1 else top1_idx
    for i in np.where(ambiguous)[0]:
        a, b = sorted([types[top1_idx[i]], types[second_idx[i]]])
        assignments[i] = f"Mixed_{a}_{b}"

    return pd.Series(assignments, index=scores_df.index, name=f"cell_type_{prefix}")


def assign_multilabel(scores_df: pd.DataFrame, prefix: str, z_thr: float = 1.0) -> pd.Series:
    """Multi-label assignment for L3 functional states (a cell can have multiple)."""
    types = [c.replace(f"{prefix}__", "") for c in scores_df.columns]
    arr = scores_df.values
    arr_z = (arr - arr.mean(axis=0, keepdims=True)) / (arr.std(axis=0, keepdims=True) + 1e-9)
    above = arr_z >= z_thr  # boolean per (cell, signature)
    labels = []
    for i in range(arr_z.shape[0]):
        flagged = [types[j] for j in np.where(above[i])[0]]
        labels.append("|".join(flagged) if flagged else "None")
    return pd.Series(labels, index=scores_df.index, name=f"cell_type_{prefix}")


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log.info("=" * 70)
    log.info("F0 — Hierarchical Cell Type Re-annotation")
    log.info("=" * 70)

    log.info(f"Loading {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  AnnData: {adata.shape}")
    log.info(f"  Existing cell_type column: present={'cell_type' in adata.obs.columns}")

    log.info(f"Loading marker hierarchy from {args.config}")
    available = set(adata.var_names)
    markers, params = load_markers(args.config, available)
    log.info(f"  L1 types loaded: {len(markers['level_1'])}")
    log.info(f"  L2 types loaded: {len(markers['level_2'])}")
    log.info(f"  L3 states loaded: {len(markers['level_3'])}")

    z_thr = float(params.get("z_score_threshold", 1.0))
    ambig_thr = float(params.get("ambiguity_threshold", 0.3))
    ctrl_size = int(params.get("ctrl_size", 50))

    # ── L2 scoring + assignment ────────────────────────────────────────────
    log.info("Scoring L2 signatures...")
    l2_scores = score_signatures(adata, markers["level_2"], prefix="L2", ctrl_size=ctrl_size)
    l2_assign = assign_argmax(l2_scores, prefix="L2", z_thr=z_thr, ambig_thr=ambig_thr)
    adata.obs["cell_type_L2"] = pd.Categorical(l2_assign.values)
    log.info(f"  L2 distribution:\n{adata.obs['cell_type_L2'].value_counts().to_string()}")

    # ── L3 scoring + multilabel ────────────────────────────────────────────
    log.info("Scoring L3 functional states (multilabel)...")
    l3_scores = score_signatures(adata, markers["level_3"], prefix="L3", ctrl_size=ctrl_size)
    l3_assign = assign_multilabel(l3_scores, prefix="L3", z_thr=z_thr)
    adata.obs["cell_type_L3_states"] = pd.Categorical(l3_assign.values)
    log.info(f"  L3 distribution (top 10):\n{adata.obs['cell_type_L3_states'].value_counts().head(10).to_string()}")

    # ── Save scores CSVs ───────────────────────────────────────────────────
    l2_scores.to_csv(args.outdir / "L2_scores.csv")
    l3_scores.to_csv(args.outdir / "L3_scores.csv")
    log.info(f"  Saved: {args.outdir}/L2_scores.csv ({l2_scores.shape})")
    log.info(f"  Saved: {args.outdir}/L3_scores.csv ({l3_scores.shape})")

    # ── Save assignments table for downstream ──────────────────────────────
    assignments = adata.obs[["cell_type", "cell_type_L2", "cell_type_L3_states"]].copy()
    if "cell_type_immune_granular" in adata.obs.columns:
        assignments["cell_type_immune_granular"] = adata.obs["cell_type_immune_granular"]
    assignments.to_csv(args.outdir / "cell_assignments.csv")
    log.info(f"  Saved: {args.outdir}/cell_assignments.csv")

    # ── Cross-tab L1 × L2 (sanity check) ───────────────────────────────────
    if "cell_type" in adata.obs.columns:
        ct = pd.crosstab(adata.obs["cell_type"], adata.obs["cell_type_L2"])
        ct.to_csv(args.outdir / "L1_vs_L2_crosstab.csv")
        log.info(f"  Saved: {args.outdir}/L1_vs_L2_crosstab.csv (shape {ct.shape})")

    # ── Summary markdown ───────────────────────────────────────────────────
    summary = [
        "# F0 — Cell Type Re-annotation Summary",
        f"\n**Total cells:** {adata.n_obs:,}",
        f"**Total genes:** {adata.n_vars}",
        f"\n## L2 granular distribution",
        adata.obs["cell_type_L2"].value_counts().to_markdown(),
        f"\n## L3 functional states (top 10)",
        adata.obs["cell_type_L3_states"].value_counts().head(10).to_markdown(),
    ]
    (args.outdir / "reannotation_summary.md").write_text("\n".join(summary))
    log.info(f"  Saved: {args.outdir}/reannotation_summary.md")

    # ── Optional: persist into zarr ────────────────────────────────────────
    if args.write_back:
        log.info("Writing back to sdata.zarr...")
        # Replace the table in-place; spatialdata supports this via write_consolidated_metadata
        sdata.tables[table_name] = adata
        # Persist updated obs (write to zarr table)
        sdata.write_element(element_name=table_name, element_type="tables", overwrite=True) \
            if hasattr(sdata, "write_element") else None
        log.info("  Done.")
    else:
        log.info("Skip persisting to zarr (--write-back not set)")

    log.info(f"\n✓ F0 complete in {time.time() - t0:.1f}s")
    log.info(f"  Outputs at: {args.outdir}")


if __name__ == "__main__":
    main()
