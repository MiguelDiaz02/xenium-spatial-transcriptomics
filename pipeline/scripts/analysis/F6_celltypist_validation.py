#!/usr/bin/env python3
"""
F6 — CellTypist multi-model annotation + concordance vs cell_type_L2 (v2)

Modelos: Human_Lung_Atlas, Human_PF_Lung, Human_IPF_Lung, Immune_All_High, Immune_All_Low
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
import matplotlib.colors as mcolors
import seaborn as sns

import spatialdata as sd
import scanpy as sc
import celltypist
from celltypist import models

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
from utils.paths import sdata_path, results_root
log = get_logger(__name__)

SDATA_PATH = sdata_path()
OUTDIR = results_root() / "02_biology/F6_celltypist_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)

MODELS = [
    "Human_Lung_Atlas.pkl",
    "Human_PF_Lung.pkl",
    "Human_IPF_Lung.pkl",
    "Immune_All_High.pkl",
    "Immune_All_Low.pkl",
]

def main():
    t0 = time.time()
    log.info("=" * 70)
    log.info("F6 — CellTypist multi-model validation")
    log.info("=" * 70)

    # Load data
    log.info(f"\n[1/4] Cargando sdata.zarr...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_key = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_key].copy()
    log.info(f"  {adata.n_obs:,} células × {adata.n_vars} genes")
    log.info(f"  cell_type_L2 disponible: {sorted(adata.obs['cell_type_L2'].unique().tolist())}")

    # X is already log1p-normalized (normalize_total 1e4 + log1p applied in pipeline)
    adata_ct = adata.copy()
    log.info(f"  X ya es log1p-normalizado (max={adata_ct.X.max():.2f}) — usando directamente")

    # Download/check models
    log.info(f"\n[2/4] Verificando modelos CellTypist (cached localmente)...")
    available = celltypist.models.get_all_models()
    for m in MODELS:
        if m in available:
            log.info(f"  ✓ {m} (cached)")
        else:
            log.warning(f"  ✗ {m} no encontrado en caché")

    # Run CellTypist for each model
    log.info(f"\n[3/4] Corriendo predicciones ({len(MODELS)} modelos)...")
    results_dict = {}
    for model_name in MODELS:
        log.info(f"\n  → {model_name}...")
        try:
            t1 = time.time()
            model = models.Model.load(model=model_name)
            predictions = celltypist.annotate(
                adata_ct,
                model=model,
                majority_voting=True,
                over_clustering="cell_type_L2",
            )
            pred_df = predictions.predicted_labels
            col_pred = f"ct_{model_name}"
            col_prob = f"ct_{model_name}_conf"
            adata.obs[col_pred] = pred_df["majority_voting"].values
            adata.obs[col_prob] = predictions.probability_matrix.max(axis=1).values
            results_dict[model_name] = pred_df["majority_voting"].values
            log.info(f"    ✓ {model_name} en {time.time()-t1:.1f}s — {pred_df['majority_voting'].nunique()} tipos predichos")
        except Exception as e:
            log.error(f"    ✗ {model_name} falló: {e}")
            import traceback; traceback.print_exc()

    # Concordance analysis
    log.info(f"\n[4/4] Concordancia CellTypist vs cell_type_L2...")
    summary = {}
    for model_name, preds in results_dict.items():
        ct_col = f"ct_{model_name}"
        if ct_col not in adata.obs.columns:
            continue

        # Confusion matrix (L2 rows, CellTypist cols) — normalized by row
        l2_types = sorted(adata.obs["cell_type_L2"].unique())
        ct_types = sorted(adata.obs[ct_col].unique())

        conf = pd.crosstab(
            adata.obs["cell_type_L2"],
            adata.obs[ct_col],
            normalize="index",
        )
        conf.to_csv(OUTDIR / f"confusion_{model_name}.csv")

        # Top CellTypist match per L2 type
        top_match = conf.idxmax(axis=1)
        top_pct = conf.max(axis=1)
        match_df = pd.DataFrame({
            "L2_type": top_match.index,
            "celltypist_match": top_match.values,
            "match_pct": top_pct.values,
        })
        match_df.to_csv(OUTDIR / f"top_match_{model_name}.csv", index=False)

        # Overall concordance: % cells where L2 maps >50% to a single CT label
        high_conf = (top_pct > 0.5).sum()
        concordance_pct = 100 * high_conf / len(l2_types)
        summary[model_name] = {
            "n_ct_types": len(ct_types),
            "n_l2_high_concordance": int(high_conf),
            "n_l2_types": len(l2_types),
            "pct_l2_with_dominant_ct": round(concordance_pct, 1),
        }
        log.info(f"  {model_name}: {high_conf}/{len(l2_types)} tipos L2 con >50% concordancia ({concordance_pct:.1f}%)")

        # Heatmap
        fig, ax = plt.subplots(figsize=(min(30, len(ct_types)*0.35)+2, min(20, len(l2_types)*0.45)+1))
        sns.heatmap(conf, ax=ax, cmap="Blues", vmin=0, vmax=1,
                    xticklabels=True, yticklabels=True, linewidths=0.3,
                    cbar_kws={"label": "Proporción"})
        ax.set_xlabel(f"CellTypist: {model_name}", fontsize=9)
        ax.set_ylabel("cell_type_L2 (v2)", fontsize=9)
        ax.set_title(f"Concordancia L2 vs {model_name}\n({concordance_pct:.1f}% tipos L2 con >50% acuerdo)", fontsize=10)
        plt.xticks(rotation=45, ha="right", fontsize=6)
        plt.yticks(fontsize=7)
        plt.tight_layout()
        plt.savefig(OUTDIR / f"heatmap_{model_name}.png", dpi=150, bbox_inches="tight")
        plt.close()
        log.info(f"    ✓ heatmap_{model_name}.png")

    # Save summary
    with open(OUTDIR / "concordance_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Save all predictions to CSV for inspection
    pred_cols = ["cell_type_L1", "cell_type_L2", "cell_type_L3", "cell_type_confidence"] + \
                [f"ct_{m}" for m in MODELS if f"ct_{m}" in adata.obs.columns] + \
                [f"ct_{m}_conf" for m in MODELS if f"ct_{m}_conf" in adata.obs.columns]
    pred_cols = [c for c in pred_cols if c in adata.obs.columns]
    adata.obs[pred_cols].to_csv(OUTDIR / "all_predictions.csv")

    elapsed = time.time() - t0
    log.info(f"\n{'='*70}")
    log.info(f"✓ F6 completo en {elapsed:.1f}s → {OUTDIR}")
    log.info(f"\nResumen de concordancia:")
    for m, s in summary.items():
        log.info(f"  {m}: {s['pct_l2_with_dominant_ct']}% L2 con >50% acuerdo ({s['n_ct_types']} tipos CT predichos)")

if __name__ == "__main__":
    main()
