#!/usr/bin/env python3
"""
F6b — Validación interna de anotación v2:
  1. cell_type_L2 vs leiden (15 clusters) — composición por cluster
  2. cell_type_L2 vs cell_type_immune_granular (11 tipos inmunes, subclustering real)
  3. Adjusted Rand Index (ARI) como métrica cuantitativa
  4. Recomendación de reconciliación
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
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import spatialdata as sd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
from utils.paths import sdata_path, results_root
log = get_logger(__name__)

SDATA_PATH = sdata_path()
OUTDIR = results_root() / "02_biology/F6_celltypist_validation"
OUTDIR.mkdir(parents=True, exist_ok=True)


def plot_heatmap_normalized(df, title, outpath, figsize=None):
    if figsize is None:
        figsize = (max(10, len(df.columns)*0.5)+2, max(8, len(df.index)*0.45)+1)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df, ax=ax, cmap="Blues", vmin=0, vmax=1,
                xticklabels=True, yticklabels=True, linewidths=0.3,
                cbar_kws={"label": "Proporción (norm. por fila)"})
    ax.set_title(title, fontsize=10)
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"  ✓ {outpath.name}")


def main():
    t0 = time.time()
    log.info("=" * 70)
    log.info("F6b — Validación interna de anotación L2 v2")
    log.info("=" * 70)

    sdata = sd.read_zarr(str(SDATA_PATH))
    table_key = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_key]
    obs = adata.obs.copy()
    log.info(f"  {len(obs):,} células cargadas")
    log.info(f"  Columnas: {list(obs.columns)}")

    results = {}

    # ─────────────────────────────────────────────────────────────
    # 1. L2 vs Leiden (15 clusters globales)
    # ─────────────────────────────────────────────────────────────
    log.info("\n[1/3] cell_type_L2 vs leiden (composición por cluster)...")
    conf1 = pd.crosstab(obs["leiden"].astype(str), obs["cell_type_L2"], normalize="index")
    conf1.to_csv(OUTDIR / "internal_L2_vs_leiden.csv")
    plot_heatmap_normalized(
        conf1,
        "Composición cell_type_L2 por cluster Leiden\n(prop. normalizada por cluster)",
        OUTDIR / "heatmap_L2_vs_leiden.png",
    )

    # Dominant L2 per Leiden cluster
    dom1 = pd.DataFrame({
        "leiden": conf1.index,
        "dominant_L2": conf1.idxmax(axis=1).values,
        "dominant_pct": conf1.max(axis=1).values,
        "n_cells": obs.groupby("leiden")["cell_type_L2"].count().reindex(conf1.index).values,
    })
    dom1.to_csv(OUTDIR / "dominant_L2_per_leiden.csv", index=False)
    log.info("\n  Tipo L2 dominante por cluster Leiden:")
    for _, r in dom1.iterrows():
        flag = "⚠️ MIXED" if r["dominant_pct"] < 0.4 else ("✓" if r["dominant_pct"] > 0.7 else "~")
        log.info(f"    Leiden {r['leiden']:>3}: {r['dominant_L2']:<25} {r['dominant_pct']*100:.1f}%  {flag}  (n={r['n_cells']:,})")

    ari_leiden = adjusted_rand_score(obs["leiden"].astype(str), obs["cell_type_L2"].astype(str))
    nmi_leiden = normalized_mutual_info_score(obs["leiden"].astype(str), obs["cell_type_L2"].astype(str))
    results["L2_vs_Leiden"] = {"ARI": round(ari_leiden, 4), "NMI": round(nmi_leiden, 4)}
    log.info(f"\n  ARI(L2, Leiden) = {ari_leiden:.4f}   NMI = {nmi_leiden:.4f}")

    # ─────────────────────────────────────────────────────────────
    # 2. L2 vs cell_type_immune_granular (11 tipos)
    # ─────────────────────────────────────────────────────────────
    log.info("\n[2/3] cell_type_L2 vs cell_type_immune_granular...")
    obs_immune = obs[obs["cell_type_immune_granular"].notna()].copy()
    obs_immune = obs_immune[obs_immune["cell_type_immune_granular"] != "Non_immune"]
    log.info(f"  Células inmunes anotadas: {len(obs_immune):,}")

    if len(obs_immune) > 0:
        conf2 = pd.crosstab(obs_immune["cell_type_immune_granular"], obs_immune["cell_type_L2"], normalize="index")
        conf2.to_csv(OUTDIR / "internal_L2_vs_immune_granular.csv")
        plot_heatmap_normalized(
            conf2,
            "cell_type_immune_granular (subclustering) vs cell_type_L2 v2\n(prop. normalizada por tipo inmune)",
            OUTDIR / "heatmap_L2_vs_immune_granular.png",
        )

        dom2 = pd.DataFrame({
            "immune_granular": conf2.index,
            "dominant_L2": conf2.idxmax(axis=1).values,
            "dominant_pct": conf2.max(axis=1).values,
            "n_cells": obs_immune.groupby("cell_type_immune_granular")["cell_type_L2"].count().reindex(conf2.index).values,
        })
        dom2.to_csv(OUTDIR / "dominant_L2_per_immune_granular.csv", index=False)
        log.info("\n  Tipo L2 dominante por tipo inmune granular (subclustering):")
        for _, r in dom2.iterrows():
            flag = "⚠️ DISCORDANTE" if r["dominant_pct"] < 0.4 else ("✓ CONCORDANTE" if r["dominant_pct"] > 0.6 else "~ PARCIAL")
            log.info(f"    {r['immune_granular']:<25} → {r['dominant_L2']:<25} {r['dominant_pct']*100:.1f}%  {flag}  (n={r['n_cells']:,})")

        ari_immune = adjusted_rand_score(obs_immune["cell_type_immune_granular"].astype(str), obs_immune["cell_type_L2"].astype(str))
        nmi_immune = normalized_mutual_info_score(obs_immune["cell_type_immune_granular"].astype(str), obs_immune["cell_type_L2"].astype(str))
        results["L2_vs_ImmuneGranular"] = {"ARI": round(ari_immune, 4), "NMI": round(nmi_immune, 4)}
        log.info(f"\n  ARI(L2, immune_granular) = {ari_immune:.4f}   NMI = {nmi_immune:.4f}")
    else:
        log.warning("  No hay células inmunes con anotación granular")

    # ─────────────────────────────────────────────────────────────
    # 3. L2 vs leiden_immune (sub-clusters inmunes)
    # ─────────────────────────────────────────────────────────────
    log.info("\n[3/3] cell_type_L2 vs leiden_immune...")
    obs_li = obs[obs["leiden_immune"].notna()].copy()
    log.info(f"  Células con leiden_immune: {len(obs_li):,}")

    if len(obs_li) > 0:
        conf3 = pd.crosstab(obs_li["leiden_immune"].astype(str), obs_li["cell_type_L2"], normalize="index")
        conf3.to_csv(OUTDIR / "internal_L2_vs_leiden_immune.csv")
        plot_heatmap_normalized(
            conf3,
            "Composición L2 v2 por leiden_immune\n(sub-clusters inmunes, prop. norm. por sub-cluster)",
            OUTDIR / "heatmap_L2_vs_leiden_immune.png",
            figsize=(16, 12),
        )
        ari_li = adjusted_rand_score(obs_li["leiden_immune"].astype(str), obs_li["cell_type_L2"].astype(str))
        nmi_li = normalized_mutual_info_score(obs_li["leiden_immune"].astype(str), obs_li["cell_type_L2"].astype(str))
        results["L2_vs_LeidenImmune"] = {"ARI": round(ari_li, 4), "NMI": round(nmi_li, 4)}
        log.info(f"  ARI(L2, leiden_immune) = {ari_li:.4f}   NMI = {nmi_li:.4f}")

    # Summary
    with open(OUTDIR / "internal_validation_summary.json", "w") as f:
        json.dump(results, f, indent=2)

    elapsed = time.time() - t0
    log.info(f"\n{'='*70}")
    log.info(f"✓ F6b completo en {elapsed:.1f}s")
    log.info(f"\nResumen ARI/NMI:")
    for k, v in results.items():
        ari = v['ARI']
        interp = "EXCELENTE" if ari > 0.6 else ("BUENA" if ari > 0.4 else ("MODERADA" if ari > 0.2 else "BAJA"))
        log.info(f"  {k}: ARI={ari:.4f} NMI={v['NMI']:.4f} → {interp}")


if __name__ == "__main__":
    main()
