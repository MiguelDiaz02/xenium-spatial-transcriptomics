#!/usr/bin/env python3
"""
F0c — CCC LIANA+ rank_aggregate sobre cell_type_L2 (26 tipos granulares)

Reemplaza phase2b_ccc_liana.py (que corría sobre 7 tipos L1 antiguos).
Mismo método, mismos parámetros, anotación granular.
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
import matplotlib.pyplot as plt

import spatialdata as sd
import scanpy as sc
import liana as li
import liana.method as lm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
log = get_logger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--groupby", default="cell_type_L2")
    p.add_argument("--n_perms", type=int, default=1000)
    p.add_argument("--min_cells", type=int, default=50)
    p.add_argument("--expr_prop", type=float, default=0.1)
    p.add_argument("--max_pval", type=float, default=0.05)
    return p.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    log.info("=" * 70)
    log.info(f"F0c — CCC LIANA+ Granular ({args.groupby})")
    log.info("=" * 70)

    t0 = time.time()
    sdata = sd.read_zarr(str(args.sdata))
    table = list(sdata.tables.keys())[0]
    adata = sdata.tables[table]
    log.info(f"  {adata.n_obs:,} células × {adata.n_vars} genes")
    log.info(f"  {args.groupby}: {sorted(adata.obs[args.groupby].dropna().unique().tolist())}")

    # Excluir Unassigned
    if "Unassigned" in adata.obs[args.groupby].cat.categories:
        keep = adata.obs[args.groupby] != "Unassigned"
        adata = adata[keep].copy()
        adata.obs[args.groupby] = adata.obs[args.groupby].cat.remove_unused_categories()
        log.info(f"  Sin Unassigned: {adata.n_obs:,} células")

    # LIANA+ rank_aggregate
    log.info(f"\n[1/3] LIANA+ rank_aggregate (n_perms={args.n_perms})...")
    lm.rank_aggregate(
        adata, groupby=args.groupby, resource_name="consensus",
        expr_prop=args.expr_prop, min_cells=args.min_cells, n_perms=args.n_perms,
        use_raw=False, verbose=True, return_all_lrs=False, seed=1337,
    )
    liana_res = adata.uns["liana_res"].copy()
    log.info(f"  Total interacciones: {len(liana_res):,}")

    # Filtrar self-interacciones + significancia con BH FDR
    # AUDIT FIX (2026-05-06): the previous logic used rank-based filtering when
    # `cellphone_pvals` was absent, which is NOT a true FDR criterion. The
    # manuscript reported "FDR < 5%" — we now apply Benjamini-Hochberg explicitly.
    from statsmodels.stats.multitest import multipletests
    log.info(f"\n[2/3] Filtrando significancia (BH FDR < {args.max_pval}) ...")
    df = liana_res[liana_res["source"] != liana_res["target"]].copy()
    log.info(f"  Sin self-interacciones: {len(df):,}")

    if "cellphone_pvals" not in df.columns:
        raise RuntimeError(
            "LIANA+ output missing `cellphone_pvals` column — cannot compute FDR. "
            "Ensure rank_aggregate was run with CellPhoneDB enabled."
        )

    pvals = df["cellphone_pvals"].fillna(1.0).values
    _, padj, _, _ = multipletests(pvals, alpha=args.max_pval, method="fdr_bh")
    df["cellphone_pvals_BH_FDR"] = padj
    df_sig = df[padj < args.max_pval].copy()
    log.info(f"  Total tests: {len(df):,}")
    log.info(f"  Significant (BH FDR < {args.max_pval}): {len(df_sig):,}")

    # Guardar
    df.to_csv(args.outdir / "liana_L2_all_interactions.csv", index=False)
    df_sig.to_csv(args.outdir / "liana_L2_significant.csv", index=False)
    log.info(f"  ✓ Guardado: liana_L2_all_interactions.csv + liana_L2_significant.csv")

    # Top hits
    if "magnitude_rank" in df_sig.columns:
        top_hits = df_sig.nsmallest(50, "magnitude_rank")
    else:
        top_hits = df_sig.head(50)
    top_hits.to_csv(args.outdir / "liana_L2_top50.csv", index=False)

    # Resumen JSON
    log.info(f"\n[3/3] Resumen y figura...")
    cell_types = sorted(adata.obs[args.groupby].unique().tolist())
    pair_counts = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)
    for _, row in df_sig.iterrows():
        if row["source"] in pair_counts.index and row["target"] in pair_counts.columns:
            pair_counts.loc[row["source"], row["target"]] += 1
    pair_counts.to_csv(args.outdir / "liana_L2_celltype_pair_counts.csv")

    # Heatmap
    fig, ax = plt.subplots(figsize=(min(20, len(cell_types)*0.5)+2, min(18, len(cell_types)*0.5)+1))
    im = ax.imshow(pair_counts.values, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(cell_types)))
    ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=7)
    ax.set_xlabel("Target (L2)")
    ax.set_ylabel("Source (L2)")
    ax.set_title(f"LIANA+ CCC: # interacciones significativas por par L2 (p<{args.max_pval})")
    plt.colorbar(im, ax=ax, label="# L-R pares")
    plt.tight_layout()
    plt.savefig(args.outdir / "liana_L2_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()

    elapsed = time.time() - t0
    summary = {
        "groupby": args.groupby,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_celltypes": len(cell_types),
        "celltypes": cell_types,
        "n_perms": args.n_perms,
        "min_cells": args.min_cells,
        "expr_prop": args.expr_prop,
        "n_interactions_total": int(len(df)),
        "n_interactions_significant": int(len(df_sig)),
        "max_pval": args.max_pval,
        "execution_seconds": round(elapsed, 1),
    }
    with open(args.outdir / "liana_L2_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log.info(f"\n✓ CCC LIANA+ granular completo en {elapsed:.1f}s → {args.outdir}")


if __name__ == "__main__":
    main()
