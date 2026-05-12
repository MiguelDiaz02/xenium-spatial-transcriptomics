#!/usr/bin/env python3
"""
F0b — DGEA con anotación granular (L2, 26 tipos)

Wilcoxon 1vRest sobre cell_type_L2 + integración Moran's I.
Sustituye/complementa Phase 1A (que corría sobre los 7 tipos L1 antiguos).

Salida:
  - dgea_L2_<celltype>.csv (un CSV por tipo)
  - dgea_L2_top50_all.csv (resumen)
  - dgea_L2_dotplot_top5.png
  - dgea_L2_summary.json
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import spatialdata as sd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
log = get_logger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--groupby", default="cell_type_L2",
                   help="Columna de agrupamiento (cell_type_L1, L2, o L3)")
    p.add_argument("--min_cells", type=int, default=100,
                   help="Mínimo de células por grupo para incluir")
    return p.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 70)
    log.info(f"F0b — DGEA Granular ({args.groupby})")
    log.info("=" * 70)

    sdata = sd.read_zarr(str(args.sdata))
    table = list(sdata.tables.keys())[0]
    adata = sdata.tables[table]
    log.info(f"  {adata.n_obs:,} células × {adata.n_vars} genes")

    # Filtrar grupos con muy pocas células
    counts = adata.obs[args.groupby].value_counts()
    keep = counts[counts >= args.min_cells].index.tolist()
    if "Unassigned" in keep:
        keep.remove("Unassigned")
    log.info(f"  Grupos con ≥{args.min_cells} células: {len(keep)}/{len(counts)}")
    log.info(f"  → {keep}")

    mask = adata.obs[args.groupby].isin(keep)
    adata_sub = adata[mask].copy()
    adata_sub.obs[args.groupby] = adata_sub.obs[args.groupby].cat.remove_unused_categories()
    log.info(f"  Subseteado: {adata_sub.n_obs:,} células")

    # Wilcoxon 1vRest
    log.info(f"\n[1/3] Wilcoxon 1-vs-rest...")
    sc.tl.rank_genes_groups(adata_sub, groupby=args.groupby,
                             method="wilcoxon", n_genes=adata_sub.n_vars,
                             key_added="wilcoxon_L2")

    # Extraer y guardar por grupo
    log.info(f"\n[2/3] Extrayendo y guardando resultados...")
    all_results = []
    for ct in adata_sub.obs[args.groupby].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_sub, group=ct, key="wilcoxon_L2")
        df["cell_type"] = ct
        df = df.sort_values("scores", ascending=False).reset_index(drop=True)
        df.to_csv(args.outdir / f"dgea_L2_{ct}.csv", index=False)
        all_results.append(df.head(50))
    top50 = pd.concat(all_results, ignore_index=True)
    top50.to_csv(args.outdir / "dgea_L2_top50_all.csv", index=False)
    log.info(f"  ✓ {len(adata_sub.obs[args.groupby].cat.categories)} CSVs por tipo + top50_all.csv")

    # Dotplot top 5
    log.info(f"\n[3/3] Dotplot top 5 marcadores por tipo...")
    top5_genes = []
    for ct in adata_sub.obs[args.groupby].cat.categories:
        df = sc.get.rank_genes_groups_df(adata_sub, group=ct, key="wilcoxon_L2")
        for g in df["names"].head(5).tolist():
            if g not in top5_genes:
                top5_genes.append(g)

    fig = sc.pl.dotplot(adata_sub, var_names=top5_genes[:60],
                         groupby=args.groupby, dendrogram=False,
                         standard_scale="var", show=False, return_fig=True,
                         figsize=(min(20, len(top5_genes)*0.3), max(8, len(keep)*0.4)))
    plt.savefig(args.outdir / "dgea_L2_dotplot_top5.png", dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"  ✓ Dotplot guardado")

    # Resumen JSON
    summary = {
        "groupby": args.groupby,
        "n_cells_input": int(adata.n_obs),
        "n_cells_used": int(adata_sub.n_obs),
        "n_groups_kept": len(keep),
        "groups_kept": keep,
        "min_cells_threshold": args.min_cells,
        "n_dgea_files": len(adata_sub.obs[args.groupby].cat.categories),
    }
    with open(args.outdir / "dgea_L2_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info(f"\n✓ DGEA granular completo → {args.outdir}")


if __name__ == "__main__":
    main()
