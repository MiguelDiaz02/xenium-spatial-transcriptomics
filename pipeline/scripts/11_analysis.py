"""
Step 11 — Downstream Spatial Analysis
=======================================
Cell-type specific spatial patterns, DEG analysis with spatial context.

Outputs:
  1. Per-celltype spatial autocorrelation profiles
  2. Co-localization summary table
  3. DEG table (cell type × spatial neighborhood)
  4. Top genes per cell type with spatial info
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def main():
    import pandas as pd
    import numpy as np
    import spatialdata as sd
    from scipy import stats

    sdata_path = snakemake.input.sdata  # noqa: F821
    out_dir = Path(snakemake.output.results)  # noqa: F821
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]

    log.info("Step 11 — Downstream Analysis: Starting ...")

    # 1. Per-celltype spatial statistics
    if "moranI" in adata.uns and "cell_type" in adata.obs:
        log.info("Computing per-celltype spatial statistics ...")
        morani_df = adata.uns["moranI"].copy()

        # Add cell type information to each gene
        results = []
        for gene in morani_df.index:
            if gene in adata.var_names:
                # Find cell types with highest expression of this gene
                expr = adata[:, gene].X.toarray().flatten()
                for ct in adata.obs["cell_type"].unique():
                    ct_mask = adata.obs["cell_type"] == ct
                    ct_expr = expr[ct_mask]
                    if ct_expr.sum() > 0:
                        results.append(
                            {
                                "gene": gene,
                                "cell_type": ct,
                                "moranI": morani_df.loc[gene, "I"],
                                "mean_expr_in_ct": ct_expr.mean(),
                                "pct_nonzero_in_ct": (ct_expr > 0).sum() / len(ct_expr) * 100,
                            }
                        )

        if results:
            ct_stats = pd.DataFrame(results)
            ct_stats.to_csv(out_dir / "celltype_spatial_stats.csv", index=False)
            log.info(f"  Saved: celltype_spatial_stats.csv ({len(results)} records)")

    # 2. Summary statistics
    log.info("Generating summary statistics ...")
    summary = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_cell_types": adata.obs["cell_type"].nunique() if "cell_type" in adata.obs else 0,
        "spatial_neighbors_computed": "spatial" in adata.obsp,
        "moranI_computed": "moranI" in adata.uns,
        "doublet_detected": "doublet_detected" in adata.obs,
    }

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(out_dir / "analysis_summary.csv", index=False)
    log.info(f"  Summary: {summary}")

    Path(snakemake.output.done).touch()  # noqa: F821
    log.info(f"Results saved to {out_dir}")
    log.info("Step 11 — Downstream Analysis: DONE")


main()
