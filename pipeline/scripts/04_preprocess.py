"""
Step 04 — Preprocessing
========================
1. Preserve raw counts in adata.layers["counts"] (required by ResolVI and DESeq2).
2. Normalize total counts per cell (target_sum from config).
3. Log1p transform.
4. Optionally select highly variable genes (for large panels; disabled for <500 genes).

After this step:
  adata.X              → normalized + log1p (used for PCA/clustering)
  adata.layers["counts"] → raw integer counts (used by scVI, ResolVI, DESeq2)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def main():
    import scipy.sparse as sp
    import scanpy as sc
    import spatialdata as sd

    sdata_path = snakemake.input.sdata       # noqa: F821
    done_path  = snakemake.output.done       # noqa: F821
    target_sum = snakemake.params.target_sum  # noqa: F821
    do_log1p   = snakemake.params.log1p      # noqa: F821
    select_hvg = snakemake.params.select_hvg  # noqa: F821
    n_hvg      = snakemake.params.n_hvg      # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.table

    # 1. Preserve raw counts
    if "counts" not in adata.layers:
        log.info("Saving raw counts to adata.layers['counts'] ...")
        adata.layers["counts"] = adata.X.copy() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    else:
        log.info("adata.layers['counts'] already present.")

    # 2. Normalize
    log.info(f"Normalizing total counts (target_sum={target_sum}) ...")
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # 3. Log1p
    if do_log1p:
        log.info("Applying log1p transform ...")
        sc.pp.log1p(adata)

    # 4. Highly variable genes (only useful for large panels)
    n_genes = adata.n_vars
    if select_hvg and n_genes > 200:
        log.info(f"Selecting top {n_hvg} highly variable genes (from {n_genes}) ...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3",
                                    layer="counts")
        adata = adata[:, adata.var["highly_variable"]].copy()
        log.info(f"After HVG selection: {adata.n_vars} genes retained.")
    else:
        log.info(f"Skipping HVG selection (panel has {n_genes} genes; select_hvg={select_hvg}).")

    sdata.table = adata
    log.info("Writing preprocessed table to zarr store ...")
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")

    Path(done_path).touch()
    log.info("Step 04 — Preprocess: DONE")


main()
