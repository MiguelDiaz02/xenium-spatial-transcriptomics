"""
Step 05 — Dimensionality Reduction & Clustering
================================================
PCA → k-NN neighbors → UMAP → Leiden clustering.
Uses the normalized+log1p expression (adata.X) from step 04.

Outputs in adata:
  adata.obsm["X_pca"]    → PCA embedding
  adata.obsm["X_umap"]   → UMAP embedding
  adata.obs["leiden"]    → cluster labels
"""


import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def main():
    import scanpy as sc
    import spatialdata as sd

    sdata_path  = snakemake.input.sdata         # noqa: F821
    done_path   = snakemake.output.done         # noqa: F821
    n_pcs       = snakemake.params.n_pcs        # noqa: F821
    n_neighbors = snakemake.params.n_neighbors  # noqa: F821
    min_dist    = snakemake.params.umap_min_dist  # noqa: F821
    resolution  = snakemake.params.resolution   # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]  # Get first (usually only) table

    # Clip n_pcs to available genes
    n_pcs = min(n_pcs, adata.n_vars - 1)
    log.info(f"PCA: {n_pcs} components on {adata.n_obs:,} cells × {adata.n_vars} genes ...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    log.info(f"Computing k-NN neighbors (k={n_neighbors}) ...")
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    log.info(f"UMAP (min_dist={min_dist}) ...")
    sc.tl.umap(adata, min_dist=min_dist)

    log.info(f"Leiden clustering (resolution={resolution}) ...")
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    log.info(f"Leiden: {n_clusters} clusters detected.")

    table_name = list(sdata.tables.keys())[0]
    sdata.tables[table_name] = adata
    log.info("Writing table with embeddings to zarr store ...")
    sdata.delete_element_from_disk(table_name)
    sdata.write_element(table_name)

    Path(done_path).touch()
    log.info("Step 05 — Reduction: DONE")


main()
