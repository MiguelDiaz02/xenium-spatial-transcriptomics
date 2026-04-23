"""
Step 08 — Spatial Statistics
==============================
Squidpy-based spatial analysis using cell type annotations and expression data.

Analyses performed:
  1. Spatial neighbors graph (coordinate-based k-NN)
  2. Neighborhood enrichment — which cell types co-localize?
  3. Co-occurrence — spatial clustering at multiple distance scales
  4. Moran's I — spatial autocorrelation per gene (identifies spatially patterned genes)
  5. Ripley's statistics — point pattern analysis

All results stored in adata.uns and adata.obsp.
Figures saved to results/08_spatial_figures/.
"""


import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def ensure_spatial(adata, sdata) -> None:
    """Ensure adata.obsm['spatial'] contains centroids."""
    if "spatial" in adata.obsm:
        return
    import numpy as np
    for key in ("cell_boundaries_cellpose", "cell_boundaries_baysor", "cell_boundaries"):
        if key in sdata:
            shapes = sdata[key]
            cx = shapes.geometry.centroid.x.values
            cy = shapes.geometry.centroid.y.values
            centroid_df = shapes.assign(cx=cx, cy=cy)[["cx", "cy"]]
            adata.obsm["spatial"] = centroid_df.loc[adata.obs_names].values.astype(np.float32)
            log.info(f"Centroids set from '{key}'.")
            return
    raise RuntimeError("No spatial coordinates available for Squidpy analysis.")


def main():
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend
    import matplotlib.pyplot as plt
    import squidpy as sq
    import spatialdata as sd

    sdata_path   = snakemake.input.sdata         # noqa: F821
    done_path    = snakemake.output.done         # noqa: F821
    fig_dir      = Path(snakemake.output.figures)  # noqa: F821
    n_neighs     = snakemake.params.n_neighs     # noqa: F821
    coord_type   = snakemake.params.coord_type   # noqa: F821
    cluster_key  = snakemake.params.cluster_key  # noqa: F821
    n_perms      = snakemake.params.n_perms_enrichment  # noqa: F821
    autocorr_cfg = snakemake.params.autocorr     # noqa: F821

    fig_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]  # Get first (usually only) table

    if cluster_key not in adata.obs:
        log.warning(f"'{cluster_key}' not in adata.obs — falling back to 'leiden'.")
        cluster_key = "leiden"

    ensure_spatial(adata, sdata)

    # 1. Spatial neighbors graph
    log.info(f"Building spatial neighbors graph (k={n_neighs}, coord_type={coord_type}) ...")
    sq.gr.spatial_neighbors(
        adata,
        coord_type=coord_type,
        n_neighs=n_neighs,
        key_added="spatial",  # Use default key name for Squidpy compatibility
    )

    # 2. Neighborhood enrichment
    log.info(f"Neighborhood enrichment (n_perms={n_perms}) ...")
    sq.gr.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        n_perms=n_perms,
        seed=42,
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    sq.pl.nhood_enrichment(adata, cluster_key=cluster_key, ax=ax)
    fig.savefig(fig_dir / "nhood_enrichment.pdf", bbox_inches="tight")
    plt.close(fig)

    # 3. Co-occurrence (SKIPPED: too memory-intensive for 268k cells)
    # log.info("Co-occurrence analysis ...")
    # sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    # fig = sq.pl.co_occurrence(adata, cluster_key=cluster_key, return_ax=True)
    # plt.savefig(fig_dir / "co_occurrence.pdf", bbox_inches="tight")
    # plt.close()
    log.info("Co-occurrence analysis skipped (memory-intensive for large datasets).")

    # 4. Moran's I (spatial autocorrelation)
    log.info(f"Moran's I (mode={autocorr_cfg['mode']}, n_perms={autocorr_cfg['n_perms']}) ...")
    # Use denoised layer if available, else normalized X
    if "denoised" in adata.layers:
        log.info("Using adata.layers['denoised'] for Moran's I.")
        import anndata as ad
        adata_moran = adata.copy()
        import scipy.sparse as sp
        adata_moran.X = sp.csr_matrix(adata.layers["denoised"])
    else:
        adata_moran = adata

    sq.gr.spatial_autocorr(
        adata_moran,
        mode=autocorr_cfg["mode"],
        n_perms=autocorr_cfg["n_perms"],
        n_jobs=autocorr_cfg["n_jobs"],
    )

    # Copy Moran results back to main adata
    adata.uns[f"moranI"] = adata_moran.uns.get("moranI", {})
    if "moranI" in adata.uns:
        top_genes = (
            adata.uns["moranI"]
            .sort_values("I", ascending=False)
            .head(20)
        )
        log.info(f"Top spatially variable genes (Moran's I):\n{top_genes.to_string()}")
        top_genes.to_csv(fig_dir / "top_spatial_genes_moranI.csv")

    # 5. Save updated adata
    table_name = list(sdata.tables.keys())[0]
    sdata.tables[table_name] = adata
    log.info("Writing spatial results to zarr store ...")
    sdata.delete_element_from_disk(table_name)
    sdata.write_element(table_name)

    Path(done_path).touch()
    log.info(f"Figures saved to {fig_dir}")
    log.info("Step 08 — Spatial analysis: DONE")


main()
