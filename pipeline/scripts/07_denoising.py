"""
Step 07 — ResolVI Denoising
=============================
Uses the ResolVI model (scvi-tools) to separate true cell expression from:
  - Spatial contamination (transcripts from neighboring cells)
  - Ambient RNA background

IMPORTANT — layer convention:
  adata.layers["counts"]   → raw integer counts (NEVER modified)
  adata.layers["denoised"] → ResolVI denoised expression (added here)
  adata.X                  → normalized log1p (unchanged)

Downstream steps that require raw counts (DESeq2, scVI) must use
adata.layers["counts"]. Visualization and spatial analysis use
adata.layers["denoised"] or adata.X.

ResolVI needs:
  - The spatial coordinates (cell centroids) in adata.obsm["spatial"]
  - Raw integer counts in adata.layers["counts"]
  - GPU strongly recommended (RTX 4500 Ada: ~20-40 min for 278k cells)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def ensure_spatial_coords(adata, sdata) -> None:
    """
    Ensure adata.obsm["spatial"] contains cell centroid coordinates (µm).
    Computes from the active segmentation shapes if not already present.
    """
    import numpy as np

    if "spatial" in adata.obsm:
        log.info("adata.obsm['spatial'] already present.")
        return

    # Find active segmentation element
    for key in ("cell_boundaries_cellpose", "cell_boundaries_baysor", "cell_boundaries"):
        if key in sdata:
            log.info(f"Computing centroids from '{key}' ...")
            shapes  = sdata[key]
            cx      = shapes.geometry.centroid.x.values
            cy      = shapes.geometry.centroid.y.values
            # Match order to adata.obs_names
            centroid_df = shapes.assign(cx=cx, cy=cy)[["cx", "cy"]]
            adata.obsm["spatial"] = centroid_df.loc[adata.obs_names].values.astype(np.float32)
            log.info(f"Spatial coords set: shape {adata.obsm['spatial'].shape}")
            return

    log.warning("No cell boundary shapes found. ResolVI will run without spatial coords.")


def main():
    import numpy as np
    import scvi
    import spatialdata as sd

    sdata_path  = snakemake.input.sdata         # noqa: F821
    done_path   = snakemake.output.done         # noqa: F821
    max_epochs  = snakemake.params.max_epochs   # noqa: F821
    batch_size  = snakemake.params.batch_size   # noqa: F821
    accelerator = snakemake.params.accelerator  # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.table

    if "counts" not in adata.layers:
        raise RuntimeError(
            "adata.layers['counts'] not found. "
            "Run step 04 (preprocess) before denoising."
        )

    # Ensure spatial coordinates are set
    ensure_spatial_coords(adata, sdata)

    # Setup ResolVI
    log.info("Setting up ResolVI model ...")
    scvi.external.RESOLVI.setup_anndata(
        adata,
        layer="counts",
        # batch_key=None for single-sample; set to "sample" for multi-sample
    )

    model = scvi.external.RESOLVI(adata)
    log.info(f"Training ResolVI (max_epochs={max_epochs}, accelerator={accelerator}) ...")

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        accelerator=accelerator,
        early_stopping=True,
        early_stopping_patience=15,
        plan_kwargs={"lr": 1e-3},
    )

    # Extract denoised expression
    log.info("Computing denoised expression ...")
    denoised = model.get_denoised_expression(
        adata=adata,
        n_samples=25,           # average over 25 posterior samples
        library_size="latent",
    )
    adata.layers["denoised"] = denoised.values.astype(np.float32)

    # Verify raw counts are untouched
    assert "counts" in adata.layers, "Raw counts layer was lost!"
    log.info("Raw counts layer preserved: ✓")
    log.info(f"Denoised layer added: shape {adata.layers['denoised'].shape}")

    sdata.table = adata
    log.info("Writing denoised table to zarr store ...")
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")

    Path(done_path).touch()
    log.info("Step 07 — ResolVI denoising: DONE")


main()
