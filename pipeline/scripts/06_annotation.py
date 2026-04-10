"""
Step 06 — Cell Type Annotation
================================
Annotates cells using CellTypist (default), manual marker scoring, or scANVI.

CellTypist method:
  - Downloads pre-trained model (cached locally after first run)
  - Applies majority voting within Leiden clusters for robustness
  - Stores result in adata.obs["cell_type"] and adata.obs["cell_type_fine"]

Manual method:
  - Scores each cell against marker gene sets from config
  - Assigns the highest-scoring label
  - Good for validation or when CellTypist model doesn't fit the tissue

Both methods always produce:
  adata.obs["cell_type"]       → coarse label (for spatial analysis)
  adata.obs["cell_type_fine"]  → fine-grained label (when available)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def annotate_celltypist(adata, model_name: str, majority_voting: bool):
    """Run CellTypist annotation."""
    import celltypist
    from celltypist import models as ct_models

    log.info(f"CellTypist: loading model '{model_name}' ...")
    # Downloads model if not cached
    model = ct_models.Model.load(model=model_name)

    log.info("Running CellTypist annotation ...")
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
        over_clustering="leiden",  # use existing Leiden clusters
    )
    adata = predictions.to_adata()

    # Standardize column names
    if "majority_voting" in adata.obs:
        adata.obs["cell_type"] = adata.obs["majority_voting"]
    elif "predicted_labels" in adata.obs:
        adata.obs["cell_type"] = adata.obs["predicted_labels"]
    if "predicted_labels" in adata.obs:
        adata.obs["cell_type_fine"] = adata.obs["predicted_labels"]

    n_types = adata.obs["cell_type"].nunique()
    log.info(f"CellTypist: {n_types} cell types annotated.")
    log.info(adata.obs["cell_type"].value_counts().to_string())
    return adata


def annotate_manual(adata, markers: dict[str, list[str]]):
    """Score cells against marker gene sets and assign labels by max score."""
    import scanpy as sc
    import numpy as np

    log.info(f"Manual annotation: scoring {len(markers)} cell types ...")
    score_keys = []
    for cell_type, genes in markers.items():
        genes_in_panel = [g for g in genes if g in adata.var_names]
        if not genes_in_panel:
            log.warning(f"  {cell_type}: no marker genes found in panel — skipping.")
            continue
        key = f"score_{cell_type}"
        sc.tl.score_genes(adata, genes_in_panel, score_name=key)
        score_keys.append((key, cell_type))
        log.info(f"  {cell_type}: scored with {len(genes_in_panel)}/{len(genes)} markers.")

    if not score_keys:
        log.warning("No marker genes matched — assigning 'Unknown' to all cells.")
        adata.obs["cell_type"] = "Unknown"
        adata.obs["cell_type_fine"] = "Unknown"
        return adata

    # Assign label of highest score
    import pandas as pd
    score_df = pd.DataFrame(
        {ct: adata.obs[key].values for key, ct in score_keys},
        index=adata.obs_names,
    )
    adata.obs["cell_type"]      = score_df.idxmax(axis=1).values
    adata.obs["cell_type_fine"] = adata.obs["cell_type"]

    n_types = adata.obs["cell_type"].nunique()
    log.info(f"Manual annotation: {n_types} cell types assigned.")
    log.info(adata.obs["cell_type"].value_counts().to_string())
    return adata


def main():
    import spatialdata as sd

    sdata_path       = snakemake.input.sdata           # noqa: F821
    done_path        = snakemake.output.done            # noqa: F821
    method           = snakemake.params.method          # noqa: F821
    ct_model         = snakemake.params.celltypist_model  # noqa: F821
    majority_voting  = snakemake.params.majority_voting  # noqa: F821
    markers          = snakemake.params.markers          # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.table

    if "leiden" not in adata.obs:
        raise RuntimeError("Leiden clusters not found. Run step 05 (reduction) first.")

    if method == "celltypist":
        adata = annotate_celltypist(adata, ct_model, majority_voting)
    elif method == "manual":
        adata = annotate_manual(adata, markers)
    else:
        raise ValueError(f"Unknown annotation method: {method!r}. "
                         "Choose 'celltypist' or 'manual'.")

    sdata.table = adata
    log.info("Writing annotated table to zarr store ...")
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")

    Path(done_path).touch()
    log.info("Step 06 — Annotation: DONE")


main()
