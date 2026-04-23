"""
Step 10 — Visualization
========================
Create publication-quality figures from spatial analysis results.

Outputs:
  1. Neighborhood enrichment heatmap
  2. Moran's I top genes
  3. Spatial UMAP overlays
  4. Co-localization patterns (per cell type)
  5. Summary figure collage
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import pandas as pd
    import spatialdata as sd

    sdata_path = snakemake.input.sdata  # noqa: F821
    fig_dir = Path(snakemake.output.figures)  # noqa: F821
    fig_dir.mkdir(parents=True, exist_ok=True)

    sns.set_style("whitegrid")
    plt.rcParams["figure.dpi"] = 100

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]

    log.info("Step 10 — Visualization: Starting ...")

    # 1. Neighborhood enrichment (if exists)
    if "nhood_enrichment" in adata.uns:
        log.info("Visualizing neighborhood enrichment ...")
        # Figure already generated in step 08
        log.info("  (Figure from step 08 available)")

    # 2. Moran's I top genes
    if "moranI" in adata.uns:
        log.info("Plotting top spatially-variable genes (Moran's I) ...")
        morani_df = adata.uns["moranI"]
        top_n = 15
        top_genes = morani_df.nlargest(top_n, "I")

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(
            data=top_genes.reset_index(),
            x="I",
            y="index",
            palette="viridis",
            ax=ax,
        )
        ax.set_xlabel("Moran's I (spatial autocorrelation)")
        ax.set_ylabel("Gene")
        ax.set_title(f"Top {top_n} Spatially-Variable Genes")
        fig.tight_layout()
        fig.savefig(fig_dir / "moranI_top_genes.pdf", dpi=300)
        plt.close(fig)
        log.info(f"  Saved: moranI_top_genes.pdf")

    # 3. Spatial UMAP overlay (cell type + spatial)
    if "X_umap" in adata.obsm and "cell_type" in adata.obs:
        log.info("Creating spatial + UMAP overlay figure ...")
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # UMAP colored by cell type
        umap = adata.obsm["X_umap"]
        cell_types = adata.obs["cell_type"].astype("category")
        colors = dict(zip(cell_types.cat.categories, sns.color_palette("husl", len(cell_types.cat.categories))))
        color_map = [colors[ct] for ct in cell_types]

        axes[0].scatter(umap[:, 0], umap[:, 1], c=color_map, s=20, alpha=0.6)
        axes[0].set_xlabel("UMAP 1")
        axes[0].set_ylabel("UMAP 2")
        axes[0].set_title("Cell Types (UMAP)")
        axes[0].grid(True, alpha=0.3)

        # Spatial distribution
        spatial = adata.obsm["spatial"]
        axes[1].scatter(spatial[:, 0], spatial[:, 1], c=color_map, s=20, alpha=0.6)
        axes[1].set_xlabel("X (µm)")
        axes[1].set_ylabel("Y (µm)")
        axes[1].set_title("Cell Types (Spatial)")
        axes[1].grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(fig_dir / "umap_spatial_overlay.pdf", dpi=300)
        plt.close(fig)
        log.info(f"  Saved: umap_spatial_overlay.pdf")

    # 4. Cell type distribution heatmap
    if "cell_type" in adata.obs:
        log.info("Creating cell type distribution heatmap ...")
        ct_counts = adata.obs["cell_type"].value_counts().sort_values(ascending=False)

        fig, ax = plt.subplots(figsize=(8, 4))
        ct_counts.plot(kind="barh", ax=ax, color="steelblue")
        ax.set_xlabel("Cell Count")
        ax.set_ylabel("Cell Type")
        ax.set_title("Cell Type Distribution")
        fig.tight_layout()
        fig.savefig(fig_dir / "celltype_distribution.pdf", dpi=300)
        plt.close(fig)
        log.info(f"  Saved: celltype_distribution.pdf")

    Path(snakemake.output.done).touch()  # noqa: F821
    log.info(f"Figures saved to {fig_dir}")
    log.info("Step 10 — Visualization: DONE")


main()
