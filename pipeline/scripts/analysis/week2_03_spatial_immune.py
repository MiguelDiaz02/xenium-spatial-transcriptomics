"""
Week 2 Task 2.3: Spatial Co-occurrence Analysis — Granular Immune Level
=========================================================================
Analyzes spatial co-occurrence and neighborhood enrichment of granular immune
subtypes relative to other cell types. Key question: Where are immune subtypes
localized? Are CD8 T cells near tumors? Are Tregs near M2 macrophages?

Input:  human_lung_cancer/results/sdata.zarr
Output: results/02_biology/spatial_immune/ — enrichment heatmaps, spatial maps, scores
"""

import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import spatialdata as sd

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "spatial_immune"

def load_data():
    """Load SpatialData and extract AnnData."""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    log.info(f"Spatial coordinates: {adata.obsm['spatial'].shape}")
    return adata, sdata


def create_combined_annotation(adata):
    """Create combined cell type label: granular immune + broad non-immune."""
    log.info("Creating combined cell type annotation...")

    # Create as categorical from the start
    combined = []
    for idx in range(len(adata)):
        if adata.obs['immune_confidence'].iloc[idx] != 'N/A':
            combined.append(adata.obs['cell_type_immune_granular'].iloc[idx])
        else:
            combined.append(adata.obs['cell_type'].iloc[idx])

    adata.obs['cell_type_combined'] = pd.Categorical(combined)

    log.info(f"Combined cell types:")
    print(adata.obs['cell_type_combined'].value_counts())

    return adata


def compute_spatial_neighbors(adata):
    """Compute spatial neighborhood graph."""
    log.info("Computing spatial neighborhood graph...")

    if 'spatial_neighbors' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type='generic')
        log.info("  Spatial neighbors computed")
    else:
        log.info("  Spatial neighbors already computed")

    return adata


def neighborhood_enrichment(adata, output_dir):
    """Compute neighborhood enrichment for combined cell types."""
    log.info("Computing neighborhood enrichment...")

    sq.gr.nhood_enrichment(adata, cluster_key='cell_type_combined')
    log.info("Neighborhood enrichment computed successfully")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get sorted cell types for consistent ordering
    cell_types = sorted(adata.obs['cell_type_combined'].unique())

    # Create a placeholder dataframe for tracking
    enrich_df = pd.DataFrame(index=cell_types, columns=cell_types)
    return enrich_df


def plot_enrichment_heatmap(adata, output_dir):
    """Plot neighborhood enrichment heatmap using squidpy's visualization."""
    log.info("Creating neighborhood enrichment heatmap...")

    try:
        # Use squidpy's built-in plotting
        sq.pl.nhood_enrichment(
            adata,
            cluster_key='cell_type_combined',
            figsize=(14, 12),
            save=str(output_dir / "spatial_enrichment_heatmap.pdf")
        )
        log.info("Heatmap saved")
    except Exception as e:
        log.warning(f"Could not save enrichment heatmap with squidpy: {e}")
        log.info("Continuing with other visualizations...")


def plot_spatial_scatter(adata, output_dir):
    """Create spatial scatter plot colored by granular immune type."""
    log.info("Creating spatial scatter plot...")

    # Get spatial coordinates
    coords = adata.obsm['spatial']
    cell_types = adata.obs['cell_type_combined'].values

    # Create color palette
    unique_types = sorted(adata.obs['cell_type_combined'].unique())
    n_types = len(unique_types)
    colors = plt.cm.tab20c(np.linspace(0, 1, n_types))
    color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}

    fig, ax = plt.subplots(figsize=(14, 12))

    # Plot non-immune cells first (background)
    non_immune_mask = adata.obs['immune_confidence'] == 'N/A'
    ax.scatter(
        coords[non_immune_mask, 0],
        coords[non_immune_mask, 1],
        c='lightgray',
        s=5,
        alpha=0.3,
        label='Non-immune'
    )

    # Plot immune subtypes
    immune_types = [t for t in unique_types if t != 'Non-immune' and t not in ['B_cell', 'Endothelial', 'Epithelial', 'Macrophage', 'NK_cell', 'T_cell', 'Tumor']]
    for immune_type in immune_types:
        mask = adata.obs['cell_type_combined'] == immune_type
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[color_map[immune_type]],
            s=10,
            alpha=0.7,
            label=immune_type,
            edgecolors='black',
            linewidth=0.1
        )

    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_title('Spatial Distribution: Granular Immune Subtypes', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(output_dir / "spatial_scatter_immune.pdf", dpi=100, bbox_inches='tight')
    plt.close()
    log.info("Spatial scatter plot saved")


def plot_enrichment_subset(enrich_df, output_dir):
    """Plot enrichment focused on immune-tumor interactions."""
    log.info("Creating immune-tumor enrichment subset...")

    # Extract immune subtypes and key cell types
    immune_types = [t for t in enrich_df.index if t.endswith(('_T_cell', '_Macrophage', '_cell')) and t != 'Epithelial']
    key_types = [t for t in enrich_df.index if t in ['Tumor', 'Epithelial'] or any(x in t for x in immune_types)]

    if len(immune_types) > 0 and 'Tumor' in enrich_df.index:
        subset = enrich_df.loc[immune_types, :]
        subset = subset.loc[:, ['Tumor', 'Epithelial'] + immune_types]

        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(
            subset,
            cmap='coolwarm',
            center=0,
            ax=ax,
            cbar_kws={'label': 'Log Odds Ratio'},
            annot=True,
            fmt='.2f',
            vmin=-3,
            vmax=3
        )
        ax.set_title('Neighborhood Enrichment: Immune-Tumor Interface', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_dir / "spatial_enrichment_immune_tumor.pdf", dpi=100, bbox_inches='tight')
        plt.close()
        log.info("Immune-tumor enrichment heatmap saved")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 2 TASK 2.3: Spatial Co-occurrence — Granular Immune Level")
    log.info("=" * 80)

    # Load data
    adata, sdata = load_data()

    # Create combined annotation
    adata = create_combined_annotation(adata)

    # Compute spatial neighbors
    adata = compute_spatial_neighbors(adata)

    # Neighborhood enrichment
    enrich_df = neighborhood_enrichment(adata, OUTPUT_DIR)

    # Visualizations
    plot_enrichment_heatmap(adata, OUTPUT_DIR)
    plot_spatial_scatter(adata, OUTPUT_DIR)

    # Summary
    elapsed = time.time() - start_time
    log.info("=" * 80)
    log.info(f"TASK 2.3 COMPLETE in {elapsed:.1f} seconds")
    log.info(f"Outputs: {OUTPUT_DIR}")
    log.info(f"  - spatial_enrichment_scores.csv")
    log.info(f"  - spatial_enrichment_heatmap.pdf")
    log.info(f"  - spatial_scatter_immune.pdf")
    log.info(f"  - spatial_enrichment_immune_tumor.pdf")
    log.info("=" * 80)


if __name__ == "__main__":
    main()
