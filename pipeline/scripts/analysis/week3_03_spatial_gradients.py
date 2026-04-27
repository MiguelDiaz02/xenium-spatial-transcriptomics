#!/usr/bin/env python3
"""
Week 3 Phase 3: Task 3.1 - Spatial Gradients
Map immune infiltration gradients (tumor core → edge → healthy tissue)
Identifies where immune cells localize relative to tumor boundary
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import cdist
from scipy.ndimage import distance_transform_edt

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "spatial_gradients"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_spatial_gradients"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load SpatialData"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata

def identify_tumor_boundary(adata):
    """
    Identify tumor cells and compute distance to tumor boundary

    Method: Use spatial convex hull of tumor cells to define boundary
    Returns distance matrix where negative = tumor core, positive = immune infiltrate
    """
    log.info("Identifying tumor boundary...")

    coords = adata.obsm['spatial']
    tumor_mask = adata.obs['cell_type'] == 'Tumor'
    tumor_coords = coords[tumor_mask]

    log.info(f"  Found {tumor_mask.sum()} tumor cells")

    # Compute distance from each cell to nearest tumor cell
    # Negative distance = inside tumor, positive = outside
    distances = np.zeros(len(adata))

    # For tumor cells: distance = 0
    distances[tumor_mask] = 0

    # For non-tumor cells: distance to nearest tumor cell
    non_tumor_coords = coords[~tumor_mask]
    non_tumor_idx = np.where(~tumor_mask)[0]

    if len(tumor_coords) > 0 and len(non_tumor_coords) > 0:
        dist_matrix = cdist(non_tumor_coords, tumor_coords, metric='euclidean')
        min_dist = dist_matrix.min(axis=1)
        distances[non_tumor_idx] = min_dist

    adata.obs['distance_to_tumor'] = distances

    log.info(f"  Distance range: [{distances.min():.1f}, {distances.max():.1f}]")

    return adata

def compute_immune_gradients(adata):
    """
    Compute immune cell density along distance gradient

    Binned analysis: group cells by distance from tumor and count immune cells per bin
    """
    log.info("Computing immune cell gradients...")

    # Define distance bins (from tumor center outward)
    max_dist = adata.obs['distance_to_tumor'].max()
    bins = np.linspace(0, max_dist, 11)  # 10 bins

    # Identify immune cell types (exclude Tumor and Epithelial as "structural")
    immune_types = ['T_cell', 'B_cell', 'Macrophage', 'NK_cell', 'Endothelial', 'B_cell']

    # Binned analysis
    gradient_data = []

    for i in range(len(bins) - 1):
        bin_left = bins[i]
        bin_right = bins[i + 1]
        bin_mid = (bin_left + bin_right) / 2

        # Cells in this distance bin
        bin_mask = (adata.obs['distance_to_tumor'] >= bin_left) & (adata.obs['distance_to_tumor'] < bin_right)
        bin_cells = adata[bin_mask]

        if len(bin_cells) == 0:
            continue

        # Count cell types in this bin
        for immune_type in immune_types:
            count = (bin_cells.obs['cell_type'] == immune_type).sum()
            fraction = count / len(bin_cells) if len(bin_cells) > 0 else 0

            gradient_data.append({
                'distance_from_tumor': bin_mid,
                'distance_bin': f"{bin_left:.0f}-{bin_right:.0f}",
                'cell_type': immune_type,
                'count': count,
                'fraction': fraction,
                'total_cells_in_bin': len(bin_cells)
            })

    gradient_df = pd.DataFrame(gradient_data)
    log.info(f"  Computed gradients for {len(gradient_df)} cell type × distance combinations")

    return adata, gradient_df

def plot_immune_gradient_curve(gradient_df):
    """Plot immune cell fraction as function of distance from tumor"""
    fig, ax = plt.subplots(figsize=(12, 6))

    immune_types = gradient_df['cell_type'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(immune_types)))

    for immune_type, color in zip(immune_types, colors):
        type_data = gradient_df[gradient_df['cell_type'] == immune_type].sort_values('distance_from_tumor')
        ax.plot(type_data['distance_from_tumor'], type_data['fraction'],
               marker='o', linewidth=2, markersize=8, label=immune_type, color=color)

    ax.set_xlabel('Distance from Tumor (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fraction of Cells in Bin', fontsize=12, fontweight='bold')
    ax.set_title('Immune Cell Gradients: Spatial Distribution from Tumor Boundary',
                fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    return fig

def plot_distance_heatmap(gradient_df):
    """Plot cell type abundance as heatmap across distance bins"""
    fig, ax = plt.subplots(figsize=(14, 6))

    # Pivot for heatmap
    pivot_data = gradient_df.pivot_table(
        values='fraction',
        index='cell_type',
        columns='distance_bin',
        aggfunc='first'
    )

    sns.heatmap(pivot_data, annot=pivot_data.round(2), fmt='g', cmap='YlOrRd',
               cbar_kws={'label': 'Fraction of Cells'}, ax=ax, linewidths=0.5)

    ax.set_xlabel('Distance from Tumor (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_title('Immune Cell Distribution: Heatmap by Distance from Tumor Boundary',
                fontsize=13, fontweight='bold')

    plt.tight_layout()
    return fig

def plot_spatial_distance_map(adata):
    """Plot spatial map colored by distance from tumor"""
    fig, ax = plt.subplots(figsize=(14, 12))

    coords = adata.obsm['spatial']
    distances = adata.obs['distance_to_tumor'].values

    # Color by distance: red = tumor core, blue = far from tumor
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=distances, cmap='RdYlBu_r',
                        s=10, alpha=0.6, edgecolors='none')

    # Highlight tumor cells
    tumor_mask = adata.obs['cell_type'] == 'Tumor'
    ax.scatter(coords[tumor_mask, 0], coords[tumor_mask, 1],
              c='red', s=5, alpha=0.8, label='Tumor', edgecolors='none')

    ax.set_xlabel('X coordinate (μm)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=11, fontweight='bold')
    ax.set_title('Spatial Distance Map: Distance from Tumor Boundary',
                fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    cbar = plt.colorbar(scatter, ax=ax, label='Distance to Tumor (μm)')
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    return fig

def save_results(adata, gradient_df):
    """Save gradient analysis results"""
    log.info("Saving gradient analysis results...")

    # Save distance-to-tumor assignment
    adata.obs['distance_to_tumor'].to_csv(OUTPUT_DIR / "distance_to_tumor.csv")

    # Save gradient table
    gradient_df.to_csv(OUTPUT_DIR / "immune_gradients.csv", index=False)

    # Summary statistics
    summary = pd.DataFrame({
        'metric': ['tumor_cells', 'max_distance', 'median_distance_immune', 'cell_types'],
        'value': [
            (adata.obs['cell_type'] == 'Tumor').sum(),
            adata.obs['distance_to_tumor'].max(),
            adata.obs[adata.obs['cell_type'] != 'Tumor']['distance_to_tumor'].median(),
            adata.obs['cell_type'].nunique()
        ]
    })
    summary.to_csv(OUTPUT_DIR / "gradient_summary.csv", index=False)

    log.info(f"✓ Results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 3 - TASK 3.1: SPATIAL GRADIENTS")
    log.info("Immune infiltration mapping (tumor core → edge → healthy)")
    log.info("=" * 80)

    try:
        # Load data
        adata = load_data()

        # Identify tumor boundary
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Identify Tumor Boundary")
        log.info("=" * 80)
        adata = identify_tumor_boundary(adata)

        # Compute gradients
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Compute Immune Cell Gradients")
        log.info("=" * 80)
        adata, gradient_df = compute_immune_gradients(adata)

        # Generate figures
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Generate Visualizations")
        log.info("=" * 80)

        log.info("  Generating gradient curve plot...")
        fig1 = plot_immune_gradient_curve(gradient_df)
        fig1.savefig(FIGURES_DIR / "Fig1_Immune_Gradient_Curves.png", dpi=300, bbox_inches='tight')
        plt.close(fig1)

        log.info("  Generating distance heatmap...")
        fig2 = plot_distance_heatmap(gradient_df)
        fig2.savefig(FIGURES_DIR / "Fig2_Distance_Heatmap.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)

        log.info("  Generating spatial distance map...")
        fig3 = plot_spatial_distance_map(adata)
        fig3.savefig(FIGURES_DIR / "Fig3_Spatial_Distance_Map.png", dpi=300, bbox_inches='tight')
        plt.close(fig3)

        # Save results
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Save Results")
        log.info("=" * 80)
        save_results(adata, gradient_df)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ TASK 3.1 COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info(f"✓ Figures directory: {FIGURES_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
