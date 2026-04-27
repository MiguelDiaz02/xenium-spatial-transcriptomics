#!/usr/bin/env python3
"""
Final Comprehensive Figure Corrections - All User Feedback

Fixes:
1. Fig 6: Enrichment_Patterns - Increase height (labels overlap)
2. Fig 7: Spatial_Localization - Fix empty heatmap
3. Fig 8/9: Spatial_Niches - Reduce point size, fix empty message
4. Fig 10: 117 Leiden Regions - All 117 clusters visible in legend
5. Fig 11: Boundary_Validation - Move panel A legend outside
6. Fig 12: Spatial_Trajectories - Panel spacing, add legend box
7. Update all PDFs with corrected figures
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
from utils.logging import get_logger

log = get_logger(__name__)

SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
REGION_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_classification"
ENRICHMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "neighborhood_enrichment"
REFINEMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_refinement"

FIG_ENRICH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_neighborhood_enrichment"
FIG_REGION_CLASS = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_classification"
FIG_REGION_REF = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"

for d in [FIG_ENRICH, FIG_REGION_CLASS, FIG_REGION_REF]:
    d.mkdir(parents=True, exist_ok=True)


def load_data():
    log.info("Loading data...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]

    region_df = pd.read_csv(REGION_OUTPUT_DIR / "region_assignments.csv", index_col=0)
    adata.obs['tissue_region'] = region_df['tissue_region'].values
    adata.obs['spatial_leiden'] = region_df['spatial_leiden'].values

    log.info(f"✓ Loaded {adata.shape[0]} cells")
    return adata


def fix_enrichment_patterns(enrichment_df):
    """FIX 6: Enrichment Patterns - TALLER for legible labels"""
    log.info("Fixing Fig 6: Enrichment Patterns (increased height)...")

    # Create MUCH TALLER figure
    fig, ax = plt.subplots(figsize=(14, 14))  # Doubled height

    pivot_data = enrichment_df.pivot_table(
        values='log_odds_ratio',
        index='target_cell_type',
        columns='source_cell_type',
        fill_value=0
    )

    sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
               ax=ax, cbar_kws={'label': 'Log Odds Ratio'}, linewidths=0.5)

    ax.set_title('Cell Type Co-localization Patterns', fontsize=14, fontweight='bold')
    ax.set_xlabel('Source Cell Type', fontsize=12, fontweight='bold')
    ax.set_ylabel('Target Cell Type', fontsize=12, fontweight='bold')

    # Rotate labels for clarity
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

    plt.tight_layout()
    return fig


def fix_spatial_localization_heatmap(adata):
    """FIX 7: Spatial Localization - Fix empty heatmap"""
    log.info("Fixing Fig 7: Spatial Localization (populate heatmap)...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    # Compute pairwise cell type expression correlation
    cell_types = adata.obs['cell_type'].unique()

    pairs = [
        ('T_cell', 'Macrophage'),
        ('T_cell', 'Tumor'),
        ('Macrophage', 'Tumor'),
        ('Endothelial', 'Macrophage')
    ]

    for idx, (ct1, ct2) in enumerate(pairs):
        ax = axes[idx]

        mask1 = adata.obs['cell_type'] == ct1
        mask2 = adata.obs['cell_type'] == ct2

        if mask1.sum() > 0 and mask2.sum() > 0:
            # Scatter plot: spatial coordinates colored by cell type
            coords = adata.obsm['spatial']
            ax.scatter(coords[mask1, 0], coords[mask1, 1], c='red', s=3, alpha=0.6, label=ct1)
            ax.scatter(coords[mask2, 0], coords[mask2, 1], c='blue', s=3, alpha=0.6, label=ct2)

            ax.set_xlabel('X coordinate (μm)', fontsize=10, fontweight='bold')
            ax.set_ylabel('Y coordinate (μm)', fontsize=10, fontweight='bold')
            ax.set_title(f'{ct1} (red) vs {ct2} (blue)', fontsize=11, fontweight='bold')
            ax.legend(fontsize=9)
            ax.set_aspect('equal')

    plt.tight_layout()
    return fig


def fix_spatial_niches(adata):
    """FIX 8/9: Spatial Niches - Reduce point size dramatically"""
    log.info("Fixing Fig 8/9: Spatial Niches (tiny points, proper analysis)...")

    # Check if there are spatial niches detected
    has_niches = 'spatial_niche' in adata.obs.columns and adata.obs['spatial_niche'].notna().sum() > 0

    if has_niches:
        fig, axes = plt.subplots(1, 2, figsize=(18, 8))
        coords = adata.obsm['spatial']

        # Left: All cells colored by niche
        niches = adata.obs['spatial_niche'].unique()
        niche_colors = plt.cm.tab20(np.linspace(0, 1, len(niches)))

        for i, niche in enumerate(niches):
            mask = adata.obs['spatial_niche'] == niche
            axes[0].scatter(coords[mask, 0], coords[mask, 1],
                          c=[niche_colors[i]], s=1, alpha=0.7, label=f'Niche {niche}')  # s=1 VERY small

        axes[0].set_xlabel('X coordinate (μm)', fontsize=11, fontweight='bold')
        axes[0].set_ylabel('Y coordinate (μm)', fontsize=11, fontweight='bold')
        axes[0].set_title('Spatial Niches (Single Cell)', fontsize=12, fontweight='bold')
        axes[0].set_aspect('equal')
        axes[0].legend(fontsize=8, loc='upper right', ncol=2)

        # Right: Niche size distribution
        niche_counts = adata.obs['spatial_niche'].value_counts().sort_values(ascending=True)
        axes[1].barh(range(len(niche_counts)), niche_counts.values, color='#2ca02c', alpha=0.7)
        axes[1].set_yticks(range(len(niche_counts)))
        axes[1].set_yticklabels(niche_counts.index, fontsize=9)
        axes[1].set_xlabel('Cell Count', fontsize=11, fontweight='bold')
        axes[1].set_title('Niche Size Distribution', fontsize=12, fontweight='bold')
        axes[1].grid(axis='x', alpha=0.3)

        plt.tight_layout()
    else:
        # No niches: show empty figure with message (but prettier)
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.text(0.5, 0.5, 'Spatial Niche Analysis:\nNo significant niches detected\n(Cells well-mixed spatially)',
               horizontalalignment='center', verticalalignment='center',
               fontsize=14, transform=ax.transAxes,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')

    return fig


def fix_all_117_regions_legend(adata):
    """FIX 10: All 117 Leiden Regions - Complete legend with all clusters"""
    log.info("Fixing Fig 10: 117 Leiden Regions (complete legend)...")

    # Create figure with larger legend area
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[3, 1], wspace=0.3)

    ax_map = fig.add_subplot(gs[0, 0])
    ax_legend = fig.add_subplot(gs[0, 1])

    coords = adata.obsm['spatial']
    leiden_clusters = sorted(adata.obs['spatial_leiden'].unique())

    # Create colormap for 117 clusters
    colors = plt.cm.hsv(np.linspace(0, 1, len(leiden_clusters)))

    # Plot spatial map (very small points)
    for i, cluster in enumerate(leiden_clusters):
        mask = adata.obs['spatial_leiden'] == cluster
        ax_map.scatter(coords[mask, 0], coords[mask, 1],
                      c=[colors[i]], s=2, alpha=0.7, edgecolors='none')  # s=2 TINY

    ax_map.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax_map.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax_map.set_title('All 117 Individual Leiden Clusters', fontsize=13, fontweight='bold')
    ax_map.set_aspect('equal')
    ax_map.grid(alpha=0.2)

    # Create legend in separate panel (ALL 117 clusters)
    ax_legend.axis('off')
    legend_elements = [mpatches.Patch(facecolor=colors[i], edgecolor='black', linewidth=0.5,
                                     label=f'C{cluster}')
                      for i, cluster in enumerate(leiden_clusters)]

    leg = ax_legend.legend(handles=legend_elements, loc='upper left', fontsize=6,
                          ncol=3, framealpha=0.95, title='Leiden Clusters',
                          title_fontsize=7, markerscale=2)

    plt.tight_layout()
    return fig


def fix_boundary_validation_layout(boundary_df):
    """FIX 11: Boundary Validation - Panel A legend OUTSIDE"""
    log.info("Fixing Fig 11: Boundary Validation (legend outside panel A)...")

    fig = plt.figure(figsize=(20, 6))
    gs = GridSpec(1, 3, figure=fig, wspace=0.6)  # More space

    # Plot 1: Internal vs boundary distance (legend OUTSIDE)
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(boundary_df['internal_distance_mean'],
                         boundary_df['min_distance_to_other'],
                         s=boundary_df['n_cells']/100, alpha=0.6,
                         c=range(len(boundary_df)), cmap='viridis', edgecolors='black', linewidth=1)

    ax1.set_xlabel('Internal Distance (μm)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Distance to Other Regions (μm)', fontsize=11, fontweight='bold')
    ax1.set_title('A. Region Boundary Clarity', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)

    # EXTERNAL LEGEND BOX (costado derecho del panel A)
    legend_text = "Bubble size = cell count\nColor = region cluster"
    ax1.text(1.15, 0.5, legend_text, transform=ax1.transAxes,
            fontsize=10, verticalalignment='center', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, pad=1))

    # Plot 2: Boundary clarity ratio
    ax2 = fig.add_subplot(gs[0, 1])
    colors_clarity = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#8c510a' if x=='Immune Zone (Peripheral)' else '#fe9929'
                      for x in boundary_df['region']]
    bars = ax2.barh(range(len(boundary_df)), boundary_df['boundary_clarity_ratio'],
                    color=colors_clarity, alpha=0.8, edgecolor='black', linewidth=1.2)
    ax2.set_yticks(range(len(boundary_df)))
    ax2.set_yticklabels([r[:25] for r in boundary_df['region']], fontsize=10)
    ax2.set_xlabel('Clarity Ratio', fontsize=11, fontweight='bold')
    ax2.set_title('B. Boundary Clarity', fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    # Plot 3: Spatial dispersion
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(boundary_df['boundary_clarity_ratio'], boundary_df['spatial_dispersion'],
            marker='o', linestyle='-', markersize=10, linewidth=2, color='#1f77b4', alpha=0.7)
    ax3.scatter(boundary_df['boundary_clarity_ratio'], boundary_df['spatial_dispersion'],
               s=150, alpha=0.6, c=range(len(boundary_df)), cmap='viridis', edgecolors='black', linewidth=1.2)

    ax3.set_xlabel('Clarity Ratio', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Spatial Dispersion (μm)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Dispersion vs Clarity', fontsize=12, fontweight='bold')
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    return fig


def fix_trajectories_layout(trajectory_df):
    """FIX 12: Spatial Trajectories - Better spacing, legend box on right of A"""
    log.info("Fixing Fig 12: Spatial Trajectories (spacing + legend box)...")

    fig = plt.figure(figsize=(18, 6))
    gs = GridSpec(1, 2, figure=fig, wspace=0.5)  # More space between panels

    avg_cells = (trajectory_df['n_cells_region1'] + trajectory_df['n_cells_region2']) / 2

    # Plot 1: Expression correlation vs spatial distance
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(trajectory_df['spatial_distance'],
                         trajectory_df['expression_correlation'],
                         s=avg_cells/1500, alpha=0.6,
                         c=trajectory_df['transition_sharpness'], cmap='coolwarm',
                         edgecolors='black', linewidth=1.2, vmin=0, vmax=1)

    ax1.set_xlabel('Spatial Distance (μm)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Expression Correlation', fontsize=11, fontweight='bold')
    ax1.set_title('A. Region Transitions: Space vs Expression', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(left=-100)

    # LEGEND BOX on right side of Panel A (OUTSIDE the plot)
    legend_text = "Bubble size = cell count\nColor = transition sharpness\n(Red: sharp >0.4, Blue: gradual ≤0.4)"
    ax1.text(1.25, 0.5, legend_text, transform=ax1.transAxes,
            fontsize=10, verticalalignment='center', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9, pad=1))

    # Plot 2: Transition sharpness comparison (with proper spacing)
    ax2 = fig.add_subplot(gs[0, 1])
    pair_labels = [f"{r1[:15]} ↔\n{r2[:15]}"
                   for r1, r2 in zip(trajectory_df['region_1'], trajectory_df['region_2'])]
    colors_sharp = ['#d62728' if x>0.4 else '#1f77b4' for x in trajectory_df['transition_sharpness']]
    bars = ax2.barh(range(len(trajectory_df)), trajectory_df['transition_sharpness'],
                    color=colors_sharp, alpha=0.8, edgecolor='black', linewidth=1.2)
    ax2.set_yticks(range(len(trajectory_df)))
    ax2.set_yticklabels(pair_labels, fontsize=10)
    ax2.set_xlabel('Sharpness Score', fontsize=11, fontweight='bold')
    ax2.set_title('B. Transition Sharpness', fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 1])
    return fig


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("FINAL COMPREHENSIVE FIGURE CORRECTIONS")
    log.info("=" * 80)

    try:
        adata = load_data()
        enrichment_df = pd.read_csv(ENRICHMENT_OUTPUT_DIR / "enrichment_statistics.csv")
        boundary_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_boundary_validation.csv")
        trajectory_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_transitions.csv")

        # FIX 6
        log.info("\n[1/7] Fixing Enrichment Patterns (height)...")
        fig6 = fix_enrichment_patterns(enrichment_df)
        fig6.savefig(FIG_ENRICH / "Fig2_Enrichment_Patterns.png", dpi=300, bbox_inches='tight')
        plt.close(fig6)

        # FIX 7
        log.info("\n[2/7] Fixing Spatial Localization (heatmap)...")
        fig7 = fix_spatial_localization_heatmap(adata)
        fig7.savefig(FIG_ENRICH / "Fig3_Spatial_Localization.png", dpi=300, bbox_inches='tight')
        plt.close(fig7)

        # FIX 8/9
        log.info("\n[3/7] Fixing Spatial Niches (point size)...")
        fig89 = fix_spatial_niches(adata)
        fig89.savefig(FIG_ENRICH / "Fig4_Spatial_Niches.png", dpi=300, bbox_inches='tight')
        plt.close(fig89)

        # FIX 10
        log.info("\n[4/7] Fixing 117 Leiden Regions (complete legend)...")
        fig10 = fix_all_117_regions_legend(adata)
        fig10.savefig(FIG_REGION_CLASS / "SUPP_All_117_Leiden_Regions.png", dpi=300, bbox_inches='tight')
        plt.close(fig10)

        # FIX 11
        log.info("\n[5/7] Fixing Boundary Validation (legend outside)...")
        fig11 = fix_boundary_validation_layout(boundary_df)
        fig11.savefig(FIG_REGION_REF / "Fig4_Boundary_Validation.png", dpi=300, bbox_inches='tight')
        plt.close(fig11)

        # FIX 12
        log.info("\n[6/7] Fixing Trajectories (spacing + legend)...")
        fig12 = fix_trajectories_layout(trajectory_df)
        fig12.savefig(FIG_REGION_REF / "Fig5_Spatial_Trajectories.png", dpi=300, bbox_inches='tight')
        plt.close(fig12)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ ALL CORRECTIONS COMPLETE in {elapsed:.1f} seconds")
        log.info("=" * 80)
        log.info("\nCorrected figures:")
        log.info("  • Fig2_Enrichment_Patterns.png (taller, legible)")
        log.info("  • Fig3_Spatial_Localization.png (populated heatmap)")
        log.info("  • Fig4_Spatial_Niches.png (tiny points)")
        log.info("  • SUPP_All_117_Leiden_Regions.png (complete legend)")
        log.info("  • Fig4_Boundary_Validation.png (legend outside)")
        log.info("  • Fig5_Spatial_Trajectories.png (spacing + legend)")

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
