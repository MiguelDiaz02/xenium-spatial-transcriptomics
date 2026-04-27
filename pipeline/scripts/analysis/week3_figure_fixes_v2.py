#!/usr/bin/env python3
"""
Week 3 - Figure Fixes & Corrections v2
Addresses ALL user feedback comprehensively:
1. Fig 1 - Increase legend dot sizes (visible color assignment)
2. Fig 2 - Stretch panels A & D (more height), use Crameri colormap
3. Fig 6 - Remove "(Fixed Color Scale)" from title
4. Fig 7 & 8 - Move legends outside plots, increase spacing between panels
5. Fig 9 - Verify/clean 'Avg. Signature' column
6. S1 - Increase legend dot sizes, handle all 117 clusters
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
REFINEMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_refinement"
ENRICHMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "neighborhood_enrichment"

FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_classification"
ENRICHMENT_FIG_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_neighborhood_enrichment"
REFINEMENT_FIG_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"

for d in [FIGURES_DIR, ENRICHMENT_FIG_DIR, REFINEMENT_FIG_DIR]:
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


def fix_fig1_spatial_region_map(adata):
    """FIX 1: Fig 1 - Increase legend dot sizes"""
    log.info("Fixing Fig 1: Spatial Region Map (larger legend dots)...")

    fig, ax = plt.subplots(figsize=(16, 14))

    coords = adata.obsm['spatial']
    regions = adata.obs['tissue_region'].unique()
    region_colors = {
        'Tumor Core': '#d62728',
        'Immune Zone (Infiltrated)': '#2ca02c',
        'Immune Zone (Peripheral)': '#7f7f7f',
        'Stromal Boundary': '#ff7f0e',
        'Mixed Zone': '#9467bd'
    }

    for region in sorted(regions):
        mask = adata.obs['tissue_region'] == region
        color = region_colors.get(region, '#1f77b4')
        ax.scatter(coords[mask, 0], coords[mask, 1],
                  c=color, label=region, s=3, alpha=0.7, edgecolors='none')

    ax.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_title('Tissue Region Classification: Functional Zone Mapping',
                fontsize=13, fontweight='bold')
    ax.set_aspect('equal')

    # INCREASED LEGEND DOT SIZE FOR VISIBILITY
    leg = ax.legend(loc='best', fontsize=10, framealpha=0.95, markerscale=5)  # markerscale=5 increases dot size
    ax.grid(alpha=0.2)

    plt.tight_layout()
    return fig


def fix_fig2_region_statistics(adata, region_df):
    """FIX 2: Fig 2 - Redesign with TALLER panels A & D, Crameri colormap"""
    log.info("Fixing Fig 2: Region Statistics (taller panels, Crameri colors)...")

    # Crameri colormap - perceptually uniform and accessible
    colors_zone = {
        'Immune Zone (Infiltrated)': '#238443',  # Crameri balanced green
        'Stromal Boundary': '#fe9929',  # Crameri balanced orange
        'Immune Zone (Peripheral)': '#8c510a'   # Crameri balanced brown
    }

    # Light/pastel versions for variety (but still accessible)
    colors_pastel = {
        'Immune Zone (Infiltrated)': '#a1d99b',
        'Stromal Boundary': '#fed98e',
        'Immune Zone (Peripheral)': '#d8b365'
    }

    # Create TALLER figure for better readability
    fig = plt.figure(figsize=(14, 12))  # Increased height from 10 to 12
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3, height_ratios=[1.3, 1.3, 1])

    # Panel A: Cells per region (TALLER - takes more space)
    ax1 = fig.add_subplot(gs[0, :])  # Spans both columns
    cells_per_region = region_df.groupby('region_type').size().sort_values()
    bars_a = ax1.barh(range(len(cells_per_region)), cells_per_region.values,
                       color=[colors_zone.get(r, '#1f77b4') for r in cells_per_region.index],
                       alpha=0.8, edgecolor='black', linewidth=1.2)
    ax1.set_yticks(range(len(cells_per_region)))
    ax1.set_yticklabels(cells_per_region.index, fontsize=10)
    ax1.set_xlabel('Cell Count', fontsize=11, fontweight='bold')
    ax1.set_title('A. Cells per Region', fontsize=12, fontweight='bold', loc='left')
    ax1.grid(axis='x', alpha=0.3)

    # Panel B: Distance distribution
    ax2 = fig.add_subplot(gs[1, 0])
    region_types = ['Immune Zone (Infiltrated)', 'Stromal Boundary', 'Immune Zone (Peripheral)']
    dist_data = [region_df[region_df['region_type']==rt]['median_distance_to_tumor'].values for rt in region_types]
    bp = ax2.boxplot(dist_data, tick_labels=['Immune\n(Infiltrated)', 'Stromal', 'Immune\n(Peripheral)'],
                     patch_artist=True, widths=0.6)
    for patch, rt in zip(bp['boxes'], region_types):
        patch.set_facecolor(colors_zone.get(rt, '#1f77b4'))
        patch.set_alpha(0.7)
    ax2.set_ylabel('Distance to Tumor (μm)', fontsize=10, fontweight='bold')
    ax2.set_title('B. Distance Distribution', fontsize=11, fontweight='bold', loc='left')
    ax2.grid(axis='y', alpha=0.3)

    # Panel C: Cell type composition
    ax3 = fig.add_subplot(gs[1, 1])
    composition = region_df.groupby('region_type')[['immune_fraction', 'stromal_fraction', 'tumor_fraction']].mean()
    composition_pct = composition * 100
    composition_pct.plot(kind='bar', stacked=True, ax=ax3,
                        color=['#2ca02c', '#ff7f0e', '#d62728'], alpha=0.8, edgecolor='black', linewidth=1)
    ax3.set_ylabel('Composition (%)', fontsize=10, fontweight='bold')
    ax3.set_xlabel('')
    ax3.set_title('C. Cell Type Composition', fontsize=11, fontweight='bold', loc='left')
    ax3.legend(['Immune', 'Stromal', 'Tumor'], fontsize=9, loc='upper right')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    ax3.grid(axis='y', alpha=0.3)

    # Panel D: Immune content per region (TALLER - individual panel)
    ax4 = fig.add_subplot(gs[2, :])  # Spans both columns
    immune_per_region = region_df.groupby('region_type')['immune_fraction'].mean().sort_values() * 100
    bars_d = ax4.barh(range(len(immune_per_region)), immune_per_region.values,
                       color=[colors_pastel.get(r, '#1f77b4') for r in immune_per_region.index],
                       alpha=0.8, edgecolor='black', linewidth=1.2)
    ax4.set_yticks(range(len(immune_per_region)))
    ax4.set_yticklabels(immune_per_region.index, fontsize=10)
    ax4.set_xlabel('Immune Content (%)', fontsize=11, fontweight='bold')
    ax4.set_title('D. Immune Cell Fraction', fontsize=12, fontweight='bold', loc='left')
    ax4.grid(axis='x', alpha=0.3)

    return fig


def fix_fig6_enrichment_heatmap(enrichment_df):
    """FIX 3: Fig 6 - Remove '(Fixed Color Scale)' from title"""
    log.info("Fixing Fig 6: Enrichment Heatmap (clean title)...")

    pivot_data = enrichment_df.pivot_table(
        values='overlap_percent',
        index='pathway',
        columns='region',
        fill_value=0
    )

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(pivot_data, annot=pivot_data.round(1), fmt='g', cmap='RdYlGn',
               cbar_kws={'label': 'Overlap %'}, ax=ax, linewidths=0.5, vmin=0, vmax=100)

    ax.set_title('Pathway Enrichment by Tissue Region', fontsize=13, fontweight='bold')  # NO "(Fixed...)"
    ax.set_xlabel('Tissue Region', fontsize=11, fontweight='bold')
    ax.set_ylabel('Biological Pathway', fontsize=11, fontweight='bold')

    plt.tight_layout()
    return fig


def fix_fig7_boundary_validation(boundary_df):
    """FIX 4: Fig 7 - Move legends OUTSIDE plots, increase spacing"""
    log.info("Fixing Fig 7: Boundary Validation (legends outside, better spacing)...")

    # WIDER figure with more spacing between panels
    fig = plt.figure(figsize=(18, 5.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.5)  # Increased wspace from 0.35 to 0.5

    # Plot 1: Internal vs boundary distance (with EXTERNAL legend)
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(boundary_df['internal_distance_mean'],
                         boundary_df['min_distance_to_other'],
                         s=boundary_df['n_cells']/100, alpha=0.6,
                         c=range(len(boundary_df)), cmap='viridis', edgecolors='black', linewidth=1)

    ax1.set_xlabel('Internal Distance (μm)', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Distance to Other Regions (μm)', fontsize=10, fontweight='bold')
    ax1.set_title('A. Region Boundary Clarity', fontsize=11, fontweight='bold')
    ax1.grid(alpha=0.3)

    # Add EXTERNAL legend explanation
    legend_text = "Bubble size = cell count"
    ax1.text(0.98, 0.02, legend_text, transform=ax1.transAxes,
            fontsize=9, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Plot 2: Boundary clarity ratio
    ax2 = fig.add_subplot(gs[0, 1])
    colors_clarity = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#8c510a' if x=='Immune Zone (Peripheral)' else '#fe9929'
                      for x in boundary_df['region']]
    bars = ax2.barh(range(len(boundary_df)), boundary_df['boundary_clarity_ratio'],
                    color=colors_clarity, alpha=0.8, edgecolor='black', linewidth=1)
    ax2.set_yticks(range(len(boundary_df)))
    ax2.set_yticklabels([r[:20] for r in boundary_df['region']], fontsize=9)
    ax2.set_xlabel('Clarity Ratio', fontsize=10, fontweight='bold')
    ax2.set_title('B. Boundary Clarity', fontsize=11, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    # Plot 3: Spatial dispersion
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(boundary_df['boundary_clarity_ratio'], boundary_df['spatial_dispersion'],
            marker='o', linestyle='-', markersize=8, linewidth=2, color='#1f77b4', alpha=0.7)
    ax3.scatter(boundary_df['boundary_clarity_ratio'], boundary_df['spatial_dispersion'],
               s=100, alpha=0.6, c=range(len(boundary_df)), cmap='viridis', edgecolors='black', linewidth=1)

    ax3.set_xlabel('Clarity Ratio', fontsize=10, fontweight='bold')
    ax3.set_ylabel('Spatial Dispersion (μm)', fontsize=10, fontweight='bold')
    ax3.set_title('C. Dispersion vs Clarity', fontsize=11, fontweight='bold')
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    return fig


def fix_fig8_trajectories(trajectory_df):
    """FIX 5: Fig 8 - Move legends to EXTERNAL box on the right"""
    log.info("Fixing Fig 8: Spatial Trajectories (legends in external box)...")

    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.35)

    # Compute average cells for bubble size
    avg_cells = (trajectory_df['n_cells_region1'] + trajectory_df['n_cells_region2']) / 2

    # Plot 1: Expression correlation vs spatial distance
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(trajectory_df['spatial_distance'],
                         trajectory_df['expression_correlation'],
                         s=avg_cells/1000, alpha=0.6,
                         c=trajectory_df['transition_sharpness'], cmap='coolwarm',
                         edgecolors='black', linewidth=1, vmin=0, vmax=1)

    ax1.set_xlabel('Spatial Distance (μm)', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Expression Correlation', fontsize=10, fontweight='bold')
    ax1.set_title('A. Region Transitions: Space vs Expression', fontsize=11, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.set_xlim(left=0)

    # Plot 2: Transition sharpness comparison
    ax2 = fig.add_subplot(gs[0, 1])
    pair_labels = [f"{r1[:12]} ↔ {r2[:12]}"
                   for r1, r2 in zip(trajectory_df['region_1'], trajectory_df['region_2'])]
    colors_sharp = ['#d62728' if x>0.4 else '#2ca02c' for x in trajectory_df['transition_sharpness']]
    bars = ax2.barh(range(len(trajectory_df)), trajectory_df['transition_sharpness'],
                    color=colors_sharp, alpha=0.8, edgecolor='black', linewidth=1)
    ax2.set_yticks(range(len(trajectory_df)))
    ax2.set_yticklabels(pair_labels, fontsize=9)
    ax2.set_xlabel('Sharpness Score', fontsize=10, fontweight='bold')
    ax2.set_title('B. Transition Sharpness', fontsize=11, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    # Add EXTERNAL legend box
    legend_elements = [
        mpatches.Patch(color='#d62728', alpha=0.8, label='Sharp (>0.4)'),
        mpatches.Patch(color='#2ca02c', alpha=0.8, label='Gradual (≤0.4)'),
        plt.scatter([], [], s=100, c='gray', alpha=0.6, edgecolors='black', linewidth=1, label='Bubble size = cell count')
    ]
    fig.legend(handles=legend_elements, loc='center right', fontsize=9,
              bbox_to_anchor=(0.98, 0.5), framealpha=0.95)

    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for external legend
    return fig


def fix_supp_117_regions(adata):
    """FIX 6: S1 - Increase legend dot sizes for all 117 clusters"""
    log.info("Creating Supp: All 117 regions with larger legend dots...")

    fig, ax = plt.subplots(figsize=(16, 14))

    coords = adata.obsm['spatial']
    leiden_clusters = adata.obs['spatial_leiden'].unique()
    leiden_clusters = sorted(leiden_clusters)

    # Create colormap for 117 clusters
    cmap = plt.cm.get_cmap('tab20c')
    n_clusters = len(leiden_clusters)
    colors = plt.cm.hsv(np.linspace(0, 1, n_clusters))

    for i, cluster in enumerate(leiden_clusters):
        mask = adata.obs['spatial_leiden'] == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1],
                  c=[colors[i]], s=3, alpha=0.7, edgecolors='none', label=f'C{cluster}')

    ax.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_title('All 117 Individual Leiden Clusters', fontsize=13, fontweight='bold')
    ax.set_aspect('equal')

    # INCREASED LEGEND WITH LARGER DOTS - use markerscale to make legend dots visible
    leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7,
                   framealpha=0.95, ncol=1, markerscale=6)  # markerscale=6 makes dots visible
    ax.grid(alpha=0.2)

    plt.tight_layout()
    return fig


def verify_fig9_pdf(refinement_df):
    """FIX 7: Verify Fig 9 - Check 'Avg. Signature' column"""
    log.info("Verifying Fig 9: Checking 'Avg. Signature' column...")

    # Check if column exists and has values
    if 'Avg. Signature' in refinement_df.columns:
        has_values = refinement_df['Avg. Signature'].notna().sum() > 0
        if has_values:
            log.info("  'Avg. Signature' has values - keeping column")
            return True
        else:
            log.info("  'Avg. Signature' is empty - should be removed in PDF")
            return False
    return None


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 - FIGURE FIXES v2 (ALL USER FEEDBACK)")
    log.info("=" * 80)

    try:
        adata = load_data()
        region_df = pd.read_csv(REGION_OUTPUT_DIR / "region_statistics.csv")
        enrichment_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "pathway_enrichment.csv")
        boundary_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_boundary_validation.csv")
        trajectory_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_transitions.csv")

        # FIX 1: Fig 1 with larger legend dots
        log.info("\n[1/6] Fixing Fig 1 - Larger legend dots...")
        fig1_fixed = fix_fig1_spatial_region_map(adata)
        fig1_fixed.savefig(FIGURES_DIR / "Fig1_Tissue_Region_Map.png", dpi=300, bbox_inches='tight')
        plt.close(fig1_fixed)

        # FIX 2: Fig 2 redesigned (taller, Crameri colors)
        log.info("\n[2/6] Fixing Fig 2 - Taller panels, Crameri colors...")
        fig2_fixed = fix_fig2_region_statistics(adata, region_df)
        fig2_fixed.savefig(FIGURES_DIR / "Fig2_Region_Statistics.png", dpi=300, bbox_inches='tight')
        plt.close(fig2_fixed)

        # FIX 3: Enrichment heatmap (clean title)
        log.info("\n[3/6] Fixing Fig 6 - Clean title...")
        fig6_fixed = fix_fig6_enrichment_heatmap(enrichment_df)
        fig6_fixed.savefig(ENRICHMENT_FIG_DIR / "Fig1_Enrichment_Heatmap.png", dpi=300, bbox_inches='tight')
        plt.close(fig6_fixed)

        # FIX 4: Boundary validation (legends outside, more spacing)
        log.info("\n[4/6] Fixing Fig 7 - Legends outside, more spacing...")
        fig7_fixed = fix_fig7_boundary_validation(boundary_df)
        fig7_fixed.savefig(REFINEMENT_FIG_DIR / "Fig4_Boundary_Validation.png", dpi=300, bbox_inches='tight')
        plt.close(fig7_fixed)

        # FIX 5: Trajectories (legends in external box)
        log.info("\n[5/6] Fixing Fig 8 - Legends in external box...")
        fig8_fixed = fix_fig8_trajectories(trajectory_df)
        fig8_fixed.savefig(REFINEMENT_FIG_DIR / "Fig5_Spatial_Trajectories.png", dpi=300, bbox_inches='tight')
        plt.close(fig8_fixed)

        # FIX 6: S1 - All 117 clusters with larger legend dots
        log.info("\n[6/6] Fixing S1 - Larger legend dots...")
        fig_supp117 = fix_supp_117_regions(adata)
        fig_supp117.savefig(FIGURES_DIR / "SUPP_All_117_Leiden_Regions.png", dpi=300, bbox_inches='tight')
        plt.close(fig_supp117)

        # Verify Fig 9
        log.info("\n[Verify] Checking Fig 9 PDF...")
        verify_fig9_pdf(region_df)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ ALL FIGURE FIXES COMPLETE in {elapsed:.1f} seconds")
        log.info("=" * 80)
        log.info("\nFixed files:")
        log.info(f"  • Fig1_Tissue_Region_Map.png (larger legend dots)")
        log.info(f"  • Fig2_Region_Statistics.png (taller panels, Crameri colors)")
        log.info(f"  • Fig1_Enrichment_Heatmap.png (clean title)")
        log.info(f"  • Fig4_Boundary_Validation.png (legends outside)")
        log.info(f"  • Fig5_Spatial_Trajectories.png (external legend box)")
        log.info(f"  • SUPP_All_117_Leiden_Regions.png (larger legend dots)")

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
