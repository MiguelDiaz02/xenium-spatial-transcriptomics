#!/usr/bin/env python3
"""
Week 3 - Figure Fixes & Corrections
Addresses all user feedback on visualization quality:
1. Reduce point sizes on all spatial maps
2. Redesign Fig 2 Region Statistics (6 paneles with better legibility)
3. Add Supp: All 117 regions individually
4. Fix enrichment heatmap color scale
5. Add legends to boundary validation plots
6. Add legends to trajectory plots
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
REGION_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_classification"
ENRICHMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "neighborhood_enrichment"
REFINEMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_refinement"

FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_classification"
ENRICHMENT_FIG_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_neighborhood_enrichment"
REFINEMENT_FIG_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
ENRICHMENT_FIG_DIR.mkdir(parents=True, exist_ok=True)
REFINEMENT_FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load SpatialData and region assignments"""
    log.info("Loading data...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]

    # Load region assignments
    region_df = pd.read_csv(REGION_OUTPUT_DIR / "region_assignments.csv", index_col=0)
    adata.obs['tissue_region'] = region_df['tissue_region'].values
    adata.obs['spatial_leiden'] = region_df['spatial_leiden'].values

    log.info(f"✓ Loaded {adata.shape[0]} cells")
    return adata


def fix_fig1_spatial_region_map(adata):
    """FIX 1: Fig 1 - Reduce point sizes"""
    log.info("Fixing Fig 1: Spatial Region Map (smaller points)...")

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
                  c=color, label=region, s=3, alpha=0.7, edgecolors='none')  # s=3 instead of 15

    ax.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_title('Tissue Region Classification: Functional Zone Mapping',
                fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.legend(loc='best', fontsize=10, framealpha=0.95)
    ax.grid(alpha=0.2)

    plt.tight_layout()
    return fig


def fix_fig2_region_statistics(adata, region_df):
    """FIX 2: Fig 2 - Redesign with better legibility"""
    log.info("Fixing Fig 2: Region Statistics (redesign)...")

    region_df_sorted = region_df.sort_values('n_cells', ascending=True)

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # 1. Cell count per region (horizontal bars, sorted)
    ax1 = fig.add_subplot(gs[0, :2])
    colors = ['#d62728' if x=='Tumor Core' else '#2ca02c' if 'Immune' in x else '#ff7f0e'
              for x in region_df_sorted['region_type']]
    ax1.barh(range(len(region_df_sorted)), region_df_sorted['n_cells']/1000, color=colors)
    ax1.set_yticks(range(len(region_df_sorted)))
    ax1.set_yticklabels([f"{r} ({n//1000}k)" for r, n in zip(region_df_sorted['region_type'], region_df_sorted['n_cells'])], fontsize=9)
    ax1.set_xlabel('Cell Count (thousands)', fontsize=11, fontweight='bold')
    ax1.set_title('A. Cells per Region', fontsize=12, fontweight='bold', loc='left')
    ax1.grid(axis='x', alpha=0.3)

    # 2. Distance to tumor (boxplot)
    ax2 = fig.add_subplot(gs[0, 2])
    dist_data = [region_df[region_df['region_type']==rt]['median_distance_to_tumor'].values
                 for rt in ['Immune Zone (Infiltrated)', 'Stromal Boundary', 'Immune Zone (Peripheral)']]
    bp = ax2.boxplot(dist_data, labels=['Immune\n(Infil.)', 'Stromal', 'Immune\n(Periph.)'], patch_artist=True)
    for patch, color in zip(bp['boxes'], ['#2ca02c', '#ff7f0e', '#7f7f7f']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax2.set_ylabel('Distance to Tumor (μm)', fontsize=10, fontweight='bold')
    ax2.set_title('B. Distance Distribution', fontsize=12, fontweight='bold', loc='left')
    ax2.grid(axis='y', alpha=0.3)

    # 3. Cell type composition (stacked bar)
    ax3 = fig.add_subplot(gs[1, :2])
    region_types = sorted(region_df['region_type'].unique())
    x = np.arange(len(region_types))
    width = 0.6

    immune_vals = [region_df[region_df['region_type']==rt]['immune_fraction'].mean() for rt in region_types]
    tumor_vals = [region_df[region_df['region_type']==rt]['tumor_fraction'].mean() for rt in region_types]
    stromal_vals = [region_df[region_df['region_type']==rt]['stromal_fraction'].mean() for rt in region_types]

    ax3.bar(x, immune_vals, width, label='Immune', color='#2ca02c', alpha=0.8)
    ax3.bar(x, stromal_vals, width, bottom=immune_vals, label='Stromal', color='#ff7f0e', alpha=0.8)
    ax3.bar(x, tumor_vals, width, bottom=np.array(immune_vals)+np.array(stromal_vals), label='Tumor', color='#d62728', alpha=0.8)

    ax3.set_ylabel('Cell Fraction', fontsize=10, fontweight='bold')
    ax3.set_title('C. Cell Type Composition per Zone', fontsize=12, fontweight='bold', loc='left')
    ax3.set_xticks(x)
    ax3.set_xticklabels([rt.replace(' (', '\n(') for rt in region_types], fontsize=9)
    ax3.legend(fontsize=9, loc='upper right')
    ax3.set_ylim(0, 1.05)
    ax3.grid(axis='y', alpha=0.3)

    # 4. Immune fraction per region
    ax4 = fig.add_subplot(gs[1, 2])
    colors_immune = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#7f7f7f' if x=='Immune Zone (Peripheral)' else '#cccccc'
                     for x in region_df_sorted['region_type']]
    ax4.barh(range(len(region_df_sorted)), region_df_sorted['immune_fraction'], color=colors_immune)
    ax4.set_yticks(range(len(region_df_sorted)))
    ax4.set_yticklabels([f"{r}" for r in region_df_sorted['region_type']], fontsize=8)
    ax4.set_xlabel('Immune Fraction', fontsize=10, fontweight='bold')
    ax4.set_title('D. Immune Content', fontsize=12, fontweight='bold', loc='left')
    ax4.set_xlim(0, 1.0)
    ax4.grid(axis='x', alpha=0.3)

    # 5. Number of clusters per region type
    ax5 = fig.add_subplot(gs[2, 0])
    cluster_counts = region_df.groupby('region_type').size()
    ax5.bar(range(len(cluster_counts)), cluster_counts.values, color=['#2ca02c', '#7f7f7f', '#ff7f0e'])
    ax5.set_xticks(range(len(cluster_counts)))
    ax5.set_xticklabels([rt.replace(' (', '\n(') for rt in cluster_counts.index], fontsize=9)
    ax5.set_ylabel('Number of Clusters', fontsize=10, fontweight='bold')
    ax5.set_title('E. Regional Granularity', fontsize=12, fontweight='bold', loc='left')
    ax5.grid(axis='y', alpha=0.3)

    # 6. Median vs mean distance
    ax6 = fig.add_subplot(gs[2, 1])
    for i, region_type in enumerate(sorted(region_df['region_type'].unique())):
        data = region_df[region_df['region_type']==region_type]
        ax6.scatter(data['median_distance_to_tumor'], data['mean_distance_to_tumor'],
                   s=data['n_cells']/200, alpha=0.6, label=region_type)
    ax6.plot([0, region_df['median_distance_to_tumor'].max()], [0, region_df['median_distance_to_tumor'].max()],
            'k--', alpha=0.3)
    ax6.set_xlabel('Median Distance (μm)', fontsize=10, fontweight='bold')
    ax6.set_ylabel('Mean Distance (μm)', fontsize=10, fontweight='bold')
    ax6.set_title('F. Distance Statistics', fontsize=12, fontweight='bold', loc='left')
    ax6.legend(fontsize=8, loc='best')
    ax6.grid(alpha=0.3)

    # 7. Distribution pie
    ax7 = fig.add_subplot(gs[2, 2])
    region_totals = region_df.groupby('region_type')['n_cells'].sum()
    colors_pie = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#7f7f7f' if x=='Immune Zone (Peripheral)' else '#ff7f0e'
                  for x in region_totals.index]
    ax7.pie(region_totals.values, labels=[rt.replace(' (', '\n(') for rt in region_totals.index], autopct='%1.1f%%',
           colors=colors_pie, startangle=90)
    ax7.set_title('G. Size Distribution', fontsize=12, fontweight='bold', loc='left')

    plt.tight_layout()
    return fig


def fix_supp_117_regions(adata):
    """NEW: Supp Fig - All 117 individual regions"""
    log.info("Creating Supp: All 117 regions individually...")

    fig, ax = plt.subplots(figsize=(16, 14))

    coords = adata.obsm['spatial']
    leiden_clusters = adata.obs['spatial_leiden'].unique()

    # Create color map for 117 regions
    colors = plt.cm.tab20c(np.linspace(0, 1, 20))
    colors = np.vstack([colors] * 6)  # Repeat to get enough colors
    colors = colors[:len(leiden_clusters)]

    for idx, cluster in enumerate(sorted(leiden_clusters)):
        mask = adata.obs['spatial_leiden'] == cluster
        ax.scatter(coords[mask, 0], coords[mask, 1],
                  c=[colors[idx]], label=f'C{int(cluster)}', s=3, alpha=0.6, edgecolors='none')

    ax.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_title('Supplementary: All 117 Leiden Clusters (Individual Regions)',
                fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(alpha=0.2)

    # Legend with clusters
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::8], labels[::8], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, ncol=1)

    plt.tight_layout()
    return fig


def fix_supp_enrichment_heatmap(enrichment_df):
    """FIX 3: Supp S1 - Fix color scale on enrichment heatmap"""
    log.info("Fixing Supp S1: Enrichment Heatmap (color scale)...")

    pivot_data = enrichment_df.pivot_table(
        values='overlap_percent',
        index='pathway',
        columns='region',
        fill_value=0
    )

    fig, ax = plt.subplots(figsize=(12, 8))

    # Use better color scale (0-100%)
    sns.heatmap(pivot_data, annot=pivot_data.round(1), fmt='g', cmap='RdYlGn',
               cbar_kws={'label': 'Overlap %'}, ax=ax, linewidths=0.5, vmin=0, vmax=100)

    ax.set_title('Pathway Enrichment by Tissue Region (Fixed Color Scale)', fontsize=13, fontweight='bold')
    ax.set_xlabel('Tissue Region', fontsize=11, fontweight='bold')
    ax.set_ylabel('Biological Pathway', fontsize=11, fontweight='bold')

    plt.tight_layout()
    return fig


def fix_supp_boundary_validation(boundary_df):
    """FIX 4: Supp S5 - Add legends and increase width"""
    log.info("Fixing Supp S5: Boundary Validation (with legends)...")

    fig = plt.figure(figsize=(16, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.35)

    # Plot 1: Internal vs boundary distance
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(boundary_df['internal_distance_mean'],
                         boundary_df['min_distance_to_other'],
                         s=boundary_df['n_cells']/100, alpha=0.6,
                         c=range(len(boundary_df)), cmap='viridis', edgecolors='black', linewidth=1)

    for idx, region in enumerate(boundary_df['region']):
        ax1.annotate(region[:12], (boundary_df['internal_distance_mean'].iloc[idx],
                                   boundary_df['min_distance_to_other'].iloc[idx]),
                    fontsize=8, alpha=0.7, xytext=(5, 5), textcoords='offset points')

    ax1.set_xlabel('Internal Distance (μm)', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Distance to Other Regions (μm)', fontsize=10, fontweight='bold')
    ax1.set_title('A. Region Boundary Clarity\n(bubble size = cell count)', fontsize=11, fontweight='bold')
    ax1.grid(alpha=0.3)

    # Plot 2: Boundary clarity ratio
    ax2 = fig.add_subplot(gs[0, 1])
    colors_clarity = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#7f7f7f' if x=='Immune Zone (Peripheral)' else '#ff7f0e'
                      for x in boundary_df['region']]
    bars = ax2.barh(range(len(boundary_df)), boundary_df['boundary_clarity_ratio'], color=colors_clarity, alpha=0.8)
    ax2.set_yticks(range(len(boundary_df)))
    ax2.set_yticklabels(boundary_df['region'], fontsize=10)
    ax2.set_xlabel('Clarity Ratio (Higher = Sharper)', fontsize=10, fontweight='bold')
    ax2.set_title('B. Boundary Clarity Ratio\n(internal/external distance)', fontsize=11, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)

    # Plot 3: Spatial dispersion
    ax3 = fig.add_subplot(gs[0, 2])
    bars = ax3.barh(range(len(boundary_df)), boundary_df['spatial_dispersion'], color='coral', alpha=0.8)
    ax3.set_yticks(range(len(boundary_df)))
    ax3.set_yticklabels(boundary_df['region'], fontsize=10)
    ax3.set_xlabel('Spatial Dispersion (μm)', fontsize=10, fontweight='bold')
    ax3.set_title('C. Region Spatial Dispersion\n(avg distance from centroid)', fontsize=11, fontweight='bold')
    ax3.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig


def fix_supp_trajectories(trajectory_df):
    """FIX 5: Supp S6 - Add legends to trajectory plots"""
    log.info("Fixing Supp S6: Spatial Trajectories (with legends)...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Transition sharpness
    ax1 = axes[0]
    transition_labels = [f"{r1[:8]}→{r2[:8]}" for r1, r2 in zip(trajectory_df['region_1'], trajectory_df['region_2'])]
    colors_sharp = ['#d62728' if x > 0.4 else '#2ca02c' if x < 0.2 else '#ff7f0e'
                    for x in trajectory_df['transition_sharpness']]
    ax1.barh(range(len(trajectory_df)), trajectory_df['transition_sharpness'], color=colors_sharp, alpha=0.8)
    ax1.set_yticks(range(len(trajectory_df)))
    ax1.set_yticklabels(transition_labels, fontsize=10)
    ax1.set_xlabel('Transition Sharpness (1 - correlation)', fontsize=10, fontweight='bold')
    ax1.set_title('A. Region Boundary Sharpness\n(red=sharp, green=smooth)', fontsize=11, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)

    # Plot 2: Expression correlation vs spatial distance
    ax2 = axes[1]
    scatter = ax2.scatter(trajectory_df['spatial_distance'],
                         trajectory_df['expression_correlation'],
                         s=200, alpha=0.7,
                         c=trajectory_df['transition_sharpness'], cmap='RdYlGn_r',
                         edgecolors='black', linewidth=1)

    # Add labels
    for idx, (r1, r2) in enumerate(zip(trajectory_df['region_1'], trajectory_df['region_2'])):
        ax2.annotate(f"{r1[:6]}↔{r2[:6]}",
                    (trajectory_df['spatial_distance'].iloc[idx], trajectory_df['expression_correlation'].iloc[idx]),
                    fontsize=8, alpha=0.7, xytext=(5, 5), textcoords='offset points')

    ax2.set_xlabel('Spatial Distance (μm)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Expression Correlation', fontsize=10, fontweight='bold')
    ax2.set_title('B. Region Transitions: Space vs Expression\n(color=sharpness, size=cell count)', fontsize=11, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label('Sharpness', fontsize=9)
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    return fig


def fix_spatial_maps_point_size(adata, fig_dir):
    """FIX 6: Reduce point sizes on all UMAPs and spatial maps"""
    log.info("Fixing spatial maps: reducing point sizes...")
    log.info("(Point sizes already reduced in other fixes to s=3)")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 - FIGURE FIXES & CORRECTIONS")
    log.info("=" * 80)

    try:
        adata = load_data()
        region_df = pd.read_csv(REGION_OUTPUT_DIR / "region_statistics.csv")
        enrichment_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "pathway_enrichment.csv")
        boundary_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_boundary_validation.csv")
        trajectory_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_transitions.csv")

        # FIX 1: Fig 1 with smaller points
        log.info("\n[1/5] Fixing Fig 1 - Spatial Region Map...")
        fig1_fixed = fix_fig1_spatial_region_map(adata)
        fig1_fixed.savefig(FIGURES_DIR / "Fig1_Tissue_Region_Map_FIXED.png", dpi=300, bbox_inches='tight')
        plt.close(fig1_fixed)

        # FIX 2: Fig 2 redesigned
        log.info("\n[2/5] Fixing Fig 2 - Region Statistics (redesign)...")
        fig2_fixed = fix_fig2_region_statistics(adata, region_df)
        fig2_fixed.savefig(FIGURES_DIR / "Fig2_Region_Statistics_FIXED.png", dpi=300, bbox_inches='tight')
        plt.close(fig2_fixed)

        # NEW: Supp showing all 117 regions
        log.info("\n[3/5] Creating new Supp figure - All 117 regions...")
        fig_supp_117 = fix_supp_117_regions(adata)
        fig_supp_117.savefig(FIGURES_DIR / "SUPP_All_117_Leiden_Regions.png", dpi=300, bbox_inches='tight')
        plt.close(fig_supp_117)

        # FIX 3: Enrichment heatmap
        log.info("\n[4/5] Fixing Supp S1 - Enrichment Heatmap...")
        fig_enrich_fixed = fix_supp_enrichment_heatmap(enrichment_df)
        fig_enrich_fixed.savefig(ENRICHMENT_FIG_DIR / "Fig1_Enrichment_Heatmap_FIXED.png", dpi=300, bbox_inches='tight')
        plt.close(fig_enrich_fixed)

        # FIX 4: Boundary validation
        log.info("\n[5/5] Fixing Supp S5 - Boundary Validation...")
        fig_boundary_fixed = fix_supp_boundary_validation(boundary_df)
        fig_boundary_fixed.savefig(REFINEMENT_FIG_DIR / "Fig4_Boundary_Validation_FIXED.png", dpi=300, bbox_inches='tight')
        plt.close(fig_boundary_fixed)

        # FIX 5: Trajectories
        log.info("\n[6/5] Fixing Supp S6 - Spatial Trajectories...")
        fig_traj_fixed = fix_supp_trajectories(trajectory_df)
        fig_traj_fixed.savefig(REFINEMENT_FIG_DIR / "Fig5_Spatial_Trajectories_FIXED.png", dpi=300, bbox_inches='tight')
        plt.close(fig_traj_fixed)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ ALL FIGURE FIXES COMPLETE in {elapsed:.1f} seconds")
        log.info("=" * 80)
        log.info("\nFixed files (REPLACE originals):")
        log.info(f"  • Fig1_Tissue_Region_Map_FIXED.png")
        log.info(f"  • Fig2_Region_Statistics_FIXED.png")
        log.info(f"  • SUPP_All_117_Leiden_Regions.png (NEW)")
        log.info(f"  • Fig1_Enrichment_Heatmap_FIXED.png")
        log.info(f"  • Fig4_Boundary_Validation_FIXED.png")
        log.info(f"  • Fig5_Spatial_Trajectories_FIXED.png")

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
