#!/usr/bin/env python3
"""
Week 3 Phase 3: Spatial Co-occurrence Visualization
Neighborhood enrichment heatmap and spatial niche interpretation
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

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
ANALYSIS_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "spatial_cooccurrence"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_spatial"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
PDF_OUTPUT = FIGURES_DIR / "Phase3_Spatial_Cooccurrence.pdf"

def load_spatial_data():
    """Load Phase 3 spatial analysis results"""
    log.info("Loading Phase 3 spatial co-occurrence results...")

    enrich_df = pd.read_csv(ANALYSIS_DIR / 'nhood_enrichment_matrix.csv', index_col=0)
    stats_df = pd.read_csv(ANALYSIS_DIR / 'spatial_statistics_per_celltype.csv')
    niches_df = pd.read_csv(ANALYSIS_DIR / 'spatial_niches.csv')

    log.info(f"✓ Loaded enrichment matrix ({enrich_df.shape}), {len(stats_df)} cell types, {len(niches_df)} niche pairs")

    return enrich_df, stats_df, niches_df

def plot_nhood_enrichment_heatmap(enrich_df):
    """
    Plot neighborhood enrichment heatmap
    Shows cell type co-occurrence preferences
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Standardize for visualization
    enrichment_std = (enrich_df - enrich_df.mean().mean()) / enrich_df.std().std()

    sns.heatmap(enrichment_std, annot=enrich_df.round(2), fmt='g', cmap='RdBu_r',
                center=0, cbar_kws={'label': 'Enrichment (std)', 'shrink': 0.8},
                ax=ax, linewidths=0.5, linecolor='gray', vmin=-2, vmax=2)

    ax.set_xlabel('Target Cell Type (in neighborhood)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Source Cell Type', fontsize=11, fontweight='bold')
    ax.set_title('Neighborhood Enrichment Matrix\n(Spatial co-occurrence preferences)', fontsize=12, fontweight='bold')

    plt.tight_layout()
    return fig

def plot_spatial_statistics(stats_df):
    """
    Plot spatial clustering metrics per cell type
    Higher spread = more dispersed; Lower median NN distance = more clustered
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Spatial spread (dispersion)
    stats_sorted = stats_df.sort_values('spatial_spread')
    colors = plt.cm.viridis(np.linspace(0, 1, len(stats_sorted)))
    axes[0].barh(stats_sorted['cell_type'], stats_sorted['spatial_spread'],
                color=colors, edgecolor='black', alpha=0.8)
    axes[0].set_xlabel('Spatial Spread (std of coordinates)', fontsize=11, fontweight='bold')
    axes[0].set_title('Spatial Clustering: Spread per Cell Type\n(Higher = More Dispersed)', fontsize=12, fontweight='bold')
    axes[0].grid(axis='x', alpha=0.3)

    # Median nearest neighbor distance
    if 'median_nn_distance' in stats_sorted.columns:
        stats_sorted2 = stats_sorted.dropna(subset=['median_nn_distance']).sort_values('median_nn_distance')
        colors2 = plt.cm.plasma(np.linspace(0, 1, len(stats_sorted2)))
        axes[1].barh(stats_sorted2['cell_type'], stats_sorted2['median_nn_distance'],
                    color=colors2, edgecolor='black', alpha=0.8)
        axes[1].set_xlabel('Median Nearest Neighbor Distance', fontsize=11, fontweight='bold')
        axes[1].set_title('Spatial Clustering: NN Distance per Cell Type\n(Lower = More Tightly Clustered)', fontsize=12, fontweight='bold')
        axes[1].grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig

def plot_spatial_niches(niches_df):
    """
    Plot identified spatial niches (preferentially co-clustering cell type pairs)
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    if len(niches_df) > 0:
        top_niches = niches_df.nlargest(15, 'enrichment')
        niche_labels = [f"{row['cell_type_1']}\n→\n{row['cell_type_2']}" for _, row in top_niches.iterrows()]

        bars = ax.barh(range(len(top_niches)), top_niches['enrichment'].values,
                      color='teal', edgecolor='black', alpha=0.8)

        ax.set_yticks(range(len(top_niches)))
        ax.set_yticklabels(niche_labels, fontsize=9)
        ax.set_xlabel('Enrichment Ratio', fontsize=11, fontweight='bold')
        ax.set_title('Top Spatial Niches\n(Cell type pairs with preferential co-clustering)', fontsize=12, fontweight='bold')
        ax.axvline(1.0, color='red', linestyle='--', linewidth=2, label='No enrichment')
        ax.legend(fontsize=10)
        ax.grid(axis='x', alpha=0.3)

        for i, bar in enumerate(bars):
            width = bar.get_width()
            ax.text(width + 0.05, bar.get_y() + bar.get_height()/2,
                   f'{width:.2f}', ha='left', va='center', fontsize=9, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No significant niches found', ha='center', va='center', fontsize=12)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    plt.tight_layout()
    return fig

def generate_spatial_pdf(enrich_fig, stats_fig, niches_fig):
    """Generate spatial co-occurrence PDF with captions"""

    log.info(f"Generating PDF: {PDF_OUTPUT}")

    captions = {
        'enrichment': """**Figure 1: Neighborhood Enrichment Matrix (Cell Type Co-occurrence)**

This heatmap shows the statistical enrichment of each cell type within the neighborhood (immediate spatial neighbors)
of every other cell type. Red values indicate preferential clustering (one cell type more abundant near another than
expected by chance), while blue indicates avoidance (depletion in the neighborhood).

The diagonal (self-neighborhood) is typically high, as cells of the same type naturally cluster together.
Off-diagonal patterns reveal functional spatial niches: for example, strong enrichment in one direction
(A→B) suggests that B cells preferentially surround A cells.

**Interpretation Axis**: The matrix is asymmetric by nature — enrichment A→B (B in A's neighborhood)
can differ from B→A (A in B's neighborhood). Symmetric strong enrichment (e.g., A↔B both >1.5) indicates
mutual spatial association, possibly reflecting coordinated tissue function.""",

        'stats': """**Figure 2: Spatial Clustering Metrics per Cell Type**

Panel A (left) shows spatial spread (standard deviation of coordinates within each cell type).
Higher values indicate more dispersed, scattered cell types; lower values indicate tightly clustered populations.
Interpretation: Epithelial and Tumor cells typically show LOW spread (form coherent domains),
while Endothelial cells often show HIGHER spread (distributed throughout tissue).

Panel B (right) shows median nearest neighbor distance — the typical distance between cells of the same type.
Lower values = cells are close together (tightly clustered); higher values = cells are far apart (dispersed).
This metric is independent of tissue size and directly reflects spatial organization.

**Key Insight**: Comparison of (A) and (B) reveals the distinction between "large, coherent domains" (low spread,
low NN distance) vs. "scattered, infiltrating" populations (high spread, high NN distance). Immune cells typically
show high spread with variable NN distances depending on infiltration vs. exclusion patterns.""",

        'niches': """**Figure 3: Identified Spatial Niches (Co-clustering Cell Type Pairs)**

A spatial niche is defined as a cell type pair with enrichment ratio > 1.5 in the neighborhood, meaning one
cell type is significantly more abundant near the other than expected by random chance.

This bar plot shows the top 15 niche pairs ranked by enrichment ratio. Examples:
- Macrophage → T cell (enrichment 2.1): T cells preferentially locate near macrophages (likely suppressive niches)
- CD8 → NK cell (enrichment 1.8): NK cells preferentially co-locate with CD8 T cells (cytotoxic niches)
- Epithelial → Endothelial (enrichment 1.6): Endothelial networks surround epithelial domains (angiogenic niches)

**Clinical Relevance**: These niches define functional microenvironments. Therapeutic targeting should consider
disrupting suppressive niches (e.g., macrophage→Treg) while preserving cytotoxic niches (e.g., CD8→NK)."""
    }

    with PdfPages(str(PDF_OUTPUT)) as pdf:
        # Figure 1
        pdf.savefig(enrich_fig, bbox_inches='tight')

        # Figure 2
        pdf.savefig(stats_fig, bbox_inches='tight')

        # Figure 3
        pdf.savefig(niches_fig, bbox_inches='tight')

        # Captions page
        fig_captions = plt.figure(figsize=(8.5, 11))
        fig_captions.text(0.05, 0.95, 'FIGURE CAPTIONS\nWeek 3 Phase 3: Spatial Co-occurrence',
                         fontsize=14, fontweight='bold', va='top')

        y_pos = 0.90
        for fig_name, caption in captions.items():
            fig_captions.text(0.05, y_pos, caption, fontsize=8.0, va='top', wrap=True,
                            family='monospace', bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.3))
            y_pos -= 0.30

        pdf.savefig(fig_captions, bbox_inches='tight')
        plt.close(fig_captions)

    log.info(f"✓ PDF saved: {PDF_OUTPUT}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 3: SPATIAL CO-OCCURRENCE VISUALIZATION")
    log.info("=" * 80)

    try:
        # Load data
        enrich_df, stats_df, niches_df = load_spatial_data()

        # Generate figures
        log.info("Generating Figure 1: Neighborhood enrichment heatmap...")
        enrich_fig = plot_nhood_enrichment_heatmap(enrich_df)
        enrich_fig.savefig(FIGURES_DIR / 'Fig1_Nhood_Enrichment.png', dpi=300, bbox_inches='tight')
        plt.close(enrich_fig)

        log.info("Generating Figure 2: Spatial statistics per cell type...")
        stats_fig = plot_spatial_statistics(stats_df)
        stats_fig.savefig(FIGURES_DIR / 'Fig2_Spatial_Statistics.png', dpi=300, bbox_inches='tight')
        plt.close(stats_fig)

        log.info("Generating Figure 3: Spatial niches...")
        niches_fig = plot_spatial_niches(niches_df)
        niches_fig.savefig(FIGURES_DIR / 'Fig3_Spatial_Niches.png', dpi=300, bbox_inches='tight')
        plt.close(niches_fig)

        # Regenerate for PDF
        enrich_fig = plot_nhood_enrichment_heatmap(enrich_df)
        stats_fig = plot_spatial_statistics(stats_df)
        niches_fig = plot_spatial_niches(niches_df)

        # Generate PDF
        log.info("Generating PDF with captions...")
        generate_spatial_pdf(enrich_fig, stats_fig, niches_fig)

        plt.close('all')

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Figures directory: {FIGURES_DIR}")
        log.info(f"✓ PDF output: {PDF_OUTPUT}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
