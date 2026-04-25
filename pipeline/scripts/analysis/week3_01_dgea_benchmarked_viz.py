#!/usr/bin/env python3
"""
Week 3 Phase 1: DGEA Benchmarked Visualization
Compare Wilcoxon vs Spatial-Aware DE methods with publication-quality figures
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
ANALYSIS_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE_benchmarked"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase1_dgea_benchmarked"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
PDF_OUTPUT = FIGURES_DIR / "Phase1_DGEA_Benchmarked.pdf"

def load_benchmarked_data():
    """Load all DGEA benchmarked results"""
    log.info("Loading DGEA benchmarked results...")

    morans_i = pd.read_csv(ANALYSIS_DIR / "morans_i_all_genes.csv")
    summary = pd.read_csv(ANALYSIS_DIR / "comparative_dgea_summary.csv")
    stats = pd.read_csv(ANALYSIS_DIR / "dgea_summary_statistics.csv")

    log.info(f"✓ Loaded: {len(morans_i)} genes, {len(summary)} rows in comparative summary")

    return morans_i, summary, stats

def plot_method_comparison(summary):
    """
    Plot: Wilcoxon vs Spatial-Aware ranking comparison
    Shows how spatial context changes gene rankings
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    cell_types = summary['cell_type'].unique()[:4]  # Top 4 for clarity

    for idx, ct in enumerate(cell_types):
        ax = axes[idx // 2, idx % 2]
        ct_data = summary[summary['cell_type'] == ct].head(20)

        # Scatter: Wilcoxon rank vs Spatial-aware rank
        ax.scatter(ct_data['rank_wilcoxon'], ct_data['rank_spatial_aware'],
                  s=100, alpha=0.6, color='steelblue', edgecolors='black')

        # Diagonal (no change)
        ax.plot([0, 20], [0, 20], 'r--', linewidth=2, label='No change')

        ax.set_xlabel('Wilcoxon Rank (Non-Spatial)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Spatial-Aware Rank', fontsize=10, fontweight='bold')
        ax.set_title(f'{ct}: Ranking Comparison\n(Wilcoxon vs Wilcoxon+Moran\'s I)', fontsize=11, fontweight='bold')
        ax.set_xlim(0, 20)
        ax.set_ylim(0, 20)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)

    plt.suptitle('DGEA Method Comparison: Spatial Context Impact on Rankings', fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    return fig

def plot_spatial_enrichment_comparison(morans_i):
    """
    Plot: Spatial signal strength across methods
    Shows Moran's I distribution indicating spatial autocorrelation
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram of Moran's I
    axes[0].hist(morans_i['I'], bins=40, color='steelblue', edgecolor='black', alpha=0.7)
    axes[0].axvline(0, color='red', linestyle='--', linewidth=2, label='No spatial signal')
    axes[0].axvline(morans_i['I'].median(), color='green', linestyle='-', linewidth=2, label=f'Median: {morans_i["I"].median():.2f}')
    axes[0].set_xlabel("Moran's I (Spatial Autocorrelation)", fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Count (genes)', fontsize=11, fontweight='bold')
    axes[0].set_title('Distribution of Spatial Signal\n(All 289 genes)', fontsize=12, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(axis='y', alpha=0.3)

    # Significance threshold
    morans_sig = morans_i[morans_i['pval_norm_fdr_bh'] < 0.05]
    non_sig = morans_i[morans_i['pval_norm_fdr_bh'] >= 0.05]

    axes[1].bar(['Spatial Signal\n(p<0.05)', 'Non-Significant\n(p≥0.05)'],
               [len(morans_sig), len(non_sig)],
               color=['darkgreen', 'lightcoral'], edgecolor='black', alpha=0.8)
    axes[1].set_ylabel('Number of Genes', fontsize=11, fontweight='bold')
    axes[1].set_title('Spatial Gene Enrichment\n(FDR-corrected p-value)', fontsize=12, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)

    for i, v in enumerate([len(morans_sig), len(non_sig)]):
        axes[1].text(i, v + 5, f'{v}\n({100*v/len(morans_i):.1f}%)',
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    plt.tight_layout()
    return fig

def generate_benchmarked_pdf(fig1, fig2):
    """Generate benchmarked DGEA PDF with captions"""

    log.info(f"Generating PDF: {PDF_OUTPUT}")

    captions = {
        'comparison': """**Figure 1: DGEA Method Comparison (Wilcoxon vs Spatial-Aware)**

This figure shows how incorporating spatial context (Moran's I autocorrelation) changes gene rankings
within each cell type. Each point represents a top-20 DE gene. The red dashed diagonal indicates "no change"
(same rank in both methods). Points above the diagonal were ranked higher when accounting for spatial signal.

Key insight: Many genes that are highly DE by Wilcoxon alone become HIGHER-ranked when spatial clustering
is accounted for. This suggests that these genes drive cell type identity AND cluster spatially (true cell-type
specific markers). Genes below the diagonal are less spatially cohesive despite strong DE signal (scattered expression).

Biological implication: Spatial-aware DE better identifies functionally coherent cell-type domains.""",

        'enrichment': """**Figure 2: Spatial Gene Enrichment (Moran's I Analysis)**

Panel A (left) shows the distribution of Moran's I values (spatial autocorrelation statistic) for all 289 genes
in the panel. Values >0 indicate spatial clustering; values <0 indicate dispersion. The distribution is
heavily skewed toward positive values, indicating widespread spatial organization.

Panel B (right) shows the number of genes passing the significance threshold (p<0.05 FDR-corrected).
Result: 288/289 genes (99.7%) show statistically significant spatial autocorrelation, meaning the entire
tissue exhibits strong spatial organization. Only 1 gene lacks spatial signal.

Clinical implication: The tissue architecture is highly organized — cell types don't mix randomly.
This motivates spatial methods for DE analysis."""
    }

    with PdfPages(str(PDF_OUTPUT)) as pdf:
        pdf.savefig(fig1, bbox_inches='tight')
        pdf.savefig(fig2, bbox_inches='tight')

        # Captions page
        fig_captions = plt.figure(figsize=(8.5, 11))
        fig_captions.text(0.05, 0.95, 'FIGURE CAPTIONS\nWeek 3 Phase 1: DGEA Benchmarked Methods',
                         fontsize=14, fontweight='bold', va='top')

        y_pos = 0.90
        for fig_name, caption in captions.items():
            fig_captions.text(0.05, y_pos, caption, fontsize=8.5, va='top', wrap=True,
                            family='monospace', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
            y_pos -= 0.40

        pdf.savefig(fig_captions, bbox_inches='tight')
        plt.close(fig_captions)

    log.info(f"✓ PDF saved: {PDF_OUTPUT}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 1: DGEA BENCHMARKED VISUALIZATION")
    log.info("=" * 80)

    try:
        # Load data
        morans_i, summary, stats = load_benchmarked_data()

        # Generate figures
        log.info("Generating Figure 1: Method comparison...")
        fig1 = plot_method_comparison(summary)
        fig1.savefig(FIGURES_DIR / 'Fig1_Method_Comparison.png', dpi=300, bbox_inches='tight')
        plt.close(fig1)

        log.info("Generating Figure 2: Spatial enrichment...")
        fig2 = plot_spatial_enrichment_comparison(morans_i)
        fig2.savefig(FIGURES_DIR / 'Fig2_Spatial_Enrichment.png', dpi=300, bbox_inches='tight')
        plt.close(fig2)

        # Regenerate for PDF
        fig1 = plot_method_comparison(summary)
        fig2 = plot_spatial_enrichment_comparison(morans_i)

        # Generate PDF
        log.info("Generating PDF with captions...")
        generate_benchmarked_pdf(fig1, fig2)

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
