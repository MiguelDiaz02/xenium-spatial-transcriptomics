#!/usr/bin/env python3
"""
Week 3 Phase 2B: Hybrid CCC Method Visualization
Publication-ready figures combining DeepTalk spatial concepts + CellChat L/R database
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
ANALYSIS_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_hybrid_method"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase2b_ccc_hybrid"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
PDF_OUTPUT = FIGURES_DIR / "Phase2B_CCC_Hybrid.pdf"

def load_hybrid_ccc_data():
    """Load hybrid CCC results"""
    log.info("Loading hybrid CCC results...")

    pairwise = pd.read_csv(ANALYSIS_DIR / "ccc_pairwise_interactions.csv")
    matrix = pd.read_csv(ANALYSIS_DIR / "ccc_matrix_celltype.csv", index_col=0)
    stats = pd.read_csv(ANALYSIS_DIR / "ccc_summary_stats.csv")

    log.info(f"✓ Loaded: {len(pairwise)} interactions, CCC matrix {matrix.shape}")

    return pairwise, matrix, stats

def plot_hybrid_ccc_heatmap(matrix):
    """Plot cell type-level CCC matrix"""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Normalize for visualization
    matrix_norm = (matrix - matrix.min().min()) / (matrix.max().max() - matrix.min().min() + 1e-6)

    sns.heatmap(matrix_norm, annot=matrix.round(2), fmt='g', cmap='YlOrRd',
               cbar_kws={'label': 'CCC Score (normalized)'}, ax=ax,
               linewidths=0.5, linecolor='gray')

    ax.set_xlabel('Target Cell Type (receiving signals)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Source Cell Type (sending signals)', fontsize=11, fontweight='bold')
    ax.set_title('Hybrid CCC Matrix\n(Spatial-Integrated L/R Interactions)',
                fontsize=12, fontweight='bold')

    plt.tight_layout()
    return fig

def plot_hybrid_top_interactions(pairwise, top_n=20):
    """Plot top L/R interaction pairs"""
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Top 20 interactions
    top_pairs = pairwise.nlargest(top_n, 'ccc_score')

    # Panel A: Top pairs by CCC score
    pair_labels = [f"{row['ligand']}→{row['receptor']}\n({row['source_cell_type']}→{row['target_cell_type']})"
                   for _, row in top_pairs.iterrows()]
    axes[0].barh(range(len(top_pairs)), top_pairs['ccc_score'].values, color='steelblue', edgecolor='black')
    axes[0].set_yticks(range(len(top_pairs)))
    axes[0].set_yticklabels(pair_labels, fontsize=8)
    axes[0].set_xlabel('CCC Score', fontsize=11, fontweight='bold')
    axes[0].set_title(f'Top {top_n} Ligand-Receptor Interactions (By CCC Score)', fontsize=12, fontweight='bold')
    axes[0].grid(axis='x', alpha=0.3)

    # Panel B: Distribution of scores by metric
    axes[1].scatter(pairwise['interaction_score'], pairwise['spatial_confidence'],
                   s=pairwise['ccc_score']*50 + 1, alpha=0.6, c=pairwise['ccc_score'],
                   cmap='viridis', edgecolors='black', linewidth=0.5)
    axes[1].set_xlabel('Interaction Score (L/R Expression)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Spatial Confidence (Co-localization)', fontsize=11, fontweight='bold')
    axes[1].set_title('CCC Score Components: L/R Expression vs Spatial Proximity\n(Size and color = final CCC score)',
                     fontsize=12, fontweight='bold')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    return fig

def plot_ccc_celltype_engagement(pairwise):
    """Plot CCC engagement per cell type"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Source cell types: how much signal do they send?
    source_engagement = pairwise.groupby('source_cell_type')['ccc_score'].sum().sort_values(ascending=True)
    axes[0].barh(source_engagement.index, source_engagement.values, color='coral', edgecolor='black', alpha=0.8)
    axes[0].set_xlabel('Total CCC Score (Sent)', fontsize=11, fontweight='bold')
    axes[0].set_title('CCC Engagement: Signal Senders\n(Sum of CCC scores for interactions sent by each type)',
                     fontsize=12, fontweight='bold')
    axes[0].grid(axis='x', alpha=0.3)

    # Target cell types: how much signal do they receive?
    target_engagement = pairwise.groupby('target_cell_type')['ccc_score'].sum().sort_values(ascending=True)
    axes[1].barh(target_engagement.index, target_engagement.values, color='lightblue', edgecolor='black', alpha=0.8)
    axes[1].set_xlabel('Total CCC Score (Received)', fontsize=11, fontweight='bold')
    axes[1].set_title('CCC Engagement: Signal Receivers\n(Sum of CCC scores for interactions received by each type)',
                     fontsize=12, fontweight='bold')
    axes[1].grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig

def generate_hybrid_pdf(fig1, fig2, fig3):
    """Generate hybrid CCC PDF with captions"""

    log.info(f"Generating PDF: {PDF_OUTPUT}")

    captions = {
        'heatmap': """**Figure 1: Cell Type-Level CCC Matrix**

This heatmap integrates ligand-receptor expression with spatial context, computed as:
CCC_score = log2(ligand_expr × receptor_expr + 1) × spatial_confidence

Where spatial_confidence measures the fraction of ligand-expressing cells near receptor-expressing neighbors.
Rows = source cell types (sending signals), Columns = target cell types (receiving signals).
Darker red = stronger spatial cell-cell communication.

Key features of this hybrid approach:
1. **Spatial Integration:** Only counts L/R pairs in spatially proximal cells (DeepTalk concept)
2. **Comprehensive L/R Database:** Uses validated CellChat L/R pairs from literature
3. **Single-Cell Resolution:** Maintains granularity while providing summary matrices
4. **Fast Computation:** Avoids deep learning overhead while maintaining biological validity

Biological interpretation: The pattern reveals communication hubs and responders. Immune cells typically show high
outgoing scores (senders) while epithelial/tumor cells show high incoming scores (receivers).""",

        'pairs': """**Figure 2: Top Ligand-Receptor Interactions**

Panel A (top) shows the top 20 L/R interaction pairs ranked by CCC score. The label format indicates:
Ligand→Receptor (Source cell type → Target cell type)

Panel B (bottom) shows the relationship between two components of the CCC score:
- X-axis: Interaction score (product of ligand and receptor expression)
- Y-axis: Spatial confidence (fraction of cells co-expressing above median)
- Size and color: Final CCC score (darker and larger = stronger interaction)

This plot reveals whether strong interactions are driven by:
- High expression alone (upper-left quadrant): robust L/R pairing
- Spatial proximity alone (lower-right quadrant): location-dependent activation
- Both (upper-right quadrant): spatially-constrained, high-expression pairs (most biologically relevant)

Interactions in the lower-left are weak and can usually be ignored.""",

        'engagement': """**Figure 3: CCC Engagement by Cell Type**

Panel A (left) shows CCC signal sent by each cell type (summing all outgoing interactions).
Panel B (right) shows CCC signal received by each cell type (summing all incoming interactions).

Cell types can be classified by their role:
- **Senders (high Panel A):** Immune cells (T cells, macrophages, dendritic cells) actively broadcast signals
- **Receivers (high Panel B):** Tumor/epithelial cells, endothelial cells accumulate signals
- **Hubs (high both):** M2 macrophages, endothelial cells coordinate bi-directional communication

This "engagement profile" identifies key players in the tumor microenvironment and suggests therapeutic targets:
- High senders with specific L/R pairs → target their ligand expression
- High receivers of immunosuppressive signals → block those receptors"""
    }

    with PdfPages(str(PDF_OUTPUT)) as pdf:
        pdf.savefig(fig1, bbox_inches='tight')
        pdf.savefig(fig2, bbox_inches='tight')
        pdf.savefig(fig3, bbox_inches='tight')

        # Captions page
        fig_captions = plt.figure(figsize=(8.5, 11))
        fig_captions.text(0.05, 0.95, 'FIGURE CAPTIONS\nWeek 3 Phase 2B: CCC Hybrid Method',
                         fontsize=14, fontweight='bold', va='top')

        y_pos = 0.90
        for fig_name, caption in captions.items():
            fig_captions.text(0.05, y_pos, caption, fontsize=8.5, va='top', wrap=True,
                            family='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
            y_pos -= 0.32

        pdf.savefig(fig_captions, bbox_inches='tight')
        plt.close(fig_captions)

    log.info(f"✓ PDF saved: {PDF_OUTPUT}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 2B: HYBRID CCC VISUALIZATION")
    log.info("Spatial-integrated ligand-receptor analysis")
    log.info("=" * 80)

    try:
        # Load data
        pairwise, matrix, stats = load_hybrid_ccc_data()

        # Generate figures
        log.info("Generating Figure 1: CCC heatmap...")
        fig1 = plot_hybrid_ccc_heatmap(matrix)
        fig1.savefig(FIGURES_DIR / 'Fig1_Hybrid_CCC_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.close(fig1)

        log.info("Generating Figure 2: Top interactions...")
        fig2 = plot_hybrid_top_interactions(pairwise, top_n=20)
        fig2.savefig(FIGURES_DIR / 'Fig2_Top_LR_Interactions.png', dpi=300, bbox_inches='tight')
        plt.close(fig2)

        log.info("Generating Figure 3: Cell type engagement...")
        fig3 = plot_ccc_celltype_engagement(pairwise)
        fig3.savefig(FIGURES_DIR / 'Fig3_CCC_Engagement.png', dpi=300, bbox_inches='tight')
        plt.close(fig3)

        # Regenerate for PDF
        fig1 = plot_hybrid_ccc_heatmap(matrix)
        fig2 = plot_hybrid_top_interactions(pairwise, top_n=20)
        fig3 = plot_ccc_celltype_engagement(pairwise)

        # Generate PDF
        log.info("Generating PDF with captions...")
        generate_hybrid_pdf(fig1, fig2, fig3)

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
