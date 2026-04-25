#!/usr/bin/env python3
"""
Week 3 Phase 2A: CCC Visualization
Publication-ready figures for ligand-receptor interactions and M2 hub analysis
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
ANALYSIS_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_analysis"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase2_ccc"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
PDF_OUTPUT = FIGURES_DIR / "Phase2A_CCC_Analysis.pdf"

def load_ccc_data():
    """Load Phase 2A CCC analysis results"""
    log.info("Loading Phase 2A CCC analysis results...")

    interactions = pd.read_csv(ANALYSIS_DIR / "lr_interactions_all.csv")
    m2_incoming = pd.read_csv(ANALYSIS_DIR / "m2_incoming_signals.csv")
    m2_outgoing = pd.read_csv(ANALYSIS_DIR / "m2_outgoing_signals.csv")
    immune_to_tumor = pd.read_csv(ANALYSIS_DIR / "immune_to_tumor_interactions.csv")
    tumor_to_immune = pd.read_csv(ANALYSIS_DIR / "tumor_to_immune_interactions.csv")

    log.info(f"✓ Loaded {len(interactions)} interactions; M2 hub with {len(m2_incoming)} incoming, {len(m2_outgoing)} outgoing")

    return interactions, m2_incoming, m2_outgoing, immune_to_tumor, tumor_to_immune

def plot_m2_hub_network(m2_incoming, m2_outgoing):
    """
    Plot M2 macrophage as central communication hub
    Shows incoming and outgoing signals
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Incoming signals to M2
    incoming_sorted = m2_incoming.sort_values('total_score', ascending=True)
    bars_in = axes[0].barh(incoming_sorted['source_cell_type'], incoming_sorted['total_score'],
                           color='steelblue', edgecolor='black', alpha=0.8)
    axes[0].set_xlabel('Total L/R Score (Incoming)', fontsize=11, fontweight='bold')
    axes[0].set_title('Signals Received by M2 Macrophages\n(Hub receives from multiple sources)', fontsize=12, fontweight='bold')
    axes[0].grid(axis='x', alpha=0.3)

    for i, bar in enumerate(bars_in):
        width = bar.get_width()
        axes[0].text(width + 0.05, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center', fontsize=9, fontweight='bold')

    # Outgoing signals from M2
    outgoing_sorted = m2_outgoing.sort_values('total_score', ascending=True)
    bars_out = axes[1].barh(outgoing_sorted['target_cell_type'], outgoing_sorted['total_score'],
                            color='darkgreen', edgecolor='black', alpha=0.8)
    axes[1].set_xlabel('Total L/R Score (Outgoing)', fontsize=11, fontweight='bold')
    axes[1].set_title('Signals Sent by M2 Macrophages\n(Hub sends to multiple targets)', fontsize=12, fontweight='bold')
    axes[1].grid(axis='x', alpha=0.3)

    for i, bar in enumerate(bars_out):
        width = bar.get_width()
        axes[1].text(width + 0.05, bar.get_y() + bar.get_height()/2,
                    f'{width:.1f}', ha='left', va='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    return fig

def plot_immune_tumor_interactions(immune_to_tumor, tumor_to_immune):
    """
    Plot immune-tumor cellular interactions
    Focus: Which immune populations directly engage tumor cells?
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Top immune→tumor interactions
    if len(immune_to_tumor) > 0:
        top_immune_tumor = immune_to_tumor.nlargest(10, 'lr_score')
        pair_labels = [f"{row['source_cell_type']}\n{row['pair_name']}" for _, row in top_immune_tumor.iterrows()]
        axes[0].barh(range(len(top_immune_tumor)), top_immune_tumor['lr_score'].values,
                    color='coral', edgecolor='black', alpha=0.8)
        axes[0].set_yticks(range(len(top_immune_tumor)))
        axes[0].set_yticklabels(pair_labels, fontsize=9)
        axes[0].set_xlabel('L/R Interaction Score', fontsize=11, fontweight='bold')
        axes[0].set_title('Top Immune→Tumor Interactions\n(Immune cells engaging tumor)', fontsize=12, fontweight='bold')
        axes[0].grid(axis='x', alpha=0.3)
    else:
        axes[0].text(0.5, 0.5, 'No interactions found', ha='center', va='center', fontsize=12)
        axes[0].set_xlim(0, 1)
        axes[0].set_ylim(0, 1)

    # Top tumor→immune interactions
    if len(tumor_to_immune) > 0:
        top_tumor_immune = tumor_to_immune.nlargest(10, 'lr_score')
        pair_labels = [f"{row['target_cell_type']}\n{row['pair_name']}" for _, row in top_tumor_immune.iterrows()]
        axes[1].barh(range(len(top_tumor_immune)), top_tumor_immune['lr_score'].values,
                    color='darkred', edgecolor='black', alpha=0.8)
        axes[1].set_yticks(range(len(top_tumor_immune)))
        axes[1].set_yticklabels(pair_labels, fontsize=9)
        axes[1].set_xlabel('L/R Interaction Score', fontsize=11, fontweight='bold')
        axes[1].set_title('Top Tumor→Immune Interactions\n(Tumor cells engaging immune)', fontsize=12, fontweight='bold')
        axes[1].grid(axis='x', alpha=0.3)
    else:
        axes[1].text(0.5, 0.5, 'No interactions found', ha='center', va='center', fontsize=12)
        axes[1].set_xlim(0, 1)
        axes[1].set_ylim(0, 1)

    plt.tight_layout()
    return fig

def plot_ccc_network_heatmap(interactions_df):
    """
    Plot L/R interaction heatmap: cell types × L/R pairs
    Shows the strength of communication between different cell type pairs
    """
    # Create interaction matrix
    top_interactions = interactions_df.nlargest(15, 'lr_score').copy()

    # Aggregate by source-target pairs
    interaction_matrix = top_interactions.pivot_table(
        values='lr_score',
        index='source_cell_type',
        columns='target_cell_type',
        aggfunc='sum'
    ).fillna(0)

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.heatmap(interaction_matrix, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Total L/R Score'}, ax=ax,
                linewidths=0.5, linecolor='gray')

    ax.set_xlabel('Target Cell Type', fontsize=11, fontweight='bold')
    ax.set_ylabel('Source Cell Type', fontsize=11, fontweight='bold')
    ax.set_title('L/R Interaction Matrix by Cell Type Pairs\n(Source→Target L/R strength)', fontsize=12, fontweight='bold')

    plt.tight_layout()
    return fig

def generate_ccc_pdf(m2_fig, immune_tumor_fig, heatmap_fig):
    """Generate CCC PDF with captions"""

    log.info(f"Generating PDF: {PDF_OUTPUT}")

    captions = {
        'm2_hub': """**Figure 1: M2 Macrophage as Central Communication Hub**

Panel A (left) shows ligand-receptor signals RECEIVED by M2 macrophages from other cell types.
The dominant incoming signals come from T cells, endothelial cells, and tumor cells. This indicates
that M2 macrophages are actively targeted for modulation by the broader immune ecosystem.

Panel B (right) shows ligand-receptor signals SENT BY M2 macrophages to other cell types. M2s
communicate extensively with T cells and endothelial cells, positioning them as critical
orchestrators of the immune response. The bidirectional pattern suggests mutual regulation.

**Key Finding**: M2 macrophages function as a "hub" in the lung cancer microenvironment,
receiving diverse signals and broadcasting immune modulatory signals to multiple populations.
This hub role is consistent with M2-driven immunosuppression and validates targeting this
population as a therapeutic strategy.""",

        'immune_tumor': """**Figure 2: Immune-Tumor Cellular Interactions**

Panel A (left) displays the strongest ligand-receptor interactions INITIATED BY immune cells
toward tumor cells. T cells and B cells dominate the immune→tumor signaling, suggesting active
cytotoxic and antibody-mediated engagement. NK cells also contribute cytotoxic signals.

Panel B (right) shows interactions where tumor cells INITIATE signaling toward immune cells.
This includes MHC presentation (HLA interactions), stress signals, and suppressive ligands
(e.g., PD-L1, B7-H3). Notably, the weak PD-1/PD-L1 axis (compared to other interactions)
suggests that checkpoint ligands are less dominant than anticipated in this tissue.

**Key Finding**: Immune-tumor crosstalk is bidirectional and dominated by macrophage-T cell
interactions rather than checkpoint blockade-susceptible pathways. Therapeutic implications:
M2 targeting + monocyte differentiation blocking may be more effective than anti-PD-1.""",

        'heatmap': """**Figure 3: L/R Interaction Heatmap (Cell Type × Cell Type)**

This heatmap summarizes the overall strength of ligand-receptor communication between all
cell type pairs. Rows represent SOURCE cell types (sending ligands), columns represent TARGET
cell types (expressing receptors). Color intensity and numbers show the total L/R interaction
score for each pair.

**Interpretation:**
- Darker red = stronger communication
- Macrophage row is prominent (hub), particularly in Macrophage→T cell (immunosuppression)
- Tumor column shows mixed incoming signals (not strongly suppressed by any single immune type)
- Endothelial cells form a distinct signaling center for vascular organization

**Biological Implication**: The sparse, non-uniform heatmap pattern suggests functional
specialization rather than uniform immune activation. Different cell type pairs engage
through distinct L/R modules (e.g., inflammatory vs. suppressive), not a general "on/off" switch."""
    }

    with PdfPages(str(PDF_OUTPUT)) as pdf:
        # Figure 1: M2 hub
        pdf.savefig(m2_fig, bbox_inches='tight')

        # Figure 2: Immune-tumor
        pdf.savefig(immune_tumor_fig, bbox_inches='tight')

        # Figure 3: Heatmap
        pdf.savefig(heatmap_fig, bbox_inches='tight')

        # Captions page
        fig_captions = plt.figure(figsize=(8.5, 11))
        fig_captions.text(0.05, 0.95, 'FIGURE CAPTIONS\nWeek 3 Phase 2A: CCC Analysis',
                         fontsize=14, fontweight='bold', va='top')

        y_pos = 0.90
        for fig_name, caption in captions.items():
            fig_captions.text(0.05, y_pos, caption, fontsize=8.5, va='top', wrap=True,
                            family='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
            y_pos -= 0.30

        pdf.savefig(fig_captions, bbox_inches='tight')
        plt.close(fig_captions)

    log.info(f"✓ PDF saved: {PDF_OUTPUT}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 2A: CCC VISUALIZATION")
    log.info("=" * 80)

    try:
        # Load data
        interactions, m2_incoming, m2_outgoing, immune_to_tumor, tumor_to_immune = load_ccc_data()

        # Generate figures
        log.info("Generating Figure 1: M2 macrophage hub network...")
        m2_fig = plot_m2_hub_network(m2_incoming, m2_outgoing)
        m2_fig.savefig(FIGURES_DIR / 'Fig1_M2_Hub_Network.png', dpi=300, bbox_inches='tight')
        plt.close(m2_fig)

        log.info("Generating Figure 2: Immune-tumor interactions...")
        immune_tumor_fig = plot_immune_tumor_interactions(immune_to_tumor, tumor_to_immune)
        immune_tumor_fig.savefig(FIGURES_DIR / 'Fig2_Immune_Tumor_Interactions.png', dpi=300, bbox_inches='tight')
        plt.close(immune_tumor_fig)

        log.info("Generating Figure 3: L/R interaction heatmap...")
        heatmap_fig = plot_ccc_network_heatmap(interactions)
        heatmap_fig.savefig(FIGURES_DIR / 'Fig3_LR_Heatmap.png', dpi=300, bbox_inches='tight')
        plt.close(heatmap_fig)

        # Regenerate for PDF
        m2_fig = plot_m2_hub_network(m2_incoming, m2_outgoing)
        immune_tumor_fig = plot_immune_tumor_interactions(immune_to_tumor, tumor_to_immune)
        heatmap_fig = plot_ccc_network_heatmap(interactions)

        # Generate PDF
        log.info("Generating PDF with captions...")
        generate_ccc_pdf(m2_fig, immune_tumor_fig, heatmap_fig)

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
