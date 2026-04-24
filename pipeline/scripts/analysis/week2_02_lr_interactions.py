"""
Week 2 Task 2.2: Ligand-Receptor Interactions — Immune × Tumor Focus
=====================================================================
Analyzes L/R interactions between immune subtypes and tumor/epithelial cells.
Uses curated L/R pairs from literature that overlap with the 289-gene panel.

Input:  human_lung_cancer/results/sdata.zarr
Output: results/02_biology/lr_immune_tumor/ — interaction tables, heatmaps, bubble plots
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
import spatialdata as sd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "lr_immune_tumor"

# Curated L/R pairs: immune checkpoint, adhesion, cytokine signaling
# Sourced from CellChatDB, NicheNet, and lung cancer literature
LR_PAIRS = [
    # T cell checkpoints
    ('PDCD1', 'CD274'),  # PD-1 / PD-L1
    ('PDCD1', 'PDCD1LG2'),  # PD-1 / PD-L2
    ('CTLA4', 'CD80'),  # CTLA-4 / B7-1
    ('CTLA4', 'CD86'),  # CTLA-4 / B7-2
    # T cell activation
    ('CD3D', 'CD28'),  # TCR / costimulation (CD28 as marker)
    ('CD8A', 'HLA-A'),  # CD8 / MHC-I
    # T-Treg interaction
    ('FOXP3', 'IL2RA'),  # Treg markers
    # Macrophage-T cell
    ('CD68', 'CD3D'),  # Mac / T cell markers
    ('TNF', 'CD3D'),  # Inflammatory signal
    # Macrophage polarization
    ('CD68', 'CD163'),  # M1/M2 markers
    ('CD68', 'ARG1'),  # M2 polarization
    # B cell maturation
    ('MS4A1', 'CD40'),  # B cell / activation
    ('MS4A1', 'CD86'),  # B cell / APC interaction
    # Adhesion molecules (tumor-immune interface)
    ('ICAM1', 'ITGAL'),  # ICAM-1 / LFA-1
    ('PECAM1', 'PECAM1'),  # PECAM-1 (homophilic)
    # Immune homing
    ('CCR7', 'CCL19'),  # DC / T cell homing
    ('CCR2', 'CCL2'),  # Macrophage recruitment
    # Cytokines (if present in panel)
    ('IFNG', 'IFNGR1'),  # IFN-γ signaling
    ('TNF', 'TNFRSF1A'),  # TNF signaling
    ('IL2', 'IL2RA'),  # IL-2 signaling
    # Tumor-immune
    ('EPCAM', 'CD3D'),  # Epithelial-T cell proximity
]

def load_data():
    """Load SpatialData and extract AnnData."""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


def filter_lr_pairs(adata):
    """Keep only L/R pairs where both genes exist in the panel."""
    available_genes = set(adata.var_names)
    valid_pairs = []

    for ligand, receptor in LR_PAIRS:
        if ligand in available_genes and receptor in available_genes:
            valid_pairs.append((ligand, receptor))
        else:
            missing = [g for g in [ligand, receptor] if g not in available_genes]
            log.debug(f"Skipping {ligand}-{receptor} (missing: {missing})")

    log.info(f"Found {len(valid_pairs)} / {len(LR_PAIRS)} L/R pairs in panel")
    return valid_pairs


def compute_lr_scores(adata, lr_pairs):
    """Compute L/R interaction scores: mean expression per cell type pair."""
    log.info("Computing L/R expression per cell type...")

    # Use broad cell types + granular immune subtypes
    adata.obs['cell_type_combined'] = adata.obs['cell_type'].astype(str)
    immune_mask = adata.obs['immune_confidence'] != 'N/A'
    adata.obs.loc[immune_mask, 'cell_type_combined'] = adata.obs.loc[immune_mask, 'cell_type_immune_granular']

    # Get X as dense (normalized, log1p)
    X = adata.X.toarray() if sparse.issparse(adata.X) else adata.X

    # Compute mean expression per cell type
    cell_types = sorted(adata.obs['cell_type_combined'].unique())
    expr_per_type = {}

    for ct in cell_types:
        mask = adata.obs['cell_type_combined'] == ct
        if mask.sum() > 0:
            expr_per_type[ct] = X[mask].mean(axis=0)

    # Create L/R interaction matrix
    interaction_scores = []

    for ligand, receptor in lr_pairs:
        ligand_idx = adata.var_names.get_loc(ligand)
        receptor_idx = adata.var_names.get_loc(receptor)

        for sender in cell_types:
            for receiver in cell_types:
                if sender == receiver:
                    continue
                l_expr = expr_per_type[sender][ligand_idx]
                r_expr = expr_per_type[receiver][receptor_idx]

                interaction_scores.append({
                    'ligand': ligand,
                    'receptor': receptor,
                    'sender': sender,
                    'receiver': receiver,
                    'ligand_expr': l_expr,
                    'receptor_expr': r_expr,
                    'interaction_score': l_expr * r_expr  # Product
                })

    df_lr = pd.DataFrame(interaction_scores)
    log.info(f"Computed scores for {len(df_lr)} cell type pairs × {len(lr_pairs)} L/R pairs")
    return df_lr


def export_lr_tables(df_lr, output_dir):
    """Export L/R interaction tables."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Exporting L/R tables to {output_dir}...")

    # Full table
    df_lr.to_csv(output_dir / "LR_all_interactions.csv", index=False)

    # Top interactions
    df_top = df_lr.nlargest(50, 'interaction_score')
    df_top.to_csv(output_dir / "LR_top50_interactions.csv", index=False)
    log.info(f"  Full table: {len(df_lr)} rows")
    log.info(f"  Top 50: {df_top.shape[0]} rows")
    log.info(f"    Top interaction: {df_top.iloc[0]['sender']} → {df_top.iloc[0]['receiver']}: "
             f"{df_top.iloc[0]['ligand']}-{df_top.iloc[0]['receptor']} (score={df_top.iloc[0]['interaction_score']:.3f})")

    return df_lr


def plot_lr_heatmap(df_lr, output_dir):
    """Create heatmap of top L/R interactions."""
    log.info("Creating L/R interaction heatmap...")

    # Get top interactions by L/R pair
    df_top_pairs = df_lr.nlargest(30, 'interaction_score')

    # Create interaction label
    df_top_pairs['LR_pair'] = df_top_pairs['ligand'] + '-' + df_top_pairs['receptor']
    df_top_pairs['sender_receiver'] = df_top_pairs['sender'] + ' → ' + df_top_pairs['receiver']

    # Pivot to matrix
    pivot = df_top_pairs.pivot_table(
        index='LR_pair',
        columns='sender_receiver',
        values='interaction_score',
        fill_value=0
    )

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(pivot, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Interaction Score'})
    ax.set_title('Top 30 Ligand-Receptor Interactions (by cell type pair)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / "LR_heatmap.pdf", dpi=100, bbox_inches='tight')
    plt.close()
    log.info("Heatmap saved")


def plot_lr_bubble(df_lr, output_dir):
    """Create bubble plot of top L/R pairs."""
    log.info("Creating L/R bubble plot...")

    df_top = df_lr.nlargest(20, 'interaction_score').copy()
    df_top['LR_pair'] = df_top['ligand'] + '-' + df_top['receptor']
    df_top['cell_pair'] = df_top['sender'] + '→' + df_top['receiver']

    fig, ax = plt.subplots(figsize=(12, 8))
    scatter = ax.scatter(
        range(len(df_top)),
        range(len(df_top)),
        s=df_top['interaction_score'] * 500,
        c=df_top['interaction_score'],
        cmap='YlOrRd',
        alpha=0.6,
        edgecolors='black',
        linewidth=1
    )

    ax.set_xticks(range(len(df_top)))
    ax.set_xticklabels(df_top['LR_pair'], rotation=45, ha='right')
    ax.set_yticks(range(len(df_top)))
    ax.set_yticklabels(df_top['cell_pair'], fontsize=9)
    ax.set_xlabel('L/R Pair')
    ax.set_ylabel('Sender → Receiver')
    ax.set_title('Top 20 L/R Interactions by Expression', fontsize=12, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Interaction Score')
    plt.tight_layout()
    plt.savefig(output_dir / "LR_bubble.pdf", dpi=100, bbox_inches='tight')
    plt.close()
    log.info("Bubble plot saved")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 2 TASK 2.2: Ligand-Receptor Interactions")
    log.info("=" * 80)

    # Load data
    adata = load_data()

    # Filter L/R pairs to those in panel
    lr_pairs = filter_lr_pairs(adata)

    # Compute interaction scores
    df_lr = compute_lr_scores(adata, lr_pairs)

    # Export and visualize
    df_lr = export_lr_tables(df_lr, OUTPUT_DIR)
    plot_lr_heatmap(df_lr, OUTPUT_DIR)
    plot_lr_bubble(df_lr, OUTPUT_DIR)

    # Summary
    elapsed = time.time() - start_time
    log.info("=" * 80)
    log.info(f"TASK 2.2 COMPLETE in {elapsed:.1f} seconds")
    log.info(f"Outputs: {OUTPUT_DIR}")
    log.info(f"  - LR_all_interactions.csv")
    log.info(f"  - LR_top50_interactions.csv")
    log.info(f"  - LR_heatmap.pdf")
    log.info(f"  - LR_bubble.pdf")
    log.info("=" * 80)


if __name__ == "__main__":
    main()
