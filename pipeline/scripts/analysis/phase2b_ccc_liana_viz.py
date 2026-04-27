#!/usr/bin/env python3
"""
Phase 2B (LIANA+): Visualizations
Publication-ready figures from LIANA+ rank_aggregate output.
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

INPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_liana"
FIG_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase2b_ccc_liana"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def fig_top_interactions(sig_df, n_top=25):
    """Horizontal bar plot of top interactions by magnitude rank."""
    log.info(f"Plotting top {n_top} interactions...")
    df = sig_df.sort_values('magnitude_rank').head(n_top).copy()
    df['interaction'] = df.apply(
        lambda r: f"{r['source']} -> {r['target']}: {r['ligand_complex']}-{r['receptor_complex']}",
        axis=1
    )
    df['neg_log_mag'] = -np.log10(df['magnitude_rank'] + 1e-10)

    fig, ax = plt.subplots(figsize=(13, 10))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(df)))
    bars = ax.barh(range(len(df)), df['neg_log_mag'].values, color=colors,
                   edgecolor='black', linewidth=0.8)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['interaction'].values, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel('-log10(magnitude_rank)  (higher = stronger interaction)',
                  fontsize=11, fontweight='bold')
    ax.set_title(f'Top {n_top} Validated Cell-Cell Interactions (LIANA+ consensus)',
                 fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig1_Top_Interactions.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_celltype_heatmap(pair_counts, pair_strength):
    """Two-panel heatmap: counts + strength."""
    log.info("Plotting cell type pair heatmaps...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    sns.heatmap(pair_counts, annot=True, fmt='d', cmap='YlOrRd',
                cbar_kws={'label': 'N interactions'},
                ax=axes[0], linewidths=0.5, square=True)
    axes[0].set_title('A. Significant Interactions per Cell Type Pair',
                      fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Target', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Source', fontsize=11, fontweight='bold')

    sns.heatmap(pair_strength, annot=True, fmt='.1f', cmap='viridis',
                cbar_kws={'label': 'Aggregate strength'},
                ax=axes[1], linewidths=0.5, square=True)
    axes[1].set_title('B. Aggregate Interaction Strength',
                      fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Target', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Source', fontsize=11, fontweight='bold')

    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig2_CellType_Heatmaps.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_dotplot(sig_df, n_top_per_group=8):
    """Dotplot: ligand-receptor x source-target, dot size=magnitude, color=specificity."""
    log.info("Plotting LR dotplot...")
    df = sig_df.copy()
    df['lr'] = df['ligand_complex'] + '-' + df['receptor_complex']
    df['st'] = df['source'] + '->' + df['target']

    top_lrs = df.groupby('lr')['magnitude_rank'].min().sort_values().head(20).index
    df_plot = df[df['lr'].isin(top_lrs)].copy()
    df_plot['neg_log_mag'] = -np.log10(df_plot['magnitude_rank'] + 1e-10)
    df_plot['neg_log_spec'] = -np.log10(df_plot['specificity_rank'] + 1e-10)

    fig, ax = plt.subplots(figsize=(14, 9))
    scatter = ax.scatter(
        df_plot['st'], df_plot['lr'],
        s=df_plot['neg_log_mag'] * 60,
        c=df_plot['neg_log_spec'],
        cmap='plasma', edgecolors='black', linewidth=0.5, alpha=0.85,
    )
    ax.set_xlabel('Source -> Target', fontsize=11, fontweight='bold')
    ax.set_ylabel('Ligand-Receptor', fontsize=11, fontweight='bold')
    ax.set_title('Top 20 L/R Pairs across Cell Type Combinations',
                 fontsize=12, fontweight='bold')
    plt.xticks(rotation=90, fontsize=8)
    plt.yticks(fontsize=9)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('-log10(specificity_rank)', fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig3_LR_Dotplot.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def fig_sender_receiver(pair_counts):
    """Bar plot: total interactions per cell type as sender vs receiver."""
    log.info("Plotting sender/receiver activity...")
    sender_total = pair_counts.sum(axis=1)
    receiver_total = pair_counts.sum(axis=0)

    fig, ax = plt.subplots(figsize=(12, 6))
    x = np.arange(len(sender_total))
    width = 0.4
    ax.bar(x - width/2, sender_total.values, width, label='Sender (outgoing)',
           color='#1f77b4', edgecolor='black')
    ax.bar(x + width/2, receiver_total.values, width, label='Receiver (incoming)',
           color='#d62728', edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(sender_total.index, fontsize=10)
    ax.set_ylabel('Number of significant interactions', fontsize=11, fontweight='bold')
    ax.set_title('Communication Activity per Cell Type', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "Fig4_Sender_Receiver.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    start = time.time()
    log.info("=" * 80)
    log.info("PHASE 2B (LIANA+) VISUALIZATIONS")
    log.info("=" * 80)

    sig_df = pd.read_csv(INPUT_DIR / "liana_significant_interactions.csv")
    pair_counts = pd.read_csv(INPUT_DIR / "celltype_pair_counts.csv", index_col=0)
    pair_strength = pd.read_csv(INPUT_DIR / "celltype_pair_strength.csv", index_col=0)

    fig_top_interactions(sig_df, n_top=25)
    fig_celltype_heatmap(pair_counts, pair_strength)
    fig_dotplot(sig_df)
    fig_sender_receiver(pair_counts)

    log.info(f"DONE in {time.time()-start:.1f}s. Output: {FIG_DIR}")


if __name__ == "__main__":
    main()
