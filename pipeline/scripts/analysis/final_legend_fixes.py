#!/usr/bin/env python3
"""
Final Legend & Layout Fixes
1. Fig2_Top_LR_Interactions - Make TALLER (more height for Y-axis)
2. Fig4_Boundary_Validation - Panel A: legend box without color + dot_size legend
3. Fig5_Spatial_Trajectories - Panel A: legend box without color + dot_size legend
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import logging

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')
log = logging.getLogger(__name__)

REFINEMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_refinement"
CCC_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_hybrid_method"

FIG_CCC = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase2b_ccc_hybrid"
FIG_REGION_REF = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"

for d in [FIG_CCC, FIG_REGION_REF]:
    d.mkdir(parents=True, exist_ok=True)


def fix_lr_interactions(ccc_df):
    """FIX 1: Top LR Interactions - MUCH TALLER for legible Y labels"""
    log.info("Fixing Fig1: Top LR Interactions (MUCH taller)...")

    # Get top 20 by CCC score
    top_20 = ccc_df.nlargest(20, 'ccc_score')

    # CREATE MUCH TALLER FIGURE
    fig, ax = plt.subplots(figsize=(14, 18))  # Very tall

    pairs = [f"({row['source_cell_type']})→({row['target_cell_type']}): {row['ligand']}-{row['receptor']}"
             for _, row in top_20.iterrows()]
    scores = top_20['ccc_score'].values

    bars = ax.barh(range(len(scores)), scores, color='#1f77b4', alpha=0.8, edgecolor='black', linewidth=1)

    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize=11)  # Larger font + more space
    ax.set_xlabel('CCC Score', fontsize=12, fontweight='bold')
    ax.set_title('Top 20 Ligand-Receptor Interactions (By CCC Score)', fontsize=13, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig


def fix_boundary_validation_proper_legend(boundary_df):
    """FIX 2: Boundary Validation - Panel A: Legend box with dot_size + color info"""
    log.info("Fixing Fig2: Boundary Validation (proper legend box)...")

    fig = plt.figure(figsize=(20, 6.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.7)  # More space

    # Plot 1: Internal vs boundary distance (legend OUTSIDE in box)
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(boundary_df['internal_distance_mean'],
                         boundary_df['min_distance_to_other'],
                         s=boundary_df['n_cells']/100, alpha=0.6,
                         c=range(len(boundary_df)), cmap='viridis',
                         edgecolors='black', linewidth=1)

    ax1.set_xlabel('Internal Distance (μm)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Distance to Other Regions (μm)', fontsize=11, fontweight='bold')
    ax1.set_title('A. Region Boundary Clarity', fontsize=12, fontweight='bold')
    ax1.grid(alpha=0.3)

    # Create PROPER legend box with bullet points (NO COLOR BACKGROUND)
    legend_lines = [
        "• Bubble size = cell count in region",
        "• Color = region cluster ID",
        "  (viridis colormap)",
        "",
        "Interpretation:",
        "• Higher Y = clearer boundary",
        "• Larger bubble = more cells"
    ]
    legend_text = "\n".join(legend_lines)

    ax1.text(1.18, 0.5, legend_text, transform=ax1.transAxes,
            fontsize=9.5, verticalalignment='center', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='none', edgecolor='black',
                     linewidth=1.5, pad=1), family='monospace')

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
               s=150, alpha=0.6, c=range(len(boundary_df)), cmap='viridis',
               edgecolors='black', linewidth=1.2)

    ax3.set_xlabel('Clarity Ratio', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Spatial Dispersion (μm)', fontsize=11, fontweight='bold')
    ax3.set_title('C. Dispersion vs Clarity', fontsize=12, fontweight='bold')
    ax3.grid(alpha=0.3)

    plt.tight_layout()
    return fig


def fix_trajectories_proper_legend(trajectory_df):
    """FIX 3: Spatial Trajectories - Panel A: Legend box with dot_size + color info"""
    log.info("Fixing Fig3: Spatial Trajectories (proper legend box)...")

    fig = plt.figure(figsize=(18, 6.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.6)  # More space

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

    # Create PROPER legend box with bullet points (NO COLOR BACKGROUND)
    legend_lines = [
        "• Bubble size = cell count",
        "  (average of adjacent regions)",
        "",
        "• Color gradient = sharpness:",
        "  Red (hot) = sharp transition",
        "  Blue (cool) = gradual transition",
        "",
        "Interpretation:",
        "• Right-upper = smooth boundary",
        "  with preserved expression",
        "• Closer to axis = abrupt change"
    ]
    legend_text = "\n".join(legend_lines)

    ax1.text(1.20, 0.5, legend_text, transform=ax1.transAxes,
            fontsize=8.5, verticalalignment='center', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='none', edgecolor='black',
                     linewidth=1.5, pad=1), family='monospace')

    # Plot 2: Transition sharpness comparison
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

    plt.tight_layout()
    return fig


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("FINAL LEGEND & LAYOUT FIXES")
    log.info("=" * 80)

    try:
        # FIX 1: Top LR Interactions
        log.info("\n[1/3] Fixing Top LR Interactions (TALLER)...")
        ccc_df = pd.read_csv(CCC_OUTPUT_DIR / "ccc_pairwise_interactions.csv")
        fig1 = fix_lr_interactions(ccc_df)
        fig1.savefig(FIG_CCC / "Fig2_Top_LR_Interactions.png", dpi=300, bbox_inches='tight')
        plt.close(fig1)

        # FIX 2: Boundary Validation
        log.info("\n[2/3] Fixing Boundary Validation (proper legend)...")
        boundary_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_boundary_validation.csv")
        fig2 = fix_boundary_validation_proper_legend(boundary_df)
        fig2.savefig(FIG_REGION_REF / "Fig4_Boundary_Validation.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)

        # FIX 3: Trajectories
        log.info("\n[3/3] Fixing Trajectories (proper legend)...")
        trajectory_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_transitions.csv")
        fig3 = fix_trajectories_proper_legend(trajectory_df)
        fig3.savefig(FIG_REGION_REF / "Fig5_Spatial_Trajectories.png", dpi=300, bbox_inches='tight')
        plt.close(fig3)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ ALL FIXES COMPLETE in {elapsed:.1f} seconds")
        log.info("=" * 80)
        log.info("\nFixed files:")
        log.info("  • Fig2_Top_LR_Interactions.png (TALLER - Y labels legible)")
        log.info("  • Fig4_Boundary_Validation.png (legend box, no color)")
        log.info("  • Fig5_Spatial_Trajectories.png (legend box, no color)")

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
