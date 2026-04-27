#!/usr/bin/env python3
"""
Compact Visual Legends for Panels A
- Fig4_Boundary_Validation: Panel A
- Fig5_Spatial_Trajectories: Panel A

Use shapes and colors instead of text for compact, visual design.
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
FIG_REGION_REF = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_refinement"

FIG_REGION_REF.mkdir(parents=True, exist_ok=True)


def fix_boundary_validation_compact_legend(boundary_df):
    """Boundary Validation - Panel A with COMPACT VISUAL legend"""
    log.info("Fixing Fig4: Boundary Validation (compact visual legend)...")

    fig = plt.figure(figsize=(20, 6.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.75)

    # Plot 1: Internal vs boundary distance (legend VERY COMPACT)
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

    # VERY COMPACT VISUAL LEGEND (vertical, minimal text)
    legend_y_start = 0.98
    legend_x = 1.15

    # Title
    ax1.text(legend_x, legend_y_start, 'Legend:', transform=ax1.transAxes,
            fontsize=10, fontweight='bold', verticalalignment='top', horizontalalignment='left')

    # Bubble sizes (visual)
    bubble_sizes = [2000, 5000, 10000]
    bubble_labels = ['Small', 'Med', 'Large']
    y_pos = legend_y_start - 0.08

    for size, label in zip(bubble_sizes, bubble_labels):
        ax1.scatter([], [], s=size/100, c='gray', alpha=0.6, edgecolors='black', linewidth=1,
                   transform=ax1.transAxes, label=f'{label}')
        ax1.text(legend_x + 0.04, y_pos, label, transform=ax1.transAxes,
                fontsize=9, verticalalignment='center', horizontalalignment='left')
        y_pos -= 0.06

    # Add small note at bottom
    ax1.text(legend_x, y_pos - 0.02, 'bubble =\ncell count', transform=ax1.transAxes,
            fontsize=8, style='italic', verticalalignment='top', horizontalalignment='left')

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


def fix_trajectories_compact_legend(trajectory_df):
    """Spatial Trajectories - Panel A with COMPACT VISUAL legend"""
    log.info("Fixing Fig5: Spatial Trajectories (compact visual legend)...")

    fig = plt.figure(figsize=(18, 6.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.65)

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

    # VERY COMPACT VISUAL LEGEND (vertical, minimal text, use colors)
    legend_y_start = 0.98
    legend_x = 1.13

    # Title
    ax1.text(legend_x, legend_y_start, 'Legend:', transform=ax1.transAxes,
            fontsize=10, fontweight='bold', verticalalignment='top', horizontalalignment='left')

    # Bubble sizes (visual)
    bubble_sizes = [30000, 60000, 100000]
    y_pos = legend_y_start - 0.08

    ax1.scatter([], [], s=bubble_sizes[0]/1500, c='gray', alpha=0.6, edgecolors='black', linewidth=1,
               transform=ax1.transAxes)
    ax1.text(legend_x + 0.03, y_pos, 'Small', transform=ax1.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='left')

    y_pos -= 0.06
    ax1.scatter([], [], s=bubble_sizes[1]/1500, c='gray', alpha=0.6, edgecolors='black', linewidth=1,
               transform=ax1.transAxes)
    ax1.text(legend_x + 0.03, y_pos, 'Med', transform=ax1.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='left')

    y_pos -= 0.06
    ax1.scatter([], [], s=bubble_sizes[2]/1500, c='gray', alpha=0.6, edgecolors='black', linewidth=1,
               transform=ax1.transAxes)
    ax1.text(legend_x + 0.03, y_pos, 'Large', transform=ax1.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='left')

    # Color gradient (sharpness)
    y_pos -= 0.08
    ax1.text(legend_x, y_pos, 'Sharpness:', transform=ax1.transAxes,
            fontsize=8, fontweight='bold', verticalalignment='top', horizontalalignment='left')

    y_pos -= 0.06
    # Red (sharp)
    rect_sharp = mpatches.Rectangle((legend_x-0.01, y_pos-0.015), 0.015, 0.03,
                                   transform=ax1.transAxes, facecolor='#d62728', edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect_sharp)
    ax1.text(legend_x + 0.025, y_pos, 'Sharp', transform=ax1.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='left')

    y_pos -= 0.05
    # Blue (gradual)
    rect_gradual = mpatches.Rectangle((legend_x-0.01, y_pos-0.015), 0.015, 0.03,
                                     transform=ax1.transAxes, facecolor='#1f77b4', edgecolor='black', linewidth=0.5)
    ax1.add_patch(rect_gradual)
    ax1.text(legend_x + 0.025, y_pos, 'Gradual', transform=ax1.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='left')

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
    log.info("COMPACT VISUAL LEGENDS FOR PANELS A")
    log.info("=" * 80)

    try:
        # FIX: Boundary Validation
        log.info("\n[1/2] Fixing Boundary Validation (compact legend)...")
        boundary_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_boundary_validation.csv")
        fig1 = fix_boundary_validation_compact_legend(boundary_df)
        fig1.savefig(FIG_REGION_REF / "Fig4_Boundary_Validation.png", dpi=300, bbox_inches='tight')
        plt.close(fig1)

        # FIX: Trajectories
        log.info("\n[2/2] Fixing Trajectories (compact legend)...")
        trajectory_df = pd.read_csv(REFINEMENT_OUTPUT_DIR / "region_transitions.csv")
        fig2 = fix_trajectories_compact_legend(trajectory_df)
        fig2.savefig(FIG_REGION_REF / "Fig5_Spatial_Trajectories.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ COMPACT LEGENDS COMPLETE in {elapsed:.1f} seconds")
        log.info("=" * 80)
        log.info("\nUpdated files:")
        log.info("  • Fig4_Boundary_Validation.png (Panel A: compact visual legend)")
        log.info("  • Fig5_Spatial_Trajectories.png (Panel A: compact visual legend)")

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
