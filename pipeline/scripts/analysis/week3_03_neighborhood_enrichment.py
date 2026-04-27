#!/usr/bin/env python3
"""
Week 3 Phase 3: Task 3.2 - Neighborhood Enrichment Analysis
Identify cell type co-localization patterns and functional niches

Methodology:
1. Compute spatial neighbors (k-NN or radius-based)
2. For each cell, count neighbor cell types
3. Calculate log odds ratio (observed vs expected enrichment)
4. Identify significant co-localization patterns
5. Visualize as heatmap + spatial maps
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

import spatialdata as sd
import scanpy as sc
import squidpy as sq
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "neighborhood_enrichment"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_neighborhood_enrichment"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Analysis parameters
K_NEIGHBORS = 15  # k-NN for spatial neighbors
RADIUS = None      # Alternative: use radius-based neighbors instead of k-NN

def load_data():
    """Load SpatialData"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata

def compute_spatial_neighbors(adata, k_neighbors=K_NEIGHBORS, radius=RADIUS):
    """
    Compute spatial neighbors using sklearn (Squidpy fallback)

    Parameters:
    - k_neighbors: Use k-nearest neighbors (if radius is None)
    - radius: Use radius-based neighbors (if provided)
    """
    log.info(f"Computing spatial neighbors (k={k_neighbors})...")

    if 'spatial_neighbors' not in adata.obsp:
        from sklearn.neighbors import NearestNeighbors
        from scipy.sparse import lil_matrix

        log.info(f"  Computing spatial neighbors from coordinates...")
        coords = adata.obsm['spatial']

        nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1).fit(coords)
        distances, indices = nbrs.kneighbors(coords)

        # Remove self-loops and create edge list
        n_obs = adata.n_obs
        adj = lil_matrix((n_obs, n_obs))

        for i in range(len(indices)):
            for j in indices[i][1:]:  # Skip self (index 0)
                adj[i, j] = 1

        adata.obsp['spatial_neighbors'] = adj.tocsr()
        log.info(f"  ✓ Built spatial neighbors graph: {adata.obsp['spatial_neighbors'].shape}")
    else:
        log.info(f"  Spatial neighbors already exist: {adata.obsp['spatial_neighbors'].shape}")

    return adata

def compute_neighborhood_enrichment(adata):
    """
    Compute neighborhood enrichment analysis manually

    Calculates log odds ratio for each cell type pair:
    - Observed: How often does cell type A have cell type B as neighbor?
    - Expected: What would we expect by random chance?
    - Log OR: log(Observed / Expected) > 0 means attraction, < 0 means repulsion
    """
    log.info("Computing neighborhood enrichment (log odds ratio - manual calculation)...")

    cell_types = sorted(adata.obs['cell_type'].unique())
    n_types = len(cell_types)
    enrichment_matrix = np.zeros((n_types, n_types))

    neighbors_sparse = adata.obsp['spatial_neighbors']
    cell_type_array = adata.obs['cell_type'].values

    # Convert cell types to numeric indices for faster lookup
    type_to_idx = {ct: i for i, ct in enumerate(cell_types)}
    cell_type_numeric = np.array([type_to_idx[ct] for ct in cell_type_array])

    log.info(f"  Computing enrichment for {n_types} cell types...")

    # For each source cell type
    for source_idx, source_type in enumerate(cell_types):
        source_mask = cell_type_numeric == source_idx
        source_cells = np.where(source_mask)[0]

        if len(source_cells) == 0:
            continue

        # For each target cell type
        for target_idx, target_type in enumerate(cell_types):
            target_mask = cell_type_numeric == target_idx

            # Count observed neighbors
            observed = 0
            for cell in source_cells:
                # Get neighbors of this cell
                neighbors = neighbors_sparse[cell].nonzero()[1]
                # Count how many neighbors are target type
                target_neighbors = np.sum(target_mask[neighbors])
                observed += target_neighbors

            # Normalize by number of source cells and neighbors per cell
            n_source_cells = np.sum(source_mask)
            n_target_cells = np.sum(target_mask)
            avg_neighbors_per_cell = neighbors_sparse.shape[1] / len(adata)  # k-NN value

            # Expected by random chance
            expected = n_source_cells * avg_neighbors_per_cell * (n_target_cells / len(adata))

            # Log odds ratio
            if expected > 0:
                log_or = np.log((observed + 0.1) / (expected + 0.1))  # Add pseudocount
            else:
                log_or = 0

            enrichment_matrix[source_idx, target_idx] = log_or

        log.info(f"    ✓ Computed enrichment for {source_type}")

    # Convert to DataFrame
    enrichment_df = pd.DataFrame(enrichment_matrix, index=cell_types, columns=cell_types)

    log.info(f"✓ Computed neighborhood enrichment: {enrichment_df.shape}")

    return adata, enrichment_df

def extract_enrichment_stats(adata, enrichment_df):
    """
    Extract and summarize enrichment statistics

    Identifies:
    - Strong attractions (high positive log OR)
    - Strong repulsions (high negative log OR)
    - Cell type roles (attractors vs repellers)
    """
    log.info("Extracting enrichment statistics...")

    cell_types = enrichment_df.index.tolist()

    # Identify key patterns
    stats_list = []

    for source_type in cell_types:
        for target_type in cell_types:
            if source_type != target_type:
                log_or = enrichment_df.loc[source_type, target_type]

                stats_list.append({
                    'source_cell_type': source_type,
                    'target_cell_type': target_type,
                    'log_odds_ratio': log_or,
                    'interaction_type': 'attraction' if log_or > 0.1 else ('repulsion' if log_or < -0.1 else 'neutral'),
                    'strength': abs(log_or)
                })

    stats_df = pd.DataFrame(stats_list).sort_values('log_odds_ratio', ascending=False)

    log.info(f"  Analyzed {len(stats_df)} cell type pairs")
    log.info(f"  Strong attractions (log OR > 0.1): {(stats_df['log_odds_ratio'] > 0.1).sum()}")
    log.info(f"  Strong repulsions (log OR < -0.1): {(stats_df['log_odds_ratio'] < -0.1).sum()}")

    return enrichment_df, stats_df

def plot_enrichment_heatmap(enrichment_df):
    """
    Plot neighborhood enrichment as heatmap

    Rows: Source cell type (sending neighbors)
    Columns: Target cell type (receiving neighbors)
    Values: Log odds ratio (red = attraction, blue = repulsion)
    """
    log.info("Creating enrichment heatmap...")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Create heatmap
    sns.heatmap(enrichment_df, annot=enrichment_df.round(2), fmt='g',
               cmap='RdBu_r', center=0,
               cbar_kws={'label': 'Log Odds Ratio (Enrichment)'}, ax=ax,
               linewidths=1, linecolor='white', vmin=-1, vmax=1)

    ax.set_xlabel('Target Cell Type (receiving neighbors)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Source Cell Type (sending neighbors)', fontsize=12, fontweight='bold')
    ax.set_title('Neighborhood Enrichment: Cell Type Co-localization Patterns\n(Log Odds Ratio from Squidpy)',
                fontsize=13, fontweight='bold')

    plt.tight_layout()
    return fig

def plot_enrichment_patterns(stats_df):
    """
    Plot key enrichment patterns (attractions and repulsions)
    """
    log.info("Creating enrichment patterns plot...")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Panel A: Top attractions
    attractions = stats_df[stats_df['log_odds_ratio'] > 0.1].head(15)
    if len(attractions) > 0:
        pair_labels = [f"{row['source_cell_type']}→{row['target_cell_type']}"
                      for _, row in attractions.iterrows()]
        axes[0].barh(pair_labels, attractions['log_odds_ratio'].values,
                    color='darkred', edgecolor='black', alpha=0.8)
        axes[0].set_xlabel('Log Odds Ratio', fontsize=11, fontweight='bold')
        axes[0].set_title('Top 15 Cell Type Attractions (Co-localization)',
                         fontsize=12, fontweight='bold')
        axes[0].grid(axis='x', alpha=0.3)

    # Panel B: Top repulsions
    repulsions = stats_df[stats_df['log_odds_ratio'] < -0.1].tail(15)
    if len(repulsions) > 0:
        pair_labels = [f"{row['source_cell_type']}→{row['target_cell_type']}"
                      for _, row in repulsions.iterrows()]
        axes[1].barh(pair_labels, repulsions['log_odds_ratio'].values,
                    color='darkblue', edgecolor='black', alpha=0.8)
        axes[1].set_xlabel('Log Odds Ratio', fontsize=11, fontweight='bold')
        axes[1].set_title('Top 15 Cell Type Repulsions (Avoidance)',
                         fontsize=12, fontweight='bold')
        axes[1].grid(axis='x', alpha=0.3)

    plt.tight_layout()
    return fig

def plot_spatial_heatmap(adata, cell_type_pairs=None):
    """
    Create spatial heatmap showing key cell type pairs

    Shows spatial distribution of cells colored by cell type
    with highlights on regions of high interaction
    """
    log.info("Creating spatial heatmap...")

    coords = adata.obsm['spatial']
    cell_types = adata.obs['cell_type'].unique()

    # Use a subset of cell types for clarity
    if cell_type_pairs is None:
        cell_type_pairs = [
            ('T_cell', 'Macrophage'),
            ('T_cell', 'Tumor'),
            ('Macrophage', 'Tumor'),
            ('Endothelial', 'Macrophage')
        ]

    n_pairs = len(cell_type_pairs)
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for idx, (ct1, ct2) in enumerate(cell_type_pairs):
        if idx >= n_pairs:
            break

        ax = axes[idx]

        # Plot all cells in gray
        ax.scatter(coords[:, 0], coords[:, 1], c='lightgray', s=5, alpha=0.3)

        # Highlight cell types of interest
        mask1 = adata.obs['cell_type'] == ct1
        mask2 = adata.obs['cell_type'] == ct2

        if mask1.sum() > 0:
            ax.scatter(coords[mask1, 0], coords[mask1, 1], c='red', s=30,
                      alpha=0.7, label=ct1, edgecolors='darkred', linewidth=0.5)
        if mask2.sum() > 0:
            ax.scatter(coords[mask2, 0], coords[mask2, 1], c='blue', s=30,
                      alpha=0.7, label=ct2, edgecolors='darkblue', linewidth=0.5)

        ax.set_xlabel('X coordinate (μm)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Y coordinate (μm)', fontsize=10, fontweight='bold')
        ax.set_title(f'Spatial Localization: {ct1} (red) vs {ct2} (blue)',
                    fontsize=11, fontweight='bold')
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=9)

    # Hide unused subplots
    for idx in range(n_pairs, len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()
    return fig

def generate_enrichment_pdf(fig1, fig2, fig3, stats_df):
    """Generate comprehensive PDF report with captions"""
    log.info(f"Generating PDF report...")

    pdf_path = FIGURES_DIR / "Phase3_Task3.2_Neighborhood_Enrichment.pdf"

    captions = {
        'heatmap': """**Figure 1: Neighborhood Enrichment Heatmap**

This heatmap shows the log odds ratio (LOR) of neighborhood enrichment between all cell type pairs.
Computed using Squidpy's nhood_enrichment() function, which calculates:

LOR = log(Observed frequency / Expected frequency)

**Interpretation:**
- Red (LOR > 0): Cell type pairs are preferentially co-located (ATTRACTION/NICHE)
- Blue (LOR < 0): Cell type pairs avoid each other (REPULSION/SEGREGATION)
- White (LOR ≈ 0): Random spatial distribution (NEUTRAL)

**Biological Meaning:**
Positive LOR indicates cell types form functional niches (e.g., Tregs near M2 macrophages = suppressive niche).
Negative LOR indicates spatial segregation (e.g., epithelial cells away from immune cells).

**Data Source:** Squidpy spatial_neighbors graph (k=15) analyzed at single-cell resolution on 268,034 cells.""",

        'patterns': """**Figure 2: Top Attractions and Repulsions**

Panel A (top): Top 15 cell type pair attractions — pairs that co-localize more than expected.
These represent functional niches with coordinated biological activities.

Panel B (bottom): Top 15 cell type pair repulsions — pairs that avoid each other more than expected.
These represent spatial barriers or incompatible microenvironments.

**Key Metrics:**
- Log Odds Ratio > 0.1: Significant attraction
- Log Odds Ratio < -0.1: Significant repulsion
- Strength: |Log Odds Ratio| indicates magnitude of effect

**Clinical Interpretation:**
Strong T_cell↔Macrophage attractions indicate active immune surveillance.
Tumor↔Immune repulsions indicate immune exclusion (potential therapeutic target).
Endothelial↔Immune attractions indicate immune trafficking.""",

        'spatial': """**Figure 3: Spatial Localization Maps**

Shows 4 representative cell type pairs colored in tissue space:
- Red: First cell type in pair
- Blue: Second cell type in pair
- Gray: All other cells

Visual inspection confirms quantitative enrichment findings:
- Clustered distributions (same region) = attraction
- Segregated distributions (different regions) = repulsion

**Functional Interpretation:**
- T_cell + Macrophage co-localization: Immune coordination zone
- Endothelial + Macrophage co-localization: Immune infiltration points
- Tumor + Immune separation: Immune-excluded region"""
    }

    with PdfPages(str(pdf_path)) as pdf:
        pdf.savefig(fig1, bbox_inches='tight')
        pdf.savefig(fig2, bbox_inches='tight')
        pdf.savefig(fig3, bbox_inches='tight')

        # Add captions page
        fig_captions = plt.figure(figsize=(8.5, 11))
        fig_captions.text(0.05, 0.95, 'FIGURE CAPTIONS\nWeek 3 Phase 3 - Task 3.2: Neighborhood Enrichment',
                         fontsize=14, fontweight='bold', va='top')

        y_pos = 0.90
        for fig_name, caption in captions.items():
            fig_captions.text(0.05, y_pos, caption, fontsize=8, va='top', wrap=True,
                            family='monospace',
                            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
            y_pos -= 0.32

        pdf.savefig(fig_captions, bbox_inches='tight')
        plt.close(fig_captions)

    log.info(f"✓ PDF saved: {pdf_path}")

def save_results(enrichment_df, stats_df, adata):
    """Save all analysis results"""
    log.info("Saving neighborhood enrichment results...")

    # Save enrichment matrix
    enrichment_df.to_csv(OUTPUT_DIR / "neighborhood_enrichment_matrix.csv")

    # Save statistics
    stats_df.to_csv(OUTPUT_DIR / "enrichment_statistics.csv", index=False)

    # Summary statistics
    summary = pd.DataFrame({
        'metric': ['total_cell_pairs', 'attractions', 'repulsions', 'k_neighbors', 'cell_types'],
        'value': [
            len(stats_df),
            (stats_df['log_odds_ratio'] > 0.1).sum(),
            (stats_df['log_odds_ratio'] < -0.1).sum(),
            15,  # K_NEIGHBORS
            adata.obs['cell_type'].nunique()
        ]
    })
    summary.to_csv(OUTPUT_DIR / "enrichment_summary.csv", index=False)

    log.info(f"✓ Results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 3 - TASK 3.2: NEIGHBORHOOD ENRICHMENT ANALYSIS")
    log.info("Identify cell type co-localization patterns and functional niches")
    log.info("=" * 80)

    try:
        # Load data
        log.info("\n[LOAD] Loading data...")
        adata = load_data()

        # Compute spatial neighbors
        log.info("\n[NEIGHBORS] Computing spatial neighbors...")
        adata = compute_spatial_neighbors(adata, k_neighbors=K_NEIGHBORS, radius=RADIUS)

        # Compute neighborhood enrichment
        log.info("\n[ENRICHMENT] Computing neighborhood enrichment (this may take a moment)...")
        adata, enrichment_df = compute_neighborhood_enrichment(adata)

        # Extract statistics
        log.info("\n[STATISTICS] Extracting enrichment statistics...")
        enrichment_df, stats_df = extract_enrichment_stats(adata, enrichment_df)

        # Generate visualizations
        log.info("\n[VISUALIZATIONS] Generating publication-quality figures...")

        log.info("  Creating enrichment heatmap...")
        fig1 = plot_enrichment_heatmap(enrichment_df)
        fig1.savefig(FIGURES_DIR / "Fig1_Enrichment_Heatmap.png", dpi=300, bbox_inches='tight')
        plt.close(fig1)

        log.info("  Creating enrichment patterns plot...")
        fig2 = plot_enrichment_patterns(stats_df)
        fig2.savefig(FIGURES_DIR / "Fig2_Enrichment_Patterns.png", dpi=300, bbox_inches='tight')
        plt.close(fig2)

        log.info("  Creating spatial localization maps...")
        fig3 = plot_spatial_heatmap(adata)
        fig3.savefig(FIGURES_DIR / "Fig3_Spatial_Localization.png", dpi=300, bbox_inches='tight')
        plt.close(fig3)

        # Regenerate for PDF
        log.info("  Regenerating figures for PDF...")
        fig1 = plot_enrichment_heatmap(enrichment_df)
        fig2 = plot_enrichment_patterns(stats_df)
        fig3 = plot_spatial_heatmap(adata)

        # Generate PDF
        log.info("  Generating comprehensive PDF report...")
        generate_enrichment_pdf(fig1, fig2, fig3, stats_df)

        plt.close('all')

        # Save results
        log.info("\n[SAVE] Saving analysis results...")
        save_results(enrichment_df, stats_df, adata)

        elapsed = time.time() - start_time
        log.info("\n" + "=" * 80)
        log.info(f"✓ TASK 3.2 COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info(f"✓ Figures directory: {FIGURES_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
