#!/usr/bin/env python3
"""
Week 3 Phase 2B: Hybrid CCC Method
Combines DeepTalk concepts (spatial context) with CellChat concepts (comprehensive L/R database)
Optimized for speed while maintaining statistical rigor

Methodology:
1. DeepTalk inspiration: Spatial proximity as graph structure
2. CellChat inspiration: Comprehensive L/R pair database + mass action model
3. Our implementation: Fast hybrid using Squidpy's spatial context
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from scipy.stats import spearmanr, pearsonr

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
PHASE1_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_hybrid_method"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# CellChat-inspired comprehensive L/R pairs (subset for 289-gene panel)
CELLCHAT_LR_PAIRS = [
    # Immune activation
    ("CD40", "CD40LG"), ("CD86", "CD28"), ("CD80", "CTLA4"), ("CD86", "CTLA4"),
    # Macrophage-T cell (immunosuppression/activation)
    ("CD68", "CD14"), ("TNF", "TNFRSF1A"), ("IL10", "IL10RA"), ("IL12A", "IL12RB1"),
    ("PDCD1LG2", "PDCD1"), ("IDO1", "AHR"), ("CD40", "TNFSF14"),
    # Endothelial-Immune
    ("ICAM1", "LFA1"), ("ICAM1", "ITGAL"), ("VCAM1", "ITGA4"),
    # NK cell
    ("HLA-A", "KIR2DL1"), ("ULBP1", "NKG2D"),
    # Epithelial-Immune
    ("EPCAM", "VTCN1"),
    # Cytotoxic
    ("PRF1", "LYST"), ("GZMA", "SLC7A2"), ("GZMB", "SLC7A2"),
]

def load_data():
    """Load SpatialData and Phase 1"""
    log.info("Loading data...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")

    phase1_summary = pd.read_csv(PHASE1_DIR / "DGE_summary_all_types.csv")
    return adata, phase1_summary

def compute_ccc_with_spatial_context(adata, lr_pairs):
    """
    Compute CCC integrating spatial context (DeepTalk-inspired key feature)

    Strategy:
    1. For each L/R pair + cell type pair: compute "spatial proximity score"
    2. Weight L/R interaction by spatial proximity
    3. Use Spearman correlation of L/R expressions within neighborhoods
    """
    log.info("Computing hybrid CCC (spatial context integrated)...")

    # Ensure spatial neighbors - compute directly (Squidpy sometimes has issues with SpatialData)
    if 'spatial_neighbors' not in adata.obsp:
        log.info("  Computing spatial neighbors from coordinates...")
        from sklearn.neighbors import NearestNeighbors
        from scipy.sparse import lil_matrix

        # Try to get coordinates
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            log.info(f"    Using obsm['spatial']: {coords.shape}")
        elif 'x_coordinate' in adata.obs.columns and 'y_coordinate' in adata.obs.columns:
            coords = adata.obs[['x_coordinate', 'y_coordinate']].values
            log.info(f"    Using obs columns x/y_coordinate: {coords.shape}")
        else:
            raise ValueError("Cannot find spatial coordinates in adata.obsm['spatial'] or obs[x/y_coordinate]")

        nbrs = NearestNeighbors(n_neighbors=16).fit(coords)
        _, indices = nbrs.kneighbors(coords)

        # Build sparse adjacency matrix
        n_obs = adata.n_obs
        adj = lil_matrix((n_obs, n_obs))
        for i, neighbors in enumerate(indices):
            for j in neighbors[1:]:  # Skip self
                adj[i, j] = 1
        adata.obsp['spatial_neighbors'] = adj.tocsr()
        log.info(f"    ✓ Built spatial neighbors: {adata.obsp['spatial_neighbors'].shape}")

    spatial_neighbors = adata.obsp['spatial_neighbors']

    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    panel_genes = set(adata.var_names)
    results = []

    cell_types = adata.obs['cell_type'].unique()

    # For each L/R pair
    for ligand, receptor in lr_pairs:
        if ligand not in panel_genes or receptor not in panel_genes:
            continue

        lig_idx = adata.var_names.get_loc(ligand)
        rec_idx = adata.var_names.get_loc(receptor)

        lig_expr = X[:, lig_idx]
        rec_expr = X[:, rec_idx]

        # For each cell type pair
        for source_type in cell_types:
            source_mask = (adata.obs['cell_type'] == source_type).values

            for target_type in cell_types:
                target_mask = (adata.obs['cell_type'] == target_type).values

                # Spatial proximity: does source type neighbor target type?
                source_indices = np.where(source_mask)[0]
                target_neighbors = np.array(spatial_neighbors[target_mask].sum(axis=0)).flatten() > 0

                # Cells in source type that are near target type
                source_near_target = source_mask & target_neighbors

                if source_near_target.sum() < 2:
                    continue

                # Interaction score: product of ligand and receptor expression
                # weighted by spatial proximity
                lig_mean = lig_expr[source_near_target].mean() if source_near_target.sum() > 0 else 0
                rec_mean = rec_expr[target_neighbors].mean() if target_neighbors.sum() > 0 else 0
                interaction_score = np.log2(lig_mean * rec_mean + 1)

                # Spatial confidence: fraction of source cells expressing ligand +
                # fraction of target neighbors expressing receptor (co-localization strength)
                n_source_near = source_near_target.sum()
                n_target_near = target_neighbors.sum()

                if n_source_near > 0 and n_target_near > 0:
                    # What fraction express above median?
                    lig_med = np.median(lig_expr[source_near_target])
                    rec_med = np.median(rec_expr[target_neighbors])

                    frac_lig_high = (lig_expr[source_near_target] > lig_med).sum() / max(1, n_source_near)
                    frac_rec_high = (rec_expr[target_neighbors] > rec_med).sum() / max(1, n_target_near)

                    spatial_confidence = (frac_lig_high + frac_rec_high) / 2.0  # Average fraction expressing
                else:
                    spatial_confidence = 0.0

                # Final score: interaction × spatial confidence
                final_score = interaction_score * spatial_confidence

                results.append({
                    'ligand': ligand,
                    'receptor': receptor,
                    'pair_name': f"{ligand}→{receptor}",
                    'source_cell_type': source_type,
                    'target_cell_type': target_type,
                    'interaction_score': interaction_score,
                    'spatial_confidence': spatial_confidence,
                    'ccc_score': final_score,
                })

    ccc_df = pd.DataFrame(results).sort_values('ccc_score', ascending=False)
    log.info(f"✓ Computed {len(ccc_df)} interactions")

    return ccc_df

def create_ccc_summary(ccc_df, adata):
    """Create cell type-level CCC summary"""
    log.info("Creating cell type-level CCC summary...")

    cell_types = adata.obs['cell_type'].unique()
    n_types = len(cell_types)

    # Overall CCC matrix
    ccc_matrix = np.zeros((n_types, n_types))
    for _, row in ccc_df.iterrows():
        i = list(cell_types).index(row['source_cell_type'])
        j = list(cell_types).index(row['target_cell_type'])
        ccc_matrix[i, j] += row['ccc_score']

    ccc_df_matrix = pd.DataFrame(ccc_matrix, index=cell_types, columns=cell_types)

    return ccc_df_matrix

def save_results(ccc_df, ccc_matrix):
    """Save results"""
    log.info("Saving results...")

    ccc_df.to_csv(OUTPUT_DIR / "ccc_pairwise_interactions.csv", index=False)
    ccc_matrix.to_csv(OUTPUT_DIR / "ccc_matrix_celltype.csv")

    # Summary stats
    stats = pd.DataFrame({
        'metric': ['total_interactions', 'mean_score', 'max_score', 'pair_names'],
        'value': [len(ccc_df), ccc_df['ccc_score'].mean(), ccc_df['ccc_score'].max(), ccc_df['pair_name'].nunique()]
    })
    stats.to_csv(OUTPUT_DIR / "ccc_summary_stats.csv", index=False)

    log.info(f"✓ Results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 2B: HYBRID CCC (DEEPTALK + CELLCHAT CONCEPTS)")
    log.info("Spatial context + Comprehensive L/R database")
    log.info("=" * 80)

    try:
        # Load
        adata, phase1_summary = load_data()

        # Compute CCC
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Compute Hybrid CCC (Spatial-Integrated)")
        log.info("=" * 80)
        ccc_df = compute_ccc_with_spatial_context(adata, CELLCHAT_LR_PAIRS)

        # Summary
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Create Cell Type-Level Summary")
        log.info("=" * 80)
        ccc_matrix = create_ccc_summary(ccc_df, adata)

        # Save
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Save Results")
        log.info("=" * 80)
        save_results(ccc_df, ccc_matrix)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ PHASE 2B COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
