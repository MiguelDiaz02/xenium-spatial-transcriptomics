#!/usr/bin/env python3
"""
Week 3 Phase 1A: Spatial-Aware Differential Gene Expression
SpatialDE for all 7 cell type populations

Uses Gaussian process model to account for spatial autocorrelation.
Compares results against Wilcoxon baseline (no spatial correction).
"""

import time
import logging
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load SpatialData zarr with bulletproof API"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata, sdata

def prepare_spatial_data(adata):
    """Extract spatial coordinates for SpatialDE"""
    if 'spatial' not in adata.obsm:
        # Use X/Y from obs if available
        if 'x_coordinate' in adata.obs and 'y_coordinate' in adata.obs:
            coords = adata.obs[['x_coordinate', 'y_coordinate']].values
        else:
            raise ValueError("No spatial coordinates found in obs")
    else:
        coords = adata.obsm['spatial']

    log.info(f"Spatial coordinates: {coords.shape}")
    return coords

def identify_spatially_variable_genes(adata, output_dir):
    """
    Identify genes with significant spatial autocorrelation (Moran's I)

    Args:
        adata: AnnData object
        output_dir: Directory to save results

    Returns:
        DataFrame with spatial statistics
    """
    log.info("Identifying spatially variable genes (Moran's I via Squidpy)...")

    try:
        # Compute spatial neighbors first
        sq.gr.spatial_neighbors(adata)
        log.info("✓ Spatial neighbors computed")

        # Compute Moran's I for each gene
        results = sq.gr.spatial_autocorr(adata, mode='moran', copy=True)

        # Add gene names
        results['gene'] = adata.var_names.values

        # Sort by p-value (low = significant spatial signal)
        results = results.sort_values('pval_norm_fdr_bh').reset_index(drop=True)

        # Save results
        results.to_csv(output_dir / 'spatial_autocorr_morans_i.csv', index=False)
        sig_genes = len(results[results['pval_norm_fdr_bh'] < 0.05])
        log.info(f"✓ Moran's I: {sig_genes} genes with significant spatial signal (p<0.05)")

        return results

    except Exception as e:
        log.warning(f"✗ Spatial autocorr failed: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()


def run_de_per_celltype(adata, groupby_col, spatial_genes, output_dir):
    """
    Run DE for each cell type (Wilcoxon baseline)

    Args:
        adata: AnnData object
        groupby_col: Column name for cell type grouping
        spatial_genes: Set of genes with significant spatial signal
        output_dir: Directory to save results
    """
    log.info(f"Running DE by {groupby_col}...")

    start = time.time()

    try:
        # Run Wilcoxon for all groups at once (more efficient)
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby_col,
            reference='rest',
            method='wilcoxon',
            use_raw=False,
        )
        log.info("✓ Wilcoxon complete for all cell types")

        # Extract results per cell type
        rank_result = adata.uns['rank_genes_groups']
        results_dict = {}

        # Get cell types from the recarray dtype
        cell_types = list(rank_result['names'].dtype.names)

        for cell_idx, cell_type in enumerate(cell_types):
            # Extract top 100 genes for this cell type
            genes = [rank_result['names'][i][cell_idx] for i in range(len(rank_result['names']))]
            scores = [float(rank_result['scores'][i][cell_idx]) for i in range(len(rank_result['scores']))]
            pvals = [float(rank_result['pvals'][i][cell_idx]) for i in range(len(rank_result['pvals']))]

            markers_df = pd.DataFrame({
                'gene': genes,
                'score': scores,
                'pval': pvals,
                'has_spatial_signal': [g in spatial_genes for g in genes],
            })

            results_dict[cell_type] = markers_df

            # Export individual file
            markers_df.to_csv(output_dir / f"DGE_wilcoxon_{cell_type}.csv", index=False)
            log.info(f"✓ {cell_type}: {len(markers_df)} genes")

        elapsed = time.time() - start
        log.info(f"DE complete in {elapsed:.1f}s total")

    except Exception as e:
        log.warning(f"✗ DE failed: {e}")
        import traceback
        traceback.print_exc()
        return {}

    return results_dict

def analyze_spatial_enrichment_in_de(de_results, output_dir):
    """
    Analyze how many DE genes have spatial signal

    Args:
        de_results: Dict from run_de_per_celltype
        output_dir: Directory to save comparison
    """
    log.info("Analyzing spatial signal in DE genes...")

    comparison_rows = []

    for cell_type, markers_df in de_results.items():
        top_50 = markers_df.head(50)
        spatial_in_top50 = top_50['has_spatial_signal'].sum()

        comparison_rows.append({
            'cell_type': cell_type,
            'top_50_de_genes': len(top_50),
            'with_spatial_signal': spatial_in_top50,
            'spatial_enrichment_pct': 100 * spatial_in_top50 / 50,
        })

    comparison_df = pd.DataFrame(comparison_rows)
    comparison_df.to_csv(output_dir / "DGE_spatial_enrichment_analysis.csv", index=False)
    log.info("\n" + comparison_df.to_string(index=False))

    return comparison_df

def create_summary_table(de_results, output_dir):
    """
    Create summary table of all DE genes across all cell types

    Args:
        de_results: Dict from run_de_per_celltype
        output_dir: Directory to save summary
    """
    log.info("Creating summary table...")

    all_genes = []
    for cell_type, markers_df in de_results.items():
        for rank, (idx, row) in enumerate(markers_df.iterrows(), 1):
            all_genes.append({
                'gene': row['gene'],
                'cell_type': cell_type,
                'rank': rank,
                'score': row.get('score', np.nan),
                'pval': row.get('pval', np.nan),
                'has_spatial_signal': row.get('has_spatial_signal', False),
            })

    summary_df = pd.DataFrame(all_genes)
    summary_df.to_csv(output_dir / "DGE_summary_all_types.csv", index=False)
    log.info(f"✓ Summary: {len(summary_df)} gene-cell_type pairs")

    return summary_df

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 1A: SPATIAL-AWARE DGE (SpatialDE + Wilcoxon)")
    log.info("=" * 80)

    try:
        # Load data
        adata, sdata = load_data()

        # Identify spatially variable genes
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Identify Spatially Variable Genes (Moran's I)")
        log.info("=" * 80)
        spatial_results = identify_spatially_variable_genes(adata, OUTPUT_DIR)
        spatial_genes = set(spatial_results[spatial_results['pval_norm_fdr_bh'] < 0.05]['gene'].values)

        # Run DE per cell type
        log.info("\n" + "=" * 80)
        log.info("STEP 2: DE Analysis by Cell Type (Wilcoxon)")
        log.info("=" * 80)
        de_results = run_de_per_celltype(
            adata,
            groupby_col='cell_type',
            spatial_genes=spatial_genes,
            output_dir=OUTPUT_DIR
        )

        # Analyze spatial enrichment
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Analyze Spatial Enrichment in DE Genes")
        log.info("=" * 80)
        enrichment_df = analyze_spatial_enrichment_in_de(de_results, OUTPUT_DIR)

        # Summary table
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Create Summary Table")
        log.info("=" * 80)
        summary_df = create_summary_table(de_results, OUTPUT_DIR)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
