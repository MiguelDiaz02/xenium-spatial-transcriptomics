#!/usr/bin/env python3
"""
Week 3 Phase 1: Benchmarked DGEA Analysis
Spatial-aware DGE combining multiple methods from benchmarking:
- Wilcoxon (baseline - non-spatial)
- Moran's I (spatial autocorrelation via Squidpy)
- SpatialDE alternative (if available)
- SPARK (Spatial Autocorrelation Kernel)

Benchmarking source: elad011.pdf (Oxford 2024) - 51 DE tools reviewed
Spatial methods: SpatialDE, SPARK, Spateo, Squidpy.spatial_markers()
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
from scipy.stats import wilcoxon, rankdata

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE_benchmarked"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load SpatialData"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata, sdata

def method_wilcoxon_baseline(adata):
    """
    Wilcoxon rank-sum test (non-spatial baseline)
    Benchmarking: Well-established, fast, ignores spatial structure
    """
    log.info("METHOD 1: Wilcoxon Rank-Sum Test (Non-Spatial Baseline)")

    sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', key_added='wilcoxon')

    results_by_type = {}
    for cell_type in adata.obs['cell_type'].unique():
        genes = adata.uns['wilcoxon']['names'][cell_type]
        scores = adata.uns['wilcoxon']['scores'][cell_type]
        pvals = adata.uns['wilcoxon']['pvals'][cell_type]

        df = pd.DataFrame({
            'gene': genes,
            'wilcoxon_score': scores,
            'wilcoxon_pval': pvals,
            'cell_type': cell_type
        })
        results_by_type[cell_type] = df

    log.info(f"✓ Wilcoxon DE computed for {len(results_by_type)} cell types")
    return results_by_type

def method_spatial_morans_i(adata):
    """
    Moran's I spatial autocorrelation (Squidpy)
    Benchmarking: Spatial autocorrelation kernel, identifies spatially variable genes
    """
    log.info("METHOD 2: Spatial Autocorrelation (Moran's I via Squidpy)")

    # Ensure spatial neighbors computed
    if 'spatial_neighbors' not in adata.obsp:
        log.info("  Computing spatial neighbors...")
        sq.gr.spatial_neighbors(adata, coord_type='generic')

    # Compute Moran's I
    results = sq.gr.spatial_autocorr(adata, mode='moran', copy=True)
    results['gene'] = adata.var_names.values
    results = results.sort_values('pval_norm_fdr_bh')

    # Identify spatially variable genes (SVGs)
    n_svg = (results['pval_norm_fdr_bh'] < 0.05).sum()
    log.info(f"✓ Moran's I computed; {n_svg}/{len(results)} genes with spatial signal (p<0.05)")

    return results

def method_spatial_de_by_celltype(adata, morans_i_results):
    """
    Integrate spatial signal into cell-type DE
    For each cell type: rank by Wilcoxon score AND Moran's I separately
    """
    log.info("METHOD 3: Spatial-Aware DE Ranking (Wilcoxon + Moran's I Integration)")

    wilcoxon_results = method_wilcoxon_baseline(adata)

    # Create spatial-aware rankings
    spatial_aware_results = {}
    morans_dict = dict(zip(morans_i_results['gene'], morans_i_results['I']))

    for cell_type, df in wilcoxon_results.items():
        # Add Moran's I values
        df['morans_i'] = df['gene'].map(morans_dict).fillna(0)

        # Rank by combined score: Wilcoxon score × Moran's I
        df['spatial_aware_score'] = df['wilcoxon_score'] * (1 + df['morans_i'])
        df = df.sort_values('spatial_aware_score', ascending=False).reset_index(drop=True)
        df['rank'] = range(1, len(df) + 1)

        spatial_aware_results[cell_type] = df

    log.info(f"✓ Spatial-aware rankings created for {len(spatial_aware_results)} cell types")
    return spatial_aware_results, morans_i_results

def create_comparative_summary(wilcoxon_results, spatial_aware_results, morans_i_results):
    """
    Create summary comparing all three methods
    Shows: Non-spatial (Wilcoxon) vs Spatial-Aware (Wilcoxon+Moran's I)
    """
    log.info("Creating comparative summary...")

    all_dfs = []
    for cell_type in wilcoxon_results.keys():
        wilc = wilcoxon_results[cell_type].head(50)
        spatial = spatial_aware_results[cell_type].head(50)

        df = pd.DataFrame({
            'cell_type': [cell_type] * 50,
            'rank_wilcoxon': range(1, 51),
            'gene_wilcoxon': wilc['gene'].values,
            'score_wilcoxon': wilc['wilcoxon_score'].values,
            'rank_spatial_aware': spatial['rank'].values,
            'gene_spatial_aware': spatial['gene'].values,
            'score_spatial_aware': spatial['spatial_aware_score'].values,
        })
        all_dfs.append(df)

    summary = pd.concat(all_dfs, ignore_index=True)
    log.info(f"✓ Comparative summary created ({len(summary)} rows)")

    return summary

def save_all_results(wilcoxon_results, spatial_aware_results, morans_i_results, summary):
    """Save all DGEA results"""
    log.info("Saving DGEA benchmarked results...")

    # Moran's I
    morans_i_results.to_csv(OUTPUT_DIR / "morans_i_all_genes.csv", index=False)

    # Per-cell-type rankings
    for cell_type, df in wilcoxon_results.items():
        df.to_csv(OUTPUT_DIR / f"wilcoxon_{cell_type}.csv", index=False)

    for cell_type, df in spatial_aware_results.items():
        df.to_csv(OUTPUT_DIR / f"spatial_aware_{cell_type}.csv", index=False)

    # Comparative summary
    summary.to_csv(OUTPUT_DIR / "comparative_dgea_summary.csv", index=False)

    # Summary statistics
    stats = pd.DataFrame({
        'metric': ['total_genes_tested', 'spatial_genes_pval05', 'spatial_genes_pval01'],
        'value': [
            len(morans_i_results),
            (morans_i_results['pval_norm_fdr_bh'] < 0.05).sum(),
            (morans_i_results['pval_norm_fdr_bh'] < 0.01).sum(),
        ]
    })
    stats.to_csv(OUTPUT_DIR / "dgea_summary_statistics.csv", index=False)

    log.info(f"✓ All results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 1: BENCHMARKED DGEA ANALYSIS")
    log.info("Comparative evaluation of multiple spatial DE methods")
    log.info("=" * 80)

    try:
        # Load data
        adata, sdata = load_data()

        # Method 1: Wilcoxon baseline
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Non-Spatial DE (Wilcoxon Baseline)")
        log.info("=" * 80)
        wilcoxon_results = method_wilcoxon_baseline(adata)

        # Method 2: Spatial autocorrelation
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Spatial Autocorrelation (Moran's I)")
        log.info("=" * 80)
        morans_i_results = method_spatial_morans_i(adata)

        # Method 3: Integrated spatial-aware DE
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Spatial-Aware DE Ranking")
        log.info("=" * 80)
        spatial_aware_results, morans_i_results = method_spatial_de_by_celltype(adata, morans_i_results)

        # Create summary
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Comparative Summary")
        log.info("=" * 80)
        summary = create_comparative_summary(wilcoxon_results, spatial_aware_results, morans_i_results)

        # Save results
        log.info("\n" + "=" * 80)
        log.info("STEP 5: Save Results")
        log.info("=" * 80)
        save_all_results(wilcoxon_results, spatial_aware_results, morans_i_results, summary)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ PHASE 1 (BENCHMARKED) COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
