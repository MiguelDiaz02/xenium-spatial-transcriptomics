#!/usr/bin/env python3
"""
Week 3 Phase 2A: Cell-Cell Communication (CCC) Analysis
Ligand-Receptor interactions with focus on Immune-Tumor and M2 macrophage hub validation
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

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
PHASE1_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Known L/R pairs from CellChatDB + FANTOM5 (curated for lung cancer immune context)
# Format: (ligand, receptor, source_cell_type, target_cell_type_optional)
KNOWN_LR_PAIRS = [
    # Macrophage M1/M2 interactions
    ("CD68", "CD14", "Macrophage", None),
    ("TNF", "TNFRSF1A", "Macrophage", None),
    ("IL12A", "IL12RB1", "Macrophage", "T_cell"),
    ("IL10", "IL10RA", "Macrophage", None),
    ("TGM2", "ITGAE", "Macrophage", "T_cell"),

    # T cell activation / suppression
    ("CD86", "CTLA4", "Macrophage", "T_cell"),
    ("CD80", "CTLA4", "Macrophage", "T_cell"),
    ("PDCD1LG2", "PDCD1", "Macrophage", "T_cell"),
    ("IDO1", "AHR", "Macrophage", "T_cell"),
    ("ICAM1", "LFA1", "Endothelial", "T_cell"),

    # NK cell interactions
    ("HLA-A", "KIR2DL1", "Tumor", "NK_cell"),
    ("HLA-A", "NKG2D", "Tumor", "NK_cell"),
    ("ULBP1", "NKG2D", "Tumor", "NK_cell"),

    # B cell / Plasma interactions
    ("CD40", "CD40LG", "B_cell", "T_cell"),
    ("CD86", "CD28", "B_cell", "T_cell"),

    # Epithelial-Immune
    ("ICAM1", "ITGAL", "Epithelial", "T_cell"),
    ("EPCAM", "VTCN1", "Epithelial", "T_cell"),

    # Cytotoxic granules
    ("PRF1", "LYST", None, None),
    ("GZMA", "SLC7A2", None, None),
    ("GZMB", "SLC7A2", None, None),
]

def load_data():
    """Load SpatialData and Phase 1 results"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Load Phase 1 results
    log.info("Loading Phase 1A DGE results...")
    phase1_summary = pd.read_csv(PHASE1_DIR / "DGE_summary_all_types.csv")
    spatial_stats = pd.read_csv(PHASE1_DIR / "spatial_autocorr_morans_i.csv")
    log.info(f"✓ Loaded Phase 1 results: {len(phase1_summary)} gene-celltype pairs")

    return adata, phase1_summary, spatial_stats

def extract_lr_pairs_from_panel(adata, known_pairs):
    """
    Extract L/R pairs that exist in the gene panel

    Args:
        adata: AnnData object
        known_pairs: List of (ligand, receptor, source, target) tuples

    Returns:
        DataFrame with available L/R pairs
    """
    log.info(f"Extracting L/R pairs from panel ({len(adata.var_names)} genes)...")

    panel_genes = set(adata.var_names)
    available_pairs = []

    for ligand, receptor, source, target in known_pairs:
        # Check if both genes exist
        if ligand in panel_genes and receptor in panel_genes:
            available_pairs.append({
                'ligand': ligand,
                'receptor': receptor,
                'source': source,
                'target': target,
                'pair_name': f"{ligand}→{receptor}"
            })

    lr_df = pd.DataFrame(available_pairs)
    log.info(f"✓ Found {len(lr_df)} L/R pairs in panel")

    return lr_df

def compute_lr_interactions(adata, lr_pairs, phase1_summary):
    """
    Compute ligand-receptor interaction scores
    Score = log2(mean_expr_ligand + mean_expr_receptor + 1)

    Args:
        adata: AnnData object
        lr_pairs: DataFrame with L/R pairs
        phase1_summary: Phase 1 DE results

    Returns:
        DataFrame with interaction scores per cell type pair
    """
    log.info("Computing L/R interaction scores...")

    # Get normalized expression matrix
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    cell_types = adata.obs['cell_type'].unique()
    results = []

    for _, row in lr_pairs.iterrows():
        ligand, receptor = row['ligand'], row['receptor']

        # Get gene indices
        try:
            ligand_idx = adata.var_names.get_loc(ligand)
            receptor_idx = adata.var_names.get_loc(receptor)
        except KeyError:
            continue

        # Compute mean expression per cell type
        expr_by_ct = {}
        for ct in cell_types:
            mask = adata.obs['cell_type'] == ct
            if mask.sum() > 0:
                ligand_mean = X[mask, ligand_idx].mean()
                receptor_mean = X[mask, receptor_idx].mean()
                expr_by_ct[ct] = {
                    'ligand_expr': ligand_mean,
                    'receptor_expr': receptor_mean,
                }

        # Create all pairwise interactions
        for source in cell_types:
            for target in cell_types:
                if source in expr_by_ct and target in expr_by_ct:
                    # L/R interaction score
                    lig_expr = expr_by_ct[source]['ligand_expr']
                    rec_expr = expr_by_ct[target]['receptor_expr']
                    score = np.log2(lig_expr + rec_expr + 1)

                    # DE enrichment from Phase 1
                    source_de = phase1_summary[
                        (phase1_summary['gene'] == ligand) &
                        (phase1_summary['cell_type'] == source)
                    ]['score'].values

                    target_de = phase1_summary[
                        (phase1_summary['gene'] == receptor) &
                        (phase1_summary['cell_type'] == target)
                    ]['score'].values

                    source_de = source_de[0] if len(source_de) > 0 else 0
                    target_de = target_de[0] if len(target_de) > 0 else 0

                    results.append({
                        'ligand': ligand,
                        'receptor': receptor,
                        'source_cell_type': source,
                        'target_cell_type': target,
                        'pair_name': row['pair_name'],
                        'lr_score': score,
                        'ligand_expr': lig_expr,
                        'receptor_expr': rec_expr,
                        'source_de': source_de,
                        'target_de': target_de,
                    })

    interactions_df = pd.DataFrame(results)

    # Filter: only keep interactions with meaningful score
    interactions_df = interactions_df[interactions_df['lr_score'] > 0.5]

    log.info(f"✓ Computed {len(interactions_df)} source-target interactions")

    return interactions_df

def identify_m2_hub(interactions_df):
    """
    Identify M2 macrophage as communication hub
    Hub = receives most incoming signals

    Args:
        interactions_df: DataFrame with all interactions

    Returns:
        DataFrame with hub analysis
    """
    log.info("Analyzing M2 macrophage as communication hub...")

    # Incoming signals to M2
    m2_targets = interactions_df[interactions_df['target_cell_type'] == 'Macrophage'].copy()
    incoming_by_source = m2_targets.groupby('source_cell_type')['lr_score'].agg(['sum', 'count', 'mean']).reset_index()
    incoming_by_source = incoming_by_source.sort_values('sum', ascending=False)
    incoming_by_source.columns = ['source_cell_type', 'total_score', 'n_interactions', 'mean_score']

    # Outgoing signals from M2
    m2_sources = interactions_df[interactions_df['source_cell_type'] == 'Macrophage'].copy()
    outgoing_by_target = m2_sources.groupby('target_cell_type')['lr_score'].agg(['sum', 'count', 'mean']).reset_index()
    outgoing_by_target = outgoing_by_target.sort_values('sum', ascending=False)
    outgoing_by_target.columns = ['target_cell_type', 'total_score', 'n_interactions', 'mean_score']

    log.info("M2 Macrophage Hub Analysis:")
    log.info(f"  Incoming signals: {len(incoming_by_source)} cell types")
    log.info(f"  Outgoing signals: {len(outgoing_by_target)} cell types")

    # Identify key pathways
    m2_pathways = m2_targets.nlargest(10, 'lr_score')[['pair_name', 'source_cell_type', 'lr_score']]

    return {
        'incoming': incoming_by_source,
        'outgoing': outgoing_by_target,
        'top_pathways': m2_pathways
    }

def validate_immune_tumor_interactions(interactions_df):
    """
    Focus on immune-tumor interactions
    Key question: Which immune populations engage tumor?

    Args:
        interactions_df: DataFrame with all interactions

    Returns:
        DataFrame with immune-tumor interactions ranked
    """
    log.info("Validating immune-tumor interactions...")

    immune_types = ['T_cell', 'B_cell', 'NK_cell', 'Macrophage', 'Endothelial']

    # Filter for immune → tumor interactions
    immune_tumor = interactions_df[
        (interactions_df['source_cell_type'].isin(immune_types)) &
        (interactions_df['target_cell_type'] == 'Tumor')
    ].copy()

    # Or tumor → immune
    tumor_immune = interactions_df[
        (interactions_df['source_cell_type'] == 'Tumor') &
        (interactions_df['target_cell_type'].isin(immune_types))
    ].copy()

    # Rank by interaction strength
    immune_tumor = immune_tumor.sort_values('lr_score', ascending=False)
    tumor_immune = tumor_immune.sort_values('lr_score', ascending=False)

    log.info(f"✓ Immune→Tumor interactions: {len(immune_tumor)}")
    log.info(f"✓ Tumor→Immune interactions: {len(tumor_immune)}")

    return {
        'immune_to_tumor': immune_tumor,
        'tumor_to_immune': tumor_immune
    }

def save_results(interactions_df, lr_pairs, m2_hub, immune_tumor):
    """Save all CCC analysis results to CSV"""
    log.info("Saving CCC analysis results...")

    # Main interaction table
    interactions_df.to_csv(OUTPUT_DIR / "lr_interactions_all.csv", index=False)

    # L/R pairs summary
    lr_pairs.to_csv(OUTPUT_DIR / "lr_pairs_available.csv", index=False)

    # M2 hub analysis
    m2_hub['incoming'].to_csv(OUTPUT_DIR / "m2_incoming_signals.csv", index=False)
    m2_hub['outgoing'].to_csv(OUTPUT_DIR / "m2_outgoing_signals.csv", index=False)
    m2_hub['top_pathways'].to_csv(OUTPUT_DIR / "m2_top_pathways.csv", index=False)

    # Immune-tumor interactions
    immune_tumor['immune_to_tumor'].to_csv(OUTPUT_DIR / "immune_to_tumor_interactions.csv", index=False)
    immune_tumor['tumor_to_immune'].to_csv(OUTPUT_DIR / "tumor_to_immune_interactions.csv", index=False)

    log.info(f"✓ Results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 2A: CELL-CELL COMMUNICATION (CCC) ANALYSIS")
    log.info("Ligand-Receptor Interactions with M2 Macrophage Hub Validation")
    log.info("=" * 80)

    try:
        # Load data
        adata, phase1_summary, spatial_stats = load_data()

        # Extract L/R pairs from panel
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Extract L/R Pairs Available in Panel")
        log.info("=" * 80)
        lr_pairs = extract_lr_pairs_from_panel(adata, KNOWN_LR_PAIRS)

        # Compute interactions
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Compute L/R Interaction Scores")
        log.info("=" * 80)
        interactions_df = compute_lr_interactions(adata, lr_pairs, phase1_summary)

        # M2 hub analysis
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Identify M2 Macrophage as Communication Hub")
        log.info("=" * 80)
        m2_hub = identify_m2_hub(interactions_df)

        # Immune-tumor interactions
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Validate Immune-Tumor Interactions")
        log.info("=" * 80)
        immune_tumor = validate_immune_tumor_interactions(interactions_df)

        # Save results
        log.info("\n" + "=" * 80)
        log.info("STEP 5: Save Results")
        log.info("=" * 80)
        save_results(interactions_df, lr_pairs, m2_hub, immune_tumor)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ PHASE 2A COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
