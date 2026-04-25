#!/usr/bin/env python3
"""
Week 3 Phase 2B: DeepTalk-Inspired CCC Analysis (GNN-based)
Graph Neural Network approach to ligand-receptor interactions with spatial context

Based on benchmarking recommendation:
- DeepTalk (Yang et al., Nature Comm 2024) — state-of-art GNN approach
- Incorporates: spatial distance, multi-sender dynamics, single-cell resolution
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.stats import zscore

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
import torch
import torch.nn as nn
from torch_geometric.data import Data
from torch_geometric.nn import GATConv, MessagePassing
from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
PHASE1_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_deeptalk_gnn"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Known L/R pairs from CellChatDB (curated list for lung cancer)
KNOWN_LR_PAIRS = [
    ("CD68", "CD14"), ("TNF", "TNFRSF1A"), ("IL12A", "IL12RB1"),
    ("IL10", "IL10RA"), ("CD86", "CTLA4"), ("CD80", "CTLA4"),
    ("PDCD1LG2", "PDCD1"), ("IDO1", "AHR"), ("ICAM1", "LFA1"),
    ("HLA-A", "KIR2DL1"), ("ULBP1", "NKG2D"), ("CD40", "CD40LG"),
    ("PRF1", "LYST"), ("GZMA", "SLC7A2"), ("GZMB", "SLC7A2"),
]

class GAT_CCC(nn.Module):
    """Graph Attention Network for CCC inference (DeepTalk-inspired)"""
    def __init__(self, in_channels, hidden_channels, out_channels, heads=4):
        super().__init__()
        self.gat1 = GATConv(in_channels, hidden_channels, heads=heads, dropout=0.2)
        self.gat2 = GATConv(hidden_channels * heads, out_channels, heads=1, dropout=0.2)

    def forward(self, x, edge_index, edge_weight=None):
        x = self.gat1(x, edge_index, edge_attr=edge_weight)
        x = torch.relu(x)
        x = self.gat2(x, edge_index, edge_attr=edge_weight)
        return x

def load_data():
    """Load SpatialData and Phase 1 results"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Load Phase 1 results
    phase1_summary = pd.read_csv(PHASE1_DIR / "DGE_summary_all_types.csv")

    return adata, phase1_summary

def build_spatial_graph(adata, k_neighbors=15):
    """
    Build spatial proximity graph using k-nearest neighbors
    DeepTalk key component: spatial context as graph structure
    """
    log.info(f"Building spatial k-NN graph (k={k_neighbors})...")

    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        coords = adata.obs[['x_coordinate', 'y_coordinate']].values

    # Compute k-NN graph
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=k_neighbors+1).fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # Remove self-loops and create edge list
    edge_list = []
    edge_weights = []
    for i in range(len(indices)):
        for j, dist in zip(indices[i][1:], distances[i][1:]):  # Skip self (index 0)
            edge_list.append([i, j])
            # Weight inversely by distance (closer = stronger)
            edge_weights.append(np.exp(-dist / np.mean(distances)))

    edge_index = torch.LongTensor(edge_list).t().contiguous()
    edge_weight = torch.FloatTensor(edge_weights)

    log.info(f"✓ Built graph with {len(edge_list)} edges")

    return edge_index, edge_weight, coords

def extract_lr_features(adata, lr_pairs):
    """
    Extract L/R pair expression features for each cell
    Input to GNN: expression vectors for ligands + receptors
    """
    log.info(f"Extracting L/R pair expression features...")

    from scipy.sparse import issparse
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    panel_genes = set(adata.var_names)
    available_pairs = []
    features_list = []

    for ligand, receptor in lr_pairs:
        if ligand in panel_genes and receptor in panel_genes:
            lig_idx = adata.var_names.get_loc(ligand)
            rec_idx = adata.var_names.get_loc(receptor)

            # Create feature vector: [ligand_expr, receptor_expr, cell_type_onehot]
            lig_expr = X[:, lig_idx:lig_idx+1]  # (n_cells, 1)
            rec_expr = X[:, rec_idx:rec_idx+1]  # (n_cells, 1)

            # Cell type encoding
            cell_types = pd.Categorical(adata.obs['cell_type'])
            n_types = len(cell_types.categories)
            ct_onehot = np.eye(n_types)[cell_types.codes]  # (n_cells, n_types)

            # Combine features
            pair_features = np.hstack([lig_expr, rec_expr, ct_onehot])
            features_list.append(pair_features)
            available_pairs.append(f"{ligand}→{receptor}")

    # Stack all pairs: (n_cells, n_features * n_pairs)
    all_features = np.hstack(features_list)

    log.info(f"✓ Extracted {len(available_pairs)} L/R pairs; feature matrix: {all_features.shape}")

    return torch.FloatTensor(all_features), available_pairs

def infer_ccc_with_gnn(adata, edge_index, edge_weight, features, available_pairs):
    """
    Infer CCC using Graph Attention Network (DeepTalk-inspired)
    Learns spatial context + L/R expression patterns
    """
    log.info("Inferring CCC with GNN (DeepTalk-inspired GAT)...")

    # Create PyG Data object
    data = Data(x=features, edge_index=edge_index, edge_attr=edge_weight)

    # Initialize model (small for efficiency)
    model = GAT_CCC(
        in_channels=features.shape[1],
        hidden_channels=64,
        out_channels=32,
        heads=4
    )

    # Simple training loop (unsupervised, learns spatial patterns)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    model.train()

    log.info("  Training GNN (100 iterations)...")
    for epoch in range(100):
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.edge_attr)

        # Loss: reconstruct L/R expression patterns
        loss = torch.nn.functional.mse_loss(out, features[:, :out.shape[1]])
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 25 == 0:
            log.info(f"    Epoch {epoch+1}: loss={loss.item():.4f}")

    # Extract learned representations
    model.eval()
    with torch.no_grad():
        ccc_scores = model(data.x, data.edge_index, data.edge_attr).numpy()

    log.info("✓ GNN training complete")

    return ccc_scores

def compute_celltype_ccc_matrix(adata, ccc_scores, available_pairs):
    """
    Aggregate single-cell CCC scores to cell type level
    DeepTalk innovation: maintains single-cell resolution while providing summary
    """
    log.info("Computing cell type-level CCC matrix...")

    cell_types = adata.obs['cell_type'].unique()
    n_types = len(cell_types)
    n_pairs = len(available_pairs)

    # Initialize matrix
    ccc_matrix = np.zeros((n_types, n_types, n_pairs))

    # Aggregate scores per cell type pair
    for i, ct_source in enumerate(cell_types):
        source_mask = adata.obs['cell_type'] == ct_source
        source_scores = ccc_scores[source_mask]

        for j, ct_target in enumerate(cell_types):
            target_mask = adata.obs['cell_type'] == ct_target

            # Use spatial context: cells near target
            if 'spatial_neighbors' not in adata.obsp:
                sq.gr.spatial_neighbors(adata, coord_type='generic')

            neighbors = adata.obsp['spatial_neighbors']
            target_neighbors = neighbors[target_mask].sum(axis=0).A1 > 0

            # Score = mean of source cells' L/R scores for cells near target
            source_near_target = source_mask.values & target_neighbors
            if source_near_target.sum() > 0:
                ccc_matrix[i, j, :] = ccc_scores[source_near_target].mean(axis=0)

    return ccc_matrix

def save_results(ccc_scores, ccc_matrix, available_pairs, adata):
    """Save GNN-inferred CCC results"""
    log.info("Saving GNN-inferred CCC results...")

    # Single-cell CCC scores
    ccc_df = pd.DataFrame(ccc_scores, columns=[f"L/R_score_{i}" for i in range(ccc_scores.shape[1])])
    ccc_df['cell_id'] = adata.obs_names
    ccc_df['cell_type'] = adata.obs['cell_type'].values
    ccc_df.to_csv(OUTPUT_DIR / "single_cell_ccc_scores.csv", index=False)

    # Cell type-level CCC summary
    cell_types = adata.obs['cell_type'].unique()
    for pair_idx, pair_name in enumerate(available_pairs):
        pair_matrix = ccc_matrix[:, :, pair_idx]
        pair_df = pd.DataFrame(pair_matrix, index=cell_types, columns=cell_types)
        pair_df.to_csv(OUTPUT_DIR / f"ccc_matrix_{pair_name.replace('→', '_to_')}.csv")

    # Overall CCC strength (aggregated across all L/R pairs)
    overall_ccc = ccc_matrix.mean(axis=2)
    overall_df = pd.DataFrame(overall_ccc, index=cell_types, columns=cell_types)
    overall_df.to_csv(OUTPUT_DIR / "ccc_matrix_overall.csv")

    log.info(f"✓ Results saved to {OUTPUT_DIR}")

def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 2B: DEEPTALK-INSPIRED CCC ANALYSIS (GNN-BASED)")
    log.info("Graph Attention Networks for spatial ligand-receptor inference")
    log.info("=" * 80)

    try:
        # Load data
        adata, phase1_summary = load_data()

        # Build spatial graph
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Build Spatial Proximity Graph")
        log.info("=" * 80)
        edge_index, edge_weight, coords = build_spatial_graph(adata, k_neighbors=15)

        # Extract L/R features
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Extract L/R Pair Features")
        log.info("=" * 80)
        features, available_pairs = extract_lr_features(adata, KNOWN_LR_PAIRS)

        # Infer CCC with GNN
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Infer CCC with Graph Attention Network")
        log.info("=" * 80)
        ccc_scores = infer_ccc_with_gnn(adata, edge_index, edge_weight, features, available_pairs)

        # Compute cell type-level matrix
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Aggregate to Cell Type Level")
        log.info("=" * 80)
        ccc_matrix = compute_celltype_ccc_matrix(adata, ccc_scores, available_pairs)

        # Save results
        log.info("\n" + "=" * 80)
        log.info("STEP 5: Save Results")
        log.info("=" * 80)
        save_results(ccc_scores, ccc_matrix, available_pairs, adata)

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
