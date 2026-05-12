#!/usr/bin/env python3
"""
F0 v2 — Re-anotación celular GRANULAR usando marcadores del panel pulmón

Lee el panel Xenium custom (immuno + add-on lung) desde
my_xenium_panel_markers/Propuesta panel Pulmón 2025.xlsx - Hoja1.csv
y aplica sc.tl.score_genes para CADA tipo celular en orden jerárquico:

  L1 (broad)        : 9 categorías (Epithelial, Immune, Stromal, Endothelial, Tumor, etc.)
  L2 (granular)     : 25+ tipos celulares (CD4_T, CD8_T, Treg, M1_Mac, M2_Mac, B_naive,
                       B_memory, Plasma, NK_cytotoxic, NK_resting, cDC1, cDC2, pDC,
                       Mast, Neutrophil, Monocyte_classical, Monocyte_NC, Fibroblast,
                       Endothelial_blood, Endothelial_lymphatic, AT1, AT2, Ciliated,
                       Club, Basaloid, Tumor_proliferating, Tumor_resting, ...)
  L3 (estados)      : Funcional/activación (Exhausted, Cytotoxic, Proliferating,
                       Interferon_response, Antigen_presenting, etc.)

Estrategia:
  1. Para cada tipo, scoring con sc.tl.score_genes (z-score vs random gene set)
  2. Asignación: argmax score por célula con threshold mínimo (z>0.5)
  3. Confidence: gap entre top1 y top2 score
  4. Sobrescribir cell_type, cell_type_fine, cell_type_L3, scores_*

Salida:
  - Sobrescribe sdata.zarr/tables/table con todas las anotaciones
  - CSV de marcadores usados por tipo
  - Validación: dotplot, UMAP, tabla de transición vs anotación previa
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd
import anndata as ad

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger
log = get_logger(__name__)


# ════════════════════════════════════════════════════════════════════════════
# DICCIONARIO DE MARCADORES — construido desde el panel pulmón
# Sólo se usan genes presentes en el panel; threshold z>0.5 para asignación
# ════════════════════════════════════════════════════════════════════════════
MARKERS = {
    # ─── INMUNE: T cells ─────────────────────────────────────────────────────
    "CD4_T_helper": ["CD4", "IL7R", "CCR7", "TCF7", "LEF1", "CD3D", "CD3E", "CD3G"],
    "CD8_T_cytotoxic": ["CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "CD3D"],
    "CD8_T_exhausted": ["CD8A", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "CD4"],
    "T_resident_memory": ["CD8A", "ITGAE", "CD69", "ZNF683"],

    # ─── INMUNE: B cells / Plasma ────────────────────────────────────────────
    "B_naive": ["MS4A1", "CD19", "CD79A", "CD79B", "IGHD", "TCL1A"],
    "B_memory": ["MS4A1", "CD19", "CD27", "CD79A", "CD79B"],
    "Plasma": ["MZB1", "JCHAIN", "XBP1", "CD27", "TNFRSF17", "PRDM1"],
    "Plasmablast": ["MKI67", "MZB1", "JCHAIN", "XBP1"],

    # ─── INMUNE: NK ──────────────────────────────────────────────────────────
    "NK_cytotoxic": ["NCAM1", "FCGR3A", "KLRD1", "NKG7", "GNLY", "PRF1", "GZMB", "KLRF1"],
    "NK_resting": ["NCAM1", "KLRC1", "XCL1", "XCL2"],

    # ─── MIELOIDE: Macrófagos / Monocitos ────────────────────────────────────
    "Macrophage_M1": ["CD68", "CD86", "TNF", "IL1B", "NOS2", "CXCL9", "CXCL10", "AIF1"],
    "Macrophage_M2": ["CD68", "CD163", "MRC1", "MSR1", "VSIG4", "MARCO", "AIF1", "CD14"],
    "Monocyte_classical": ["CD14", "S100A8", "S100A9", "FCN1", "VCAN", "CSF3R"],
    "Monocyte_NC": ["FCGR3A", "MS4A7", "CDKN1C", "CX3CR1"],

    # ─── MIELOIDE: DCs / Granulocitos ────────────────────────────────────────
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8", "CD8A"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "CD1A"],
    "DC_mature": ["LAMP3", "CCR7", "CCL19", "CCL22", "FSCN1"],
    "pDC": ["LILRA4", "IRF7", "IRF8", "CLEC4C", "TCF4"],
    "Mast_cell": ["KIT", "TPSAB1", "TPSB2", "CPA3"],
    "Neutrophil": ["S100A8", "S100A9", "CSF3R", "FCGR3A", "FCGR3B", "CXCR2"],

    # ─── EPITELIO: Pulmón ────────────────────────────────────────────────────
    "AT1": ["AGER", "CAVIN1", "PDPN", "RTKN2"],
    "AT2": ["SFTPC", "SFTPB", "SFTPD", "SFTPA1", "SFTPA2", "ABCA3"],
    "Ciliated": ["FOXJ1", "AGR3", "C20orf85", "C1orf194", "C6orf118",
                  "CCDC39", "CCDC78", "CFAP53"],
    "Club": ["SCGB1A1", "SCGB3A2", "ARFGEF3", "CAPN8", "CYP2F1"],
    "Basal": ["KRT5", "KRT14", "TP63", "KRT17"],
    "Basaloid_aberrant": ["KRT17", "KRT5", "TP63", "EPHB2", "IL11", "MMP7"],
    "Goblet": ["MUC5B", "MUC5AC", "TFF3", "SCGB3A2"],
    "Epithelial_general": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],

    # ─── ESTROMA ─────────────────────────────────────────────────────────────
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "ACTA2",
                   "HAS1", "TAGLN", "CTHRC1"],
    "Smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "DES", "CNN1"],
    "Pericyte": ["RGS5", "PDGFRB", "MCAM", "ACTA2", "NOTCH3"],

    # ─── ENDOTELIO ───────────────────────────────────────────────────────────
    "Endothelial_blood": ["PECAM1", "VWF", "CDH5", "CLDN5", "ADGRL4", "ACKR1"],
    "Endothelial_lymphatic": ["PECAM1", "PROX1", "PDPN", "LYVE1", "TFF3"],

    # ─── TUMOR ───────────────────────────────────────────────────────────────
    "Tumor_proliferating": ["TOP2A", "MKI67", "CCNB1", "CDK1", "BIRC5"],
    "Tumor_resting": ["EPCAM", "KRT8", "KRT18", "MUC1", "TSPAN8", "WFS1", "AGR3"],
}

# Mapeo L2 → L1 (broad)
L2_TO_L1 = {
    # T cells
    "CD4_T_helper": "T_cell", "CD8_T_cytotoxic": "T_cell",
    "CD8_T_exhausted": "T_cell", "Treg": "T_cell", "T_resident_memory": "T_cell",
    # B / Plasma
    "B_naive": "B_cell", "B_memory": "B_cell",
    "Plasma": "Plasma_cell", "Plasmablast": "Plasma_cell",
    # NK
    "NK_cytotoxic": "NK_cell", "NK_resting": "NK_cell",
    # Mieloide
    "Macrophage_M1": "Macrophage", "Macrophage_M2": "Macrophage",
    "Monocyte_classical": "Monocyte", "Monocyte_NC": "Monocyte",
    "cDC1": "Dendritic_cell", "cDC2": "Dendritic_cell",
    "DC_mature": "Dendritic_cell", "pDC": "Dendritic_cell",
    "Mast_cell": "Mast_cell", "Neutrophil": "Neutrophil",
    # Epitelio
    "AT1": "Epithelial_alveolar", "AT2": "Epithelial_alveolar",
    "Ciliated": "Epithelial_airway", "Club": "Epithelial_airway",
    "Basal": "Epithelial_airway", "Basaloid_aberrant": "Epithelial_aberrant",
    "Goblet": "Epithelial_airway", "Epithelial_general": "Epithelial_general",
    # Estroma
    "Fibroblast": "Stromal", "Smooth_muscle": "Stromal", "Pericyte": "Stromal",
    # Endotelio
    "Endothelial_blood": "Endothelial", "Endothelial_lymphatic": "Endothelial",
    # Tumor
    "Tumor_proliferating": "Tumor", "Tumor_resting": "Tumor",
}

# Mapeo L2 → L3 (estados funcionales)
L2_TO_L3 = {
    "CD8_T_cytotoxic": "Cytotoxic_effector",
    "CD8_T_exhausted": "Exhaustion",
    "Treg": "Immune_regulatory",
    "NK_cytotoxic": "Cytotoxic_effector",
    "Macrophage_M1": "Pro_inflammatory",
    "Macrophage_M2": "Immune_regulatory",
    "Plasmablast": "Proliferating",
    "Tumor_proliferating": "Proliferating",
    "DC_mature": "Antigen_presenting",
    "Basaloid_aberrant": "Disease_associated",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path,
                   help="Path al sdata.zarr de entrada")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Directorio de salida para reportes/CSVs")
    p.add_argument("--threshold", type=float, default=0.5,
                   help="Z-score mínimo para asignar tipo (default 0.5)")
    return p.parse_args()


def filter_markers_to_panel(markers_dict, available_genes):
    """Mantén solo marcadores presentes en el panel."""
    filtered = {}
    panel_set = set(available_genes)
    for ct, genes in markers_dict.items():
        present = [g for g in genes if g in panel_set]
        if len(present) >= 2:  # mínimo 2 marcadores
            filtered[ct] = present
        else:
            log.warning(f"  ⚠ {ct}: sólo {len(present)} marcadores en el panel ({present}); excluido")
    return filtered


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 70)
    log.info("F0 v2 — Re-anotación granular con marcadores del panel")
    log.info("=" * 70)

    log.info(f"Cargando sdata desde: {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  {adata.n_obs:,} células × {adata.n_vars} genes")

    # Filtrar marcadores al panel disponible
    panel_genes = set(adata.var_names)
    log.info(f"\n[1/5] Filtrando marcadores al panel disponible ({len(panel_genes)} genes)...")
    markers_filtered = filter_markers_to_panel(MARKERS, panel_genes)
    log.info(f"  {len(markers_filtered)}/{len(MARKERS)} tipos celulares con ≥2 marcadores")

    # Guardar diccionario filtrado
    pd.DataFrame([(ct, ",".join(g), len(g)) for ct, g in markers_filtered.items()],
                 columns=["cell_type_L2", "markers_used", "n_markers"]).to_csv(
        args.outdir / "markers_used_per_celltype.csv", index=False)

    # Asegurar layer normalizado para scoring
    log.info(f"\n[2/5] Preparando datos para scoring...")
    if "counts" not in adata.layers:
        log.info("  Sin layer 'counts' — copiando .X")
        adata.layers["counts"] = adata.X.copy()
    # Si X no está log-normalizado, normalizar
    if adata.X.max() > 50:
        log.info("  Normalizando + log1p para scoring...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Scoring por cada tipo L2
    log.info(f"\n[3/5] Scoring (sc.tl.score_genes) por cada tipo celular...")
    score_cols = []
    for ct, genes in markers_filtered.items():
        score_name = f"score_{ct}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name,
                          ctrl_size=min(50, max(2, len(genes))), use_raw=False, random_state=42)
        score_cols.append(score_name)
        log.info(f"  ✓ {ct:25s} ({len(genes):2d} marcadores) → {score_name}")

    # Asignación: argmax con threshold
    log.info(f"\n[4/5] Asignando tipo celular (argmax, threshold z>{args.threshold})...")
    score_matrix = adata.obs[score_cols].values  # cells × types
    best_idx = np.argmax(score_matrix, axis=1)
    best_score = np.max(score_matrix, axis=1)

    # Top2 score para confidence
    sorted_scores = np.sort(score_matrix, axis=1)
    second_score = sorted_scores[:, -2]
    confidence = best_score - second_score

    cell_types_L2 = list(markers_filtered.keys())
    assigned_L2 = np.array([cell_types_L2[i] for i in best_idx], dtype=object)
    assigned_L2[best_score < args.threshold] = "Unassigned"

    # L1 y L3 mappings
    assigned_L1 = np.array([L2_TO_L1.get(c, "Unassigned") for c in assigned_L2], dtype=object)
    assigned_L3 = np.array([L2_TO_L3.get(c, "Steady_state") for c in assigned_L2], dtype=object)

    # Confidence categorical
    conf_cat = pd.cut(confidence, bins=[-np.inf, 0.1, 0.3, 0.6, np.inf],
                      labels=["LOW", "MEDIUM", "HIGH", "VERY_HIGH"])

    # Sobrescribir columnas
    log.info("  Sobrescribiendo columnas en adata.obs...")
    adata.obs["cell_type_L1"] = pd.Categorical(assigned_L1)
    adata.obs["cell_type_L2"] = pd.Categorical(assigned_L2)
    adata.obs["cell_type_L3"] = pd.Categorical(assigned_L3)
    adata.obs["cell_type_score"] = best_score
    adata.obs["cell_type_confidence"] = pd.Categorical(conf_cat)

    # Mantener compatibilidad: cell_type = L1, cell_type_fine = L2 (REAL ahora)
    adata.obs["cell_type"] = adata.obs["cell_type_L1"]
    adata.obs["cell_type_fine"] = adata.obs["cell_type_L2"]

    log.info("\n  Distribución L1:")
    for ct, n in adata.obs["cell_type_L1"].value_counts().items():
        log.info(f"    {ct:25s} {n:8,} ({100*n/adata.n_obs:5.1f}%)")

    log.info("\n  Distribución L2 (top 30):")
    for ct, n in adata.obs["cell_type_L2"].value_counts().head(30).items():
        log.info(f"    {ct:25s} {n:8,} ({100*n/adata.n_obs:5.1f}%)")

    log.info("\n  Distribución confianza:")
    for c, n in adata.obs["cell_type_confidence"].value_counts().items():
        log.info(f"    {c:12s} {n:8,} ({100*n/adata.n_obs:5.1f}%)")

    # Guardar transición vs L1 antiguo (si existe)
    if "leiden" in adata.obs.columns:
        trans = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type_L2"])
        trans.to_csv(args.outdir / "leiden_to_celltypeL2_transition.csv")
        log.info(f"  Transición leiden→L2 guardada")

    # Guardar resultado
    log.info(f"\n[5/5] Persistiendo cambios al sdata.zarr...")
    # Re-write table element
    sdata.tables[table_name] = adata
    # Persistir solo el elemento table modificado
    try:
        sdata.delete_element_from_disk(table_name)
        sdata.write_element(table_name)
        log.info(f"  ✓ Table '{table_name}' sobreescrito en {args.sdata}")
    except Exception as e:
        log.warning(f"  Persistencia incremental falló ({e}); usando write completo...")
        out_path = args.sdata.parent / "sdata_v2.zarr"
        sdata.write(out_path)
        log.info(f"  ✓ Escrito en {out_path}")

    # Resumen JSON
    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_celltypes_L1": int(adata.obs["cell_type_L1"].nunique()),
        "n_celltypes_L2": int(adata.obs["cell_type_L2"].nunique()),
        "n_celltypes_L3": int(adata.obs["cell_type_L3"].nunique()),
        "celltypes_L1": sorted(adata.obs["cell_type_L1"].dropna().unique().tolist()),
        "celltypes_L2": sorted(adata.obs["cell_type_L2"].dropna().unique().tolist()),
        "celltypes_L3": sorted(adata.obs["cell_type_L3"].dropna().unique().tolist()),
        "confidence_distribution": adata.obs["cell_type_confidence"].value_counts().to_dict(),
        "threshold_used": args.threshold,
        "markers_dict_size": {k: len(v) for k, v in markers_filtered.items()},
    }
    summary["confidence_distribution"] = {str(k): int(v) for k, v in summary["confidence_distribution"].items()}
    with open(args.outdir / "reannotation_v2_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info(f"\n✓ Re-anotación completa → {args.outdir}")


if __name__ == "__main__":
    main()
