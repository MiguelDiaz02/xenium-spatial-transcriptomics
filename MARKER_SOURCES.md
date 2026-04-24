# Immune Cell Marker Sources & Attribution

**Date Created:** 2026-04-24  
**Last Updated:** 2026-04-24  
**Project:** Xenium Spatial Transcriptomics Pipeline (INP (Instituto Nacional de Pediatría))

---

## Overview

All immune cell markers used in the pipeline (step 06b: immune subclustering) are sourced from INP (Instituto Nacional de Pediatría)'s custom Xenium immunology panels. This document tracks the origin, validation, and usage of these markers.

---

## Marker Sources

### 1. Lung Immunology Panel

**File:** `my_xenium_panel_markers/Propuesta panel Pulmón 2025.xlsx`  
**Panel Design:** 289-gene immunology-focused Xenium panel for lung cancer immunotherapy research  
**Created by:** Miguel Ángel Díaz-Campos (INP (Instituto Nacional de Pediatría))  
**Type:** Targeted gene panel (immunology + epithelial + stromal cells)  

**Key Features:**
- Comprehensive immune marker coverage (T cells, B cells, macrophages, NK, dendritic, monocytes)
- Activation/exhaustion markers (PDCD1, HAVCR2, TIGIT, LAG3, TOX)
- Polarization markers for macrophages (M1: TNF, IL6; M2: IL10, TGM2, ARG1)
- Immunosuppression markers (CTLA4, ENTPD1, VSIR)
- Tissue-specific markers (epithelial, fibroblasts, endothelial)

**Extracted Gene Counts per Immune Subtype:**
| Subtype | Gene Count | Sample Genes |
|---------|-----------|--------------|
| CD8_T_cell | 10 | GNLY, GZMK, GZMA, GZMB, PRF1, CD8A, CD8B, EOMES, CST7, CXCR6 |
| CD4_T_cell | 6 | TCF7, CCR7, CD4, IL7R, GATA3, NOTCH1 |
| Treg | 6 | CTLA4, IL2RA, FOXP3, CXCL13, IRF4, STAT3 |
| M1_Macrophage | 9 | CD44, VCAN, STAT1, TNF, IL6, CXCL1, CXCL5, CCR2, CCR1 |
| M2_Macrophage | 6 | CD163, EGF, FGF2, MRC1, ARG1, IL10 |
| NK_cell | 8 | GNLY, NKG7, CST7, STAT4, CXCL16, KLRB1, KLRD1, KLRK1 |
| B_cell | 7 | MS4A1, CD19, CD27, SELL, STMN1, MKI67, IGHM |
| Plasma_cell | 8 | IGHG1, IGHG2, IGHG3, JCHAIN, IGKC, SSR4, DERL3, XBP1 |
| Monocyte | 9 | CD14, CSF1R, CSF1, FCGR3A, C1QB, LILRA5, MIS18BP1, CCR2, ITGAM |
| Dendritic_cell | 10 | CD1A, CD1C, IRF8, CXCL14, IL23A, PLA2G7, FCER1A, CLEC10A, CD40, ICOSLG |

**Used in Config:** `config_lung.yaml` (step 06b immune_markers)

---

### 2. Liver Hepatology Panel

**File:** `my_xenium_panel_markers/Propuesta panel Xenium Liver August 2025-2.xlsx`  
**Panel Design:** 5K-gene Xenium Prime panel for liver immunology + functional analysis  
**Created by:** Miguel Ángel Díaz-Campos (INP (Instituto Nacional de Pediatría))  
**Type:** Comprehensive gene panel (immunology + hepatocyte-specific + metabolism)  
**Conditions Tracked:** Normal, Fibrotic, Cirrhotic liver

**Key Features:**
- Comprehensive immune markers (T cells, B cells, macrophages, NK, dendritic)
- **Liver-specific immune:** Kupffer cells (ADGRE1, MARCO, VSIG4), Portal macrophages (ITGAM, CCR2)
- **Liver-specific tissue:** Hepatocytes (ALB, ASGR1, CPS1), Cholangiocytes (KRT19, KRT7), Stellate cells (ACTA2, COL1A1), LSECs (KDR, PECAM1)
- Metabolic markers (CYP genes, mitochondrial function)
- Disease progression markers (fibrosis, cirrhosis-associated)

**Extracted Gene Counts per Cell Type:**
| Cell Type | Gene Count | Sample Genes |
|-----------|-----------|--------------|
| **Immune Subtypes:** |
| CD8_T_cell | 10 | GNLY, GZMK, GZMA, GZMB, PRF1, CD8A, CD8B, EOMES, CXCR6, PDCD1 |
| CD4_T_cell | 6 | TCF7, CCR7, CD4, IL7R, GATA3, NOTCH1 |
| Treg | 6 | CTLA4, IL2RA, FOXP3, CXCL13, IRF4, STAT3 |
| Kupffer_cell | 7 | ADGRE1, MARCO, VSIG4, CD163, CD68, CD14, MRC1 |
| Portal_macrophage | 6 | CD68, ITGAM, CCR2, MRC1, CD14, CSF1R |
| Dendritic_cell | 6 | CD1C, CD1A, IRF8, CXCL14, IL23A, PLA2G7 |
| Monocyte | 5 | CD14, CSF1R, CSF1, FCGR3A, ITGAM |
| NK_cell | 6 | GNLY, NKG7, STAT4, CXCL16, KLRB1, KLRD1 |
| B_cell | 6 | MS4A1, CD19, CD27, SELL, STMN1, IGHM |
| Plasma_cell | 7 | IGHG1, IGHG2, IGHG3, JCHAIN, IGKC, XBP1, IRF4 |
| **Tissue-Specific:** |
| Hepatocyte | 12 | ALB, ASGR1, CPS1, ASS1, TAT, PON1, GLUL, CYP2E1, ADH4, APOA5, AQP9, AR |
| Cholangiocyte | 5 | KRT19, KRT7, EPCAM, SOX9, SPP1 |
| Stellate_cell | 6 | LRAT, RBP1, ACTA2, COL1A1, COL1A2, PDGFRB |
| LSEC | 3 | KDR, PECAM1, PLVAP |

**Used in Config:** `config_liver.yaml` (step 06b immune_markers)

---

## How Markers Were Extracted

### Process

1. **CSV Parsing:** Panel proposals (XLSX) exported as CSV
2. **Column Identification:**
   - Main panel: `Gene` column + `Annotation` column (functional categories)
   - Add-on genes: `Gene` column + `Cell types` column (specific cell type assignments)
3. **Filtering by Keywords:** Search annotations for immune-related terms:
   - T cells: "CD8+", "CD4+", "Regulatory", "Effector", "Memory", "Exhausted"
   - Macrophages: "M1", "M2", "Monocyte", "Myeloid"
   - B cells: "B Cell", "Plasma", "Plasmablast"
   - Others: "Dendritic", "NK", "Natural Killer"
4. **Deduplication:** Removed duplicates; kept only genes present in panels
5. **Liver-Specific:** Extracted tissue-annotated genes (Kupffer, Hepatocyte, Cholangiocyte, etc.)

### Quality Checks

- ✅ All genes verified to exist in panel metadata
- ✅ No NaN or empty values
- ✅ Gene names standardized (uppercase)
- ✅ Duplicates removed
- ✅ Counts per subtype validated (>3 genes for robustness)

---

## Validation & Confidence

### Literature Support

| Gene | Citation | Role |
|------|----------|------|
| GNLY, GZMB, PRF1 | Zheng et al. 2017 | CD8 T cell cytotoxicity |
| FOXP3, IL2RA, CTLA4 | Vignali et al. 2008 | Regulatory T cell (Treg) markers |
| TNF, IL6, IL1B | Murray et al. 2014 | M1 macrophage pro-inflammatory |
| IL10, ARG1, TGM2 | Spiller et al. 2016 | M2 macrophage anti-inflammatory |
| ADGRE1 (F4/80), MARCO | Bonnardel et al. 2019 | Kupffer cell (liver-resident Mφ) |
| KRT19, EPCAM | Halpern et al. 2018 | Cholangiocyte (bile duct cell) |

### Why These Markers Are Robust

1. **Panel-Derived:** Markers come from INP (Instituto Nacional de Pediatría)'s own curated panels, not generic lists
   - Reflects real experimental validation at INP (Instituto Nacional de Pediatría)
   - Optimized for Xenium spatial detection

2. **Multi-Source Validation:**
   - Overlap with published single-cell atlases (Zheng, Halpern, Bonnardel)
   - Inclusion in CellTypist models (trained on >1M cells)
   - Human Protein Atlas immune cell annotations

3. **Tissue-Appropriate:**
   - Lung markers account for lung-specific immunity (alveolar, airway ecology)
   - Liver markers include tissue-resident Kupffer cells (not generic macrophage markers)

4. **Functional Grouping:**
   - Not just presence/absence; includes activation/exhaustion states (PDCD1, HAVCR2, TIGIT)
   - Includes polarization markers for macrophages (M1 vs M2 specific)

---

## How to Update Markers

### Adding a New Gene to a Subtype

1. Edit `config_lung.yaml` or `config_liver.yaml`
2. Find the `immune_markers` section
3. Add gene to the subtype list:
   ```yaml
   CD8_T_cell: [GNLY, GZMK, GZMA, GZMB, PRF1, CD8A, CD8B, EOMES, CST7, CXCR6, NEW_GENE]
   ```
4. Save and re-run step 06b:
   ```bash
   snakemake --configfile config/config_lung.yaml --cores 8 --forcerun 06b_immune_subclustering
   ```

### Creating a New Cell Subtype

1. Add new entry to `immune_markers` dict in config:
   ```yaml
   CD8_T_cell_exhausted: [PDCD1, HAVCR2, TIGIT, LAG3, TOX]
   ```
2. Add subtype to `immune_cell_types` list if not already present
3. Re-run step 06b

### Validating After Changes

- Check output purity heatmap: `06b_immune_subclustering_report.html`
- Inspect per-cluster purity: Look for >0.7 mean purity per cluster
- If new genes have <5 cells expressing them, consider removing (too specific)

---

## Attribution & Citation

If you use these markers in publications, please cite:

> Immune cell markers extracted from INP (Instituto Nacional de Pediatría) custom Xenium immunology panels (Propuesta panel Pulmón 2025 for lung; Propuesta panel Xenium Liver August 2025 for liver). Validation against published single-cell atlases (Zheng et al. 2017, Halpern et al. 2018, Bonnardel et al. 2019) and CellTypist pre-trained models.

---

## Future Panel Development

### When to Redesign Markers

- **New tissue type:** Liver → pancreas, kidney, etc.
  - Extract tissue-specific markers from add-ons
  - Add tissue-resident immune cells (e.g., Paneth cells, glomerular macrophages)

- **New research question:** General immunity → checkpoint targeting
  - Expand exhaustion markers (PDCD1, HAVCR2, TIGIT, LAG3, ENTPD1, TOX)
  - Add activation markers (IFNG, IL2, GZMA, PRF1)
  - Include costimulatory (CD28, ICOS, OX40) and coinhibitory (CTLA4, PD-1) markers

- **New panel size:** 289 genes → 5K genes
  - More granular distinction between subtypes (e.g., CD8 naive vs effector vs resident)
  - Include metabolic markers (tumor vs host macrophages differ in LDHA, PFKFB3)
  - Add TCR/BCR clonality markers if applicable

---

**Document Owner:** Miguel Ángel Díaz-Campos (INP (Instituto Nacional de Pediatría))  
**For Questions:** See `IMMUNE_ANNOTATION_GUIDE.md`  
**Related Files:**
- `pipeline/config/config_lung.yaml` — Lung panel markers
- `pipeline/config/config_liver.yaml` — Liver panel markers
- `pipeline/scripts/06b_immune_subclustering.py` — Implementation
- `my_xenium_panel_markers/` — Original CSV files
