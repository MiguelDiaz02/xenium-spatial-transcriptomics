# Dataset Overview

## Overview

These datasets are provided as part of the *Post-Xenium In Situ Applications: Immunofluorescence, H&E, Visium v2, and Visium HD Technical Note (CG000709, Rev C)*.

- **Experiment 1:** Xenium In Situ Gene Expression (Xenium v1) data for adult human lung adenocarcinoma tissue (FFPE) using the Xenium Human Lung Gene Expression Panel with nuclear expansion.
- **Experiment 2:** Xenium Prime 5K In Situ Gene Expression with cell segmentation data for human lung adenocarcinoma tissue (FFPE) using the Xenium Prime 5K Human Pan Tissue and Pathways Panel.

The corresponding Visium HD Spatial Gene Expression datasets that were performed after these Xenium experiments are available separately.

---

## Metadata

| Category | Value |
|----------|------|
| Product | In Situ Gene Expression |
| Additional Applications | Cell Segmentation Staining |
| Software | Xenium Onboard Analysis |
| Instrument | Xenium Analyzer |
| Species | Human |
| Anatomical Entity | Lung |
| Preservation Method | FFPE |
| Disease State | Lung cancer |
| Cells or Nuclei | Cells |
| Biomaterial Type | Specimen from Organism |
| Donor Count | 1 |
| Date Published | 2024-11-06 |

---

## How to View Data

Interactively explore data with Xenium Explorer by downloading the Xenium Output Bundle (or subset). The bundle includes:

- `experiment.xenium`
- `gene_panel.json`
- `morphology_focus/`
- multi-file OME-TIFF files
- `analysis_summary.html`
- `cells.zarr.zip`
- `cell_feature_matrix.zarr.zip`
- `transcripts.zarr.zip`
- `analysis.zarr.zip`

Refer to the *Getting Started with Xenium Explorer* guide for details.

---

## Biomaterials

FFPE-preserved tissue was obtained from Avaden Biosciences  
(Lung Cancer, Invasive Acinar Adenocarcinoma, IB, T2a N0 MX).

---

## Tissue Preparation

Tissues were prepared following the *Xenium In Situ for FFPE - Tissue Preparation Handbook (CG000578)*.

- **Experiment 1:** Probe hybridization, washing, ligation, and amplification followed the Xenium In Situ Gene Expression User Guide (CG000582).
- **Experiment 2:** Same steps plus cell segmentation staining using the Xenium Prime protocol (CG000760).

Post-instrument processing followed the Xenium In Situ Gene Expression Post-Xenium Analyzer H&E Staining protocol (CG000613).

---

## Gene Panels

- **Xenium Human Lung Gene Expression Panel:** Designed to identify epithelial, stromal, endothelial, secretory, immune, and tumor cells across healthy, cancer, and fibrotic lung.
  
- **Xenium Prime 5K Human Pan Tissue and Pathways Panel:** Enables comprehensive cell type and state identification. Includes:
  - Canonical signaling pathways
  - Developmental biology genes
  - Immuno-oncology markers
  - Well-characterized biomedical genes

---

## Xenium Analyzer

- Protocol followed: *Xenium Analyzer User Guide (CG000584)*
- On-instrument analysis: Xenium Onboard Analysis v3.0.0

---

## Metrics

| Metric | Experiment 1 (Xenium v1) | Experiment 2 (Xenium Prime 5K) |
|--------|------------------------|-------------------------------|
| Median transcripts per cell | 58 | 242 |
| Cells detected | 278,659 | 278,328 |
| Nuclear transcripts per 100 µm² | 162.2 | 674.5 |
| High quality decoded transcripts | 30,871,547 | 149,372,743 |
| Region area (µm²) | 72,022,935.6 | 71,085,050.3 |

---

# Additional Information

## Key Metrics

| Metric | Value |
|--------|------|
| Number of cells detected | 278,659 |
| Median transcripts per cell | 58 |
| Nuclear transcripts per 100 µm² | 162.2 |
| Total high quality decoded transcripts | 30,871,547 |

---

## Sample Region Summary

| Field | Value |
|------|------|
| Region name | Human_Lung_Cancer_FFPE |
| Slide ID | N/A |
| Cassette name | N/A |
| Preparation method | ffpe |

---

## Run Information

| Field | Value |
|------|------|
| Run name | Human Lung Cancer Xenium v1 |
| Cell segmentation | Nuclei (DAPI) |
| Run start time | Apr 27, 2024, 00:05 GMT |
| Region area (µm²) | 72,022,935.6 |
| Total cell area (µm²) | 30,605,608.8 |

---

## Software

| Field | Value |
|------|------|
| Instrument software version | R&D |
| Analysis version | xenium-3.0.0.15 |
| Instrument serial number | R&D |

---

## Panel Specification

| Field | Value |
|------|------|
| Panel name | Xenium Human Lung Gene Expression |
| Design ID | hLung_v1.1 |
| Created by | 10x Genomics |
| Date created | 10/11/23 |
| Panel type | Predesigned |
| Tissue type | Human Lung |
| Chemistry version | Xenium v1 |
| Number of target genes (RNA) | 289 |
