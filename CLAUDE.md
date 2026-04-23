# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Modular Snakemake pipeline for end-to-end 10x Genomics Xenium spatial transcriptomics analysis.
Designed as a pilot with human lung cancer (Xenium v1, 289 genes) and built to scale to liver
(Xenium Prime 5K) via config swap only — no code changes required.

## Setup

```bash
cd pipeline
conda env create -f envs/xenium_pipeline.yaml
conda activate xenium_pipeline

# Verify GPU (required for Cellpose + ResolVI)
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

Baysor (optional, transcript-based segmentation) requires Julia:
```bash
curl -sSL https://install.julialang.org | sh
julia -e 'using Pkg; Pkg.add("Baysor")'
pip install baysorpy
```

## Running the Pipeline

```bash
cd pipeline

# Dry run — check DAG without executing
snakemake --configfile config/config_lung.yaml -n

# Full run (8 cores, uses GPU via resources)
snakemake --configfile config/config_lung.yaml --cores 8

# Resume from a specific step (e.g., after fixing QC thresholds)
snakemake --configfile config/config_lung.yaml --cores 8 --forcerun qc

# Visualize the DAG
snakemake --configfile config/config_lung.yaml --dag | dot -Tpdf > dag.pdf

# Run for liver (future) — only the config changes
snakemake --configfile config/config_liver.yaml --cores 8
```

## Architecture

```
pipeline/
├── Snakefile              # Orchestrator; auto-selects segmentation at DAG-build time
├── config/
│   ├── config_lung.yaml   # Active config (289-gene lung demo)
│   └── config_liver.yaml  # Future (5K panel, different CellTypist model)
├── rules/                 # One .smk per pipeline step
├── scripts/               # One .py per step; use snakemake.input/output/params
│   └── utils/             # Shared: io.py (read/write zarr), logging.py (rich)
└── envs/
    └── xenium_pipeline.yaml
```

### Pipeline Steps

| Step | Rule | Key output |
|------|------|------------|
| 01 | `ingest` | `sdata.zarr` (SpatialData with H&E registered) |
| 02 | `segmentation` | `cell_boundaries_cellpose` or `cell_boundaries_baysor` in zarr |
| 02b | `compare_segmentation` | `02_segmentation_comparison.html` |
| 03 | `qc` | Filtered table + `03_qc_report.html` |
| 04 | `preprocess` | `layers["counts"]` (raw) + normalized `X` |
| 05 | `reduction` | `obsm["X_pca"]`, `obsm["X_umap"]`, `obs["leiden"]` |
| 06 | `annotation` | `obs["cell_type"]`, `obs["cell_type_fine"]` |
| 07 | `denoising` | `layers["denoised"]` (raw counts preserved) |
| 08 | `spatial` | Moran's I, neighborhood enrichment, co-occurrence |
| 09 | `downstream` | Doublet scores, pseudobulk DE CSVs |

### Auto-segmentation logic (Snakefile, runs at DAG-build time)

Reads `experiment.xenium` directly from inside the ZIP before any rule executes:
- `median_tx/cell >= 100` **AND** `panel_genes >= 500` → **Baysor**
- Otherwise → **Cellpose**

Current dataset (289 genes, 58 tx/cell): **Cellpose**
Future liver 5K dataset (~242 tx/cell, ~5000 genes): **Baysor**

Override with `segmentation.method: "cellpose"` or `"baysor"` in config.

### Single SpatialData zarr store pattern

All steps read from and write back to `results/sdata.zarr`.
Each script only writes the elements it modifies (using `sdata.write_element()`).
Done files (`results/NN_stepname.done`) track completion for Snakemake.

### Step 09 (Downstream) Performance Notes

**Doublet detection (Scrublet)** is computationally intensive, especially for large datasets.

- **Small datasets (<50k cells)**: ~5–10 minutes
- **Medium datasets (50k–100k)**: ~15–30 minutes
- **Large datasets (>100k cells, e.g., 268k)**: 45–90+ minutes

**Optimizations applied** (as of 2026-04-23):
1. Reduced PCA components from 30 → 15 for datasets > 100k cells (trades marginal accuracy for speed)
2. Timeout increased to 3600s (1 hour) and memory to 32 GB
3. Progress logging added to monitor execution

**Tuning for future datasets:**
- If a dataset is >200k cells and Scrublet still times out, reduce PCA further (e.g., `n_pca = 10`)
- For rapid prototyping, set `doublet_detection: false` in config and add doublets later as an optional refinement step

### Layer conventions (critical)

| Layer | Content | Used by |
|-------|---------|---------|
| `adata.X` | Normalized + log1p | PCA, UMAP, CellTypist, Squidpy |
| `adata.layers["counts"]` | Raw integer counts | scVI, ResolVI, DESeq2 |
| `adata.layers["denoised"]` | ResolVI output | Visualization, spatial autocorr |

**Never use `layers["denoised"]` for DESeq2 or scVI — always use `layers["counts"]`.**

## Adapting to Liver (5K panel)

1. Copy `config/config_lung.yaml` → `config/config_liver.yaml`
2. Update paths to liver data files
3. Change `annotation.celltypist_model` to `Liver_Human_PIP`
4. Adjust `qc` thresholds (5K panels have higher counts; raise `max_transcripts`)
5. Set `segmentation.method: auto` (will auto-select Baysor for 5K)
6. Set `downstream.pseudobulk: true` if using multiple liver samples

No rule or script changes required.
