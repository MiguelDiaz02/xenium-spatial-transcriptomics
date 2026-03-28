# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a collection of spatial transcriptomics (ST) analysis resources focused on 10x Genomics Xenium data. The primary codebase is `xenium_code_knowledge/Xenium_5k_analysis_pipeline/`, a reproducible Python pipeline (`recode_st` package). The other subdirectories contain reference codebases, notebooks, and documentation.

## Primary Project: Xenium_5k_analysis_pipeline

### Setup

```bash
cd xenium_code_knowledge/Xenium_5k_analysis_pipeline

# Conda (recommended)
conda env create -f environment.yml
conda activate xenium_5k_venv
pip install --no-build-isolation --no-deps -e .

# Or venv
python -m venv recode_st && source recode_st/bin/activate
pip install -r requirements.txt && pip install -e .

# Dev dependencies
make install-dev
```

### Common Commands

```bash
# Run tests
pytest tests/
pytest -v -p no:warnings --cov=src --cov-report=html --doctest-modules

# Run a single test file
pytest tests/test_config.py

# Docs server
make serve   # mkdocs at localhost:8000

# Lint (ruff)
ruff check src/
ruff format src/
```

### Architecture

The pipeline is in `src/recode_st/` and follows a sequential analysis workflow:

1. **Entry & Config** — `__main__.py` orchestrates the pipeline; `config.py` uses Pydantic for typed configuration
2. **Ingest** — `format_data.py` converts Xenium output → SpatialData Zarr store; `subsample_data.py` for downsampling
3. **QC** — `qc.py` computes and filters quality metrics
4. **Dimensionality reduction & clustering** — `dimension_reduction.py` (PCA/UMAP/Leiden)
5. **Annotation** — `annotate.py` for cell-type labeling
6. **Integration** — `integrate_scvi.py` (scVI-tools), `integrate_ingest.py` (scanpy ingest)
7. **Spatial analysis** — `spatial_statistics.py` (Squidpy), `muspan.py` / `ms_spatial_graph.py` / `ms_spatial_stat.py` (MuSpAn)
8. **Visualization** — `view_images.py` for spatial image overlays
9. **Downstream** — `drug2cell.py`, `doublet_identification.py`, `psuedobulk.py`, `denoise_resolvi.py`

Key data format: **SpatialData** (Zarr-backed). All intermediate and final results are stored as `SpatialData` objects.

### Dependencies to know

- `spatialdata` / `spatialdata-io` — core data format
- `scanpy` + `anndata` — single-cell analysis
- `squidpy` — spatial statistics
- `torch` — required by scVI and ResolVI
- `MuSpAn` — optional; controls spatial graph/statistics modules

## Other Codebases (reference/learning)

| Directory | Language | Purpose |
|-----------|----------|---------|
| `xenium_code_knowledge/sopa/` | Python | Technology-agnostic spatial omics pipeline with CLI (`sopa`) and Snakemake integration |
| `xenium_code_knowledge/spatialdata_xenium_explorer/` | Python | Convert SpatialData to Xenium Explorer format |
| `xenium_code_knowledge/Xenium_benchmarking/` | Python | Benchmarking 25+ datasets; `xb` module |
| `xenium_code_knowledge/XeniumIO/` | R | Bioconductor package for reading Xenium data |
| `xenium_code_knowledge/XeniumSpatialAnalysis/` | R | Spatial gradient and clustering analysis |
| `xenium_code_knowledge/liver_ped_map/` | Python/R | Pediatric liver atlas + IFALD spatial analysis |
| `xenium_code_knowledge/analysis_guides/` | Python/R | Official 10x Genomics tutorial notebooks |

### sopa CLI (for reference)
```bash
pip install sopa
sopa --help
# Snakemake pipeline: sopa/workflow/
```
