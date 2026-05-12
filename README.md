# Xenium spatial transcriptomics: multi-method consensus pipeline

A reproducible Snakemake + conda pipeline for targeted *in situ* spatial transcriptomics (10x Genomics Xenium), built around a **multi-tool consensus framework**: every analytical stage runs ≥2 independent tools and retains only findings concordant across methods. Pilot demonstration on a publicly available human lung adenocarcinoma dataset; engineered for portability to multi-sample cohorts.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Citation

If you use this pipeline, please cite the manuscript (currently under review at *Nature Communications*) and this repository — see [`CITATION.cff`](CITATION.cff).

## What the pipeline does

| Stage | Tools (consensus) | Output |
|---|---|---|
| Cell-type annotation | Immune subclustering + Leiden + marker scoring | 26 granular cell types, ARI=0.742 (internal consistency) |
| Spatially variable genes | Hotspot + nnSVG + Moran's I | 288/289 consensus SVGs |
| Spatial domains | Banksy + Novae | 14+10 partition; 35.4% per-cell consensus, ARI=0.15 |
| Cell-cell communication | LIANA+ (5 methods) + Spacia (Bayesian MIL) | 20/30 dual-validated interactions |
| Pseudotime | Slingshot + Diffusion Pseudotime | Macrophage activation continuum (r=0.478) |

**Key biological finding (pilot):** ADAM17 → MUC1 hub across 6 myeloid/stromal senders → tumor cells, Spacia p_adj up to 10⁻⁷¹. Both M1- and M2-like macrophages participate (non-directional polarization continuum).

## Pipeline status

- **Pilot (1 sample, 289 genes):** complete and audited. See `manuscript/` for the manuscript PDF and supplementary information.
- **Cohort port (P0 + P1):** infrastructure complete and tested (71/71 pytest + end-to-end smoke test). Ready to receive the TBDs cohort TMAs (2 slides, 28 donors, 439-gene panel each). See `PIPELINE_PORTABILITY_CHECKLIST.md` for the full audit trail.

## Quick start

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/MiguelDiaz02/xenium-spatial-transcriptomics.git
cd xenium-spatial-transcriptomics

# 2. Create the main conda env
conda env create -f pipeline/envs/xenium_pipeline.yaml
conda activate xenium_pipeline

# 3. (Optional) the R env and Spacia env
conda env create -f pipeline/envs/xenium_R_analysis.yaml
conda env create -f pipeline/envs/spacia.yaml

# 4. Dry-run the workflow
cd pipeline
snakemake --configfile config/config_lung.yaml -n

# 5. Run end-to-end smoke test (synthetic mini-cohort, ~16 s)
conda run -n xenium_pipeline python tests/smoke_test_p1.py
```

## Repository layout

```
proyecto_demo_xenium/
├── manuscript/                  # LaTeX manuscript (EN + ES) + cover letters
├── pipeline/
│   ├── Snakefile                # main workflow
│   ├── config/                  # cohort + dataset configs
│   │   ├── config_lung.yaml     # pilot lung-cancer pipeline params
│   │   ├── cohort_TBDs.yaml     # multi-sample TMA cohort manifest
│   │   └── markers/             # per-tissue cell-type marker YAMLs
│   ├── envs/                    # conda environment yamls (per stage)
│   └── scripts/
│       ├── 01_ingest.py ... 12_export.py   # core pipeline stages
│       ├── analysis/            # F0..F8 audited analysis modules
│       ├── integration/         # P1 multi-sample (concat, scVI, pseudobulk × 2, LIANA+ multi)
│       └── utils/               # paths.py, cohort.py
├── tests/
│   ├── test_p0_portability.py   # 44 tests
│   ├── test_p1_multisample.py   # 27 tests
│   ├── smoke_test_p1.py         # end-to-end synthetic-cohort smoke test
│   └── fixtures/synthetic_cohort.py
├── human_lung_cancer/           # pilot data + results (NOT in this repo; available via 10x CG000709 Rev C)
├── PIPELINE_PORTABILITY_CHECKLIST.md
├── LICENSE                      # MIT
├── CITATION.cff
└── README.md
```

## Data

- **Pilot:** Public 10x Genomics dataset CG000709 Rev C (human lung adenocarcinoma, Xenium v1 chemistry, hLung_v1.1 panel, 289 genes, 268,034 cells after QC). Download from [10x Genomics](https://www.10xgenomics.com/datasets) and place under `human_lung_cancer/raw/`.
- **TBDs cohort (forthcoming):** 2 Xenium TMA slides, 28 donors (16 lung + 12 liver), 439-gene panels per organ. See `pipeline/config/cohort_TBDs.yaml` for the donor manifest.

## Reproducing the pilot

```bash
cd pipeline
snakemake --configfile config/config_lung.yaml --cores 8
```

End-to-end runtime: ~12 hours on a workstation (8 cores, 64 GB RAM). The Spacia stage is the bottleneck (~6 h) and runs in a containerized environment — see `pipeline/envs/spacia.yaml` and the Spacia Docker image (`spacia:audited`).

## Methodological transparency

This pipeline went through a formal audit (documented in `manuscript/manuscript.tex` Supplementary Note S1). Five issues were identified and resolved before publication:

1. tradeSeq applied to embeddings cast as integer counts → restricted to count-based macrophage subset.
2. Slingshot start/end clusters set to biologically implausible configuration → replaced with macrophage M1-rooted single lineage.
3. LIANA+ "FDR<5%" threshold was actually a magnitude-rank cutoff → corrected to BH-FDR over CellPhoneDB p-values.
4. Spacia per-test p-values reported as validation without cross-test correction → Bonferroni applied across the 30 submitted interactions.
5. Software-version drift between code and manuscript → all versions pinned in `pipeline/envs/*.yaml` and reported verbatim in Methods.

## Cohort port (P0 + P1)

The codebase is engineered to scale from the single-sample pilot to multi-sample TMA cohorts without rewriting analysis logic. See `PIPELINE_PORTABILITY_CHECKLIST.md` for the gates (B1–B5 + M1–M5) and `tests/smoke_test_p1.py` for the end-to-end validation.

## Authors

- **Miguel Ángel Díaz-Campos** — Instituto Nacional de Pediatría (INP), México. ([mdiazc161@unam.edu](mailto:mdiazc161@unam.edu))
- **Alfredo de Jesús Rodríguez Gómez** (corresponding) — INP / Instituto de Investigaciones Biomédicas, UNAM. ([alfredo.rodriguez@iibiomedicas.unam.mx](mailto:alfredo.rodriguez@iibiomedicas.unam.mx))

## License

[MIT](LICENSE).
