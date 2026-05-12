# Pipeline Portability Checklist — TBDs cohort (lung + liver fibrosis)

**Date:** 2026-05-11 (cohort numbers updated)
**Source dataset (this manuscript):** 1 sample, 268,034 cells, 289 genes, lung adenocarcinoma
**Target dataset (incoming):** 10 slides — **4 lung + 6 liver** — TBDs (dyskeratosis congenita) and non-TBDs controls.
**Panels:** 2 distinct Xenium panels (lung-specific and liver-specific), 439 genes each (389 base + 50 add-on per organ).

**Cohort design (locked 2026-05-11):**
- Lung: 4 slides → {healthy_TBD, healthy_nonTBD, fibrotic_TBD, fibrotic_nonTBD}
- Liver: 6 slides → 3 healthy + 3 cirrhotic, split across TBD / non-TBD status
- Condition contrasts: `healthy` vs `fibrotic_nonTBD` (disease baseline) and `fibrotic_nonTBD` vs `fibrotic_TBD` (TBD-specific signature on top of fibrosis)

This checklist enumerates **every hardcoded, sample-specific, or panel-specific assumption** in the pilot pipeline and the change required to run cleanly on the new cohort. Action items are sorted by priority.

---

## P0 — Blockers (must fix before any new run)

### B1. Hardcoded absolute paths in every script ✅ DONE (2026-05-08)
**Files:** `pipeline/scripts/analysis/F0_reannotation_v3.py`, `F0b_*`, `F0c_*`, `F1_*`, `F2_*`, `F3_*`, `F5_*`, `F8_*`, `Fmanuscript_figures.py`, `F_master_figure.py`
**Resolution:** Created `pipeline/scripts/utils/paths.py` (Python) + `paths.R` (R). Resolution order:
1. `XENIUM_PROJECT_ROOT` env var (explicit override)
2. Walk up from script location
Dataset selection via `XENIUM_DATASET` (default: `human_lung_cancer`). Spacia binary via `SPACIA_PATH` (default: `pipeline/external/Spacia/spacia.py`). All 15 Python scripts and 3 R scripts migrated. Shell scripts F3b/F3d migrated. **No script contains `/home/mdiaz` literal anymore.**

### B2. Cell-type marker list embedded in `F0_reannotation_v3.py` ✅ DONE (2026-05-08)
**Resolution:** Moved markers, L1/L3 hierarchy mappings, immune-granular candidates, non-immune scoring set, Leiden-direct table and annotation parameters to YAML config:
- `pipeline/config/markers/lung_pilot.yaml` — canonical for the lung adenocarcinoma pilot (matches v3 results in `sdata.zarr`).
- `pipeline/config/markers/lung_TBDs.yaml` — TEMPLATE for TBDs lung fibrosis cohort (no tumor markers; AT2_dysfunctional/Myofibroblast added; verify against panel manifest before use).
- `pipeline/config/markers/liver_TBDs.yaml` — TEMPLATE for TBDs liver fibrosis (zonal hepatocytes, hepatic stellate quiescent→activated, Kupffer, LSEC, portal fibroblasts).

`F0_reannotation_v3.py` now accepts `--markers <yaml>`; default is `lung_pilot.yaml`. Backward compatible: same numbers as the canonical run when `--markers` defaults are used.

### B3. Panel-size-dependent thresholds ✅ DONE (2026-05-08)
**Resolution:** Added `qc.max_cell_area_um2` to `config_lung.yaml`. Spacia params (mcmc, n_cells_subsample, dist_cutoff, top_n_pairs) lifted from F3 hardcoded constants to `analysis.f3_ccc.spacia_validation`. F3 now reads from config with CLI override (`--config`, `--mcmc_params`, `--n_cells`, `--dist_cutoff`, `--top_n`). Removed `spacia_repo: /home/mdiaz/tools/Spacia` from config (now resolved via `SPACIA_PATH` env var or `pipeline/external/Spacia/`). Corrected MCMC iter count from `5000` (stale config) to `50000,20000,100,2` (audit-correct, matches manuscript Note S1.5). LIANA params (n_perms, min_cells, expr_prop) already CLI-parameterized in F0c.

**Reference table (panel-size-dependent thresholds):**

| Where | Current value | Justification on 289 panel | Action for 439 panel |
|---|---|---|---|
| QC: `min_counts ≥ 5` per cell | scanpy default | OK on 289 | Re-evaluate; consider scaling to ~7–8 to maintain similar percentile filtering on 439 |
| QC: `cell_area ≤ 1500 µm²` | inherited from 10x defaults | Sample-specific tissue morphology | TBDs lung tissue may have different morphology (fibrotic remodeling); reset by histogram |
| LIANA+ `expr_prop = 0.1` | floor for L-R expression | OK for 289 (loose floor) | OK for 439 too |
| LIANA+ `min_cells = 50` per cluster | minimum cell-type abundance | OK | OK |
| Spacia `N_CELLS = 5000` per type | tractability cap | Conservative | Can stay; may relax for rarer types |
| Spacia `DIST_CUTOFF = 30 µm` | juxtacrine range | Cell-cell biophysics | Same |
| Spacia `MCMC = 50000,20000,100,2` | publishable | Same | Same |

**Effort:** 1 h to expose these as config; 2 h to revalidate on first TBDs sample.

### B4. NicheDE, hdWGCNA, SCENIC, CellTypist all gated by gene count ✅ DONE (2026-05-08)
**Resolution:** Documented in manuscript **Supplementary Note S1.1** ("Methods discarded due to gene-panel constraints"). Section explicitly enumerates the four tools, why each fails at 289 genes, and that the 439-gene TBDs panel is still below the typical lower bound (~1,000 genes). Reviewers and collaborators have a transparent record.

### B5. Spacia tool path hardcoded ✅ DONE (2026-05-08)
**Resolution:** `paths.py::spacia_tool_path()` resolves Spacia binary in this order:
1. `SPACIA_PATH` env var (explicit override).
2. `pipeline/external/Spacia/spacia.py` (if Spacia is cloned in-repo as a submodule).
F3_spacia_ccc_validation.py + F3b/F3d shell scripts all use this resolution. **Docker image already built (`spacia:audited`, 4.4 GB)** — see manuscript Note S1.5; recommended path for reproducibility on new machines. Singularity recipe (`tools/Spacia/singularity_build.def`) exists but image not yet built.

---

## P1 — Multi-sample integration (NEW for TBDs cohort)

### M1. Multi-sample integration ✅ DONE (2026-05-11)
**Resolution:** Built cohort infrastructure for the TBDs TMA design (2 TMA slides, 28 donors, ~115 cores, 2 organs, 3 condition levels). New components:
- `pipeline/config/cohort_TBDs.yaml` — donor manifest + condition contrasts (lung: TBD vs control; liver: TBD vs control, TBD vs Alc_Cirh, Alc_Cirh vs control).
- `pipeline/scripts/utils/cohort.py` — typed loader (`load_cohort`, `expand_cores`) that walks donors → cores at runtime via `obs.core_id` on the per-TMA sdata.
- `pipeline/scripts/utils/paths.py` extensions (`cohort_root`, `sample_sdata_path`, `cohort_yaml_path`).
- `integration/concat_samples.py` — reads per-TMA sdata, expands cores, concats with `obs.{sample_id,subject_id,organ,tma_slide,condition,tbd_status,control_subtype}` and namespaced `obsm[spatial_<sid>]`.
- `integration/scvi_integrate.py` — scVI with `batch_key=sample_id` + `subject_id` as categorical covariate (correction for core-level batch without collapsing donor signal).
- `integration/pseudobulk_decoupler.py` — `decoupler.get_pseudobulk` + PyDESeq2 per (cell_type, donor); contrasts pulled from cohort YAML.
- `integration/pseudobulk_muscat.R` — `muscat::aggregateData` + edgeR/DESeq2/limma trio (R counterpart).
- `integration/pseudobulk_consensus.py` — joins both methods, Fisher-combined p, Jaccard concordance, Venn figures, consensus classification (strict_consensus / majority / single / discordant).

Statistical unit = donor (subject_id); spatial unit = core (sample_id). Pseudo-replication addressed by aggregating cores per donor before DE.

### M2. Multi-sample CCC (LIANA+) ✅ DONE (2026-05-11)
**Resolution:** `integration/F0c_ccc_liana_multisample.py` runs `liana.mt.rank_aggregate` per donor, stacks results, then Mann-Whitney U + BH-FDR for condition contrasts on `magnitude_rank` per (sender, receiver, L–R) tuple. Outputs per-donor TSVs + per-organ aggregated + per-contrast significance tables.

### M3. Spacia per-core + meta-analysis ✅ SCAFFOLD (2026-05-11)
**Resolution:** `integration/spacia_meta.py` — full interface, dry-run path executable (writes `plan.json` listing ~700 expected jobs); Stouffer/Fisher meta-analysis helpers (`meta_across_cores`, `meta_across_donors`) implemented. `--execute` path documented inline; needs spacia conda env + F3 wiring (~1 day).

### M4. Cohort pseudotime ✅ SCAFFOLD (2026-05-11)
**Resolution:** `integration/F5_pseudotime_cohort.R` — dry-run executable; full Slingshot + `tradeSeq::conditionTest` algorithm laid out inline as commented blocks ready for activation once `xenium_R_analysis` env is verified on the cohort h5ad.

### M5. Spatial domain harmonization ✅ SCAFFOLD (2026-05-11)
**Resolution:** `integration/F1_novae_crosssample.py` — two execution modes (Novae multi-slide and Banksy + post-hoc matching). Dry-run prints the cohort plan; execute paths sketched with required env (GPU + novae ≥0.4).

### Smoke test (end-to-end P1 validation) ✅ DONE (2026-05-11)
**Resolution:** Built synthetic mini-cohort fixture (`tests/fixtures/synthetic_cohort.py`, 2 mini-TMAs / 10 donors / 20 cores / 5000 cells / 48 genes selected for LIANA L–R coverage) plus runner (`tests/smoke_test_p1.py`) that exercises concat → pseudobulk_decoupler → pseudobulk_consensus → LIANA+ multisample + M3/M5 dry-runs as real subprocess CLIs.

3 production API bugs found and fixed during the smoke test:
- **decoupler 2.x**: `dc.pp.pseudobulk` replaces `dc.get_pseudobulk`; column renamed `psbulk_n_cells` → `psbulk_cells`. Fixed in `pseudobulk_decoupler.py` with defensive column lookup.
- **LIANA `use_raw`**: explicit `use_raw=False` so LIANA reads from `adata.X` (log1p) instead of uninitialized `.raw`.
- **LIANA panel-resource mismatch**: added `_filter_resource_to_panel()` to pre-filter the L–R resource to genes available in the panel. Required for sparse Xenium panels (289/439 genes) where >98% of the consensus resource is absent.

Numerical validation: PyDESeq2 recovered injected condition effects (TGFB1 log2FC=2.05, padj=8.2e-70 in lung AT1 TBD-vs-control matches the simulated 3× upregulation). LIANA detected ADAM17→MUC1 (lung-cancer pilot hub, expected by panel design).

Run cadence: re-run the smoke test whenever `scvi-tools`, `decoupler`, `liana`, or `pydeseq2` are updated, to catch API breaks before they hit real data.

### M0. The pilot pipeline had zero multi-sample logic (legacy issue ABOVE has been resolved by M1–M5)
**Issue:** Every script reads a single `sdata.zarr` corresponding to one tissue section. With 10 samples (4 lung + 6 liver, 2 distinct Xenium panels), the pipeline must support per-sample QC + integration with panel-aware feature handling.  
**Required new components:**
- `pipeline/scripts/integration/concat_samples.py` — concatenate per-sample SpatialData into one multi-sample object with `obs.sample_id`, `obs.condition` (`fibrotic` vs `healthy`), `obs.organ` (`lung` vs `liver`), `obs.subject_id`.
- `pipeline/scripts/integration/scvi_integrate.py` — scVI batch correction (already in env: `scvi-tools 1.4.2`). Use `sample_id` as `batch_key`. Train with `n_layers=2`, `n_latent=30`, `gene_likelihood='nb'`.
- Sample-aware DGE: `F0b` rewritten to use `sc.tl.rank_genes_groups` with `pts=True` and explicit per-sample replicates, OR (preferred) pseudobulk with `decoupler`/`muscat` aggregation per (cell_type, sample) before DE testing — addresses pseudo-replication directly.

**Effort:** 8 h — substantial new code.

### M2. CCC (LIANA+) does not currently aggregate across samples
**Issue:** LIANA+ has built-in support for multi-sample analysis via `liana.method.rank_aggregate(adata, sample_key='sample_id', ...)`. The pilot script does not use this.  
**Fix:** Update `F0c` to accept `--sample_key` parameter; when set, run sample-level CCC and aggregate via `liana.mt.adata_to_views()`.  
**Effort:** 3 h.

### M3. Spacia is single-sample by design
**Issue:** Spacia operates on one `(counts, metadata)` pair per run.  
**Fix:** Run Spacia per sample, then aggregate with mixed-effects meta-analysis across the 6 fibrotic samples; report consistency (number of samples in which each interaction is validated).  
**Effort:** 4 h.

### M4. Pseudotime requires careful multi-sample handling
**Issue:** Pooling all samples blurs trajectories; per-sample trajectories may not align.  
**Fix:** Approach for TBDs cohort:
1. Annotate per-sample, integrate with scVI.
2. For pseudotime: subset to the cell type of interest (e.g., fibroblasts → myofibroblasts in lung; hepatic stellate → activated stellate in liver) across **all 12 samples**.
3. Compute trajectory on integrated latent (DPT and Slingshot, as in F5b/F5c).
4. Test for fibrotic-vs-healthy differences along the trajectory using condition as covariate in `tradeSeq::conditionTest`.

**Effort:** 6 h.

### M5. Batch effects in spatial domain inference
**Issue:** Banksy and Novae are computed per-section. Cross-sample harmonization of "domain identity" is not trivial because tissue architecture differs across samples.  
**Fix:** For Banksy: run per-sample and post-hoc match domains by marker similarity. For Novae: use the model in cross-sample mode (Novae supports it natively when given multiple samples).  
**Effort:** 3 h + Novae cross-sample mode requires GPU re-training.

---

## P2 — Disease-specific adjustments (TBDs vs lung adenocarcinoma)

### D1. Trajectory expectations change radically
- **Pilot (lung cancer):** trajectory framed as TME gradient (tumor → parenchyma).
- **TBDs lung fibrosis:** the relevant trajectories are **AT2 → AT1 differentiation arrest** (telomere-mediated stem-cell exhaustion) and **fibroblast → myofibroblast activation** (`COL1A1`, `ACTA2`, `FAP` upregulation).
- **TBDs liver fibrosis:** **quiescent → activated hepatic stellate cell** trajectory (`LRAT` loss, `ACTA2` gain).

**Action:** Pre-register the trajectories of interest in `docs/TBDs_analysis_plan.md` to prevent post-hoc lineage selection.

### D2. CCC hubs of biological interest will differ
- **Pilot:** ADAM17--MUC1 (myeloid--tumor immunosuppression).
- **TBDs lung:** TGF-$\beta$ family (fibrotic activation), CXCL12--CXCR4 (immune recruitment), WNT (epithelial dysfunction).
- **TBDs liver:** PDGF-BB (stellate activation), TGF-$\beta$, lactate-mediated.

**Action:** Confirm the 439-gene panel includes these axes (TGFβ1, TGFBR1/2, COL1A1, ACTA2, PDGFB, PDGFRB, CXCL12, CXCR4 should all be present in an "immune + fibrosis" panel; verify against the panel manifest when received).

### D3. No tumor cells in TBDs samples
**Issue:** Several scripts (`F5_pseudotime_slingshot.R`, master figure caption, manuscript narrative) reference tumor cells. None of these will be present in TBDs samples.  
**Fix:** All TBDs analyses must use TBDs-appropriate trajectory anchors and CCC anchor cells. Treat the lung-cancer pilot as a stand-alone manuscript; do NOT carry tumor-centric framing into TBDs analysis.

---

## P3 — Code quality / repo hygiene (before public release)

### R1. Stale output directories in `02_biology/`
- `ccc_liana/`, `ccc_liana_L2_v2/`, `dgea_L2_v2/`, `reannotation/`, `reannotation_v2/`, `sdata_DESACTUALIZADO_NO_USAR.zarr` — should be removed before public release of the GitHub repo.
- Keep only the v3 / canonical paths.
- **Effort:** 30 min.

### R2. Versioning convention
**Issue:** Some scripts have `_v2`, `_v3` suffixes; some don't. Mixed naming.  
**Fix:** Adopt one of: (a) git tags for versions; (b) `vN` suffix consistently; (c) move all "v2" to `archive/` directory. Recommend (a) + (c).

### R3. Test coverage = 0 outside `recode_st`
**Issue:** None of the `pipeline/scripts/analysis/*` files have unit tests. This is excusable for one-off pilot scripts but risky for a production multi-sample pipeline.  
**Fix:** Add at minimum integration tests on a tiny synthetic dataset for each analysis stage. Use `pytest` + a small fixture under `tests/fixtures/`. Critical for confidence that the pipeline runs end-to-end on new data without silent failures.  
**Effort:** 2 days for full coverage; 4 h for a smoke test of every script.

### R4. Logging is inconsistent
**Issue:** Some scripts use `utils.logging.get_logger`, others use `logging.basicConfig` directly, others use `print` and `cat()`.  
**Fix:** Standardize on `utils.logging` across both Python and a parallel R helper.

### R5. Color schemes not validated for colorblind safety
**Issue:** Mixed use of matplotlib defaults, `Set3` (qualitative palette known to be poor on red-green colorblind), and Crameri (colorblind-safe).  
**Fix:** Audit every plot; standardize on Crameri palettes (`crameri.batlow` for sequential, `crameri.roma` for diverging, `crameri.vik` for divergent symmetric, `crameri.batlowS` for categorical) and remove `Set3`.  
**Effort:** 4 h.

### R6. Snakemake DAG completeness
**Issue:** Many of the new audit-fix scripts (`Faudit_recompute_stats.py`, `F5b_*`, `F5c_*`) are not yet in the Snakemake workflow.  
**Fix:** Add rules to `pipeline/Snakefile`. The DAG should include the entire audited pipeline, not just the original Phase A–E rules.  
**Effort:** 2 h.

---

## P4 — Reproducibility for public release

### G1. GitHub repository organization
Current state of `https://github.com/MiguelDiaz02/xenium-spatial-transcriptomics`: mixed bibliography + early code attempts + final pipeline. Per user spec, the formal release should contain:
- `pipeline/` — Snakemake workflow + scripts
- `envs/` — conda YAMLs (one per stage env)
- `config/` — sample configs (lung pilot + TBDs template)
- `docs/` — methods.md, marker-lists, this checklist
- `tests/` — minimal test harness
- `README.md` — install + quickstart
- `LICENSE`, `CITATION.cff`

Recommendation: create `release/v0.1` branch, copy only the curated subset, tag once stable.

### G2. Container/Singularity image for Spacia
**Issue:** Current setup requires manual Spacia clone + custom env.  
**Fix:** Use the provided `singularity_build.def` to build a self-contained image; document as the recommended path for the public release.

### G3. Data deposit
**Pilot data:** already public (10x CG000709 Rev C). No action needed.  
**TBDs data:** unclear public status. Plan now: if data is to remain private/in-house, the public repo can still reference the pipeline + provide a synthetic toy dataset; if public, deposit raw `.tar.gz` + processed `sdata.zarr` to GEO or Zenodo.

### G4. CITATION.cff + Zenodo DOI
For Nat Methods or Cell Reports Methods, a Zenodo DOI for the code at submission time is expected. Set up via GitHub-Zenodo integration once the repo is cleaned.

---

## P5 — Documentation deltas

### Doc1. `CLAUDE.md` (root and project)
**Update needed:**
- Reference this checklist.
- Note that `F5_pseudotime_slingshot.R` is DEPRECATED (replaced by F5b + F5c).
- Reference the new audit corrected stats path: `human_lung_cancer/results/02_biology/audit_corrected/`.

### Doc2. `SENDA_DORADA.md`
- Add Phase G: "Audit + corrections (2026-05-06)" with the recompute outputs and the M1↔M2 sub-trajectory.
- Mark Phase E as "superseded by F5b/F5c"; keep Phase E entry for history.

### Doc3. New: `docs/methods_excluded.md`
Document the four excluded methods (NicheDE, hdWGCNA, SCENIC, CellTypist) with the empirical reason for each (gene-count threshold, model coverage, etc.). This is what readers and reviewers will look for if they ask "why didn't you use X?".

---

## Summary — readiness gate for TBDs cohort

| Gate | Effort | Status |
|---|---|---|
| B1 Path lifting | 2 h | ✅ done (2026-05-08) |
| B2 Marker externalization | 6 h | ✅ done (2026-05-08) |
| B3 Threshold parameterization | 3 h | ✅ done (2026-05-08) |
| B4 NicheDE/hdWGCNA/SCENIC documentation | 30 min | ✅ done (manuscript Note S1.1) |
| B5 Spacia tool path / Singularity | 2 h | ✅ done — env-var + Docker image (Singularity recipe pending) |
| M1 Multi-sample concat + scVI | 8 h | ❌ blocking |
| M2 LIANA+ multi-sample | 3 h | 🟡 needed |
| M3 Spacia per-sample meta-analysis | 4 h | 🟡 needed |
| M4 Pseudotime cohort handling | 6 h | 🟡 needed |
| Other (R1–R6, G1–G4) | 1–2 days | 🟢 cleanup |

**P0 (B1–B5) status: COMPLETE.** Pipeline now portable to a new dataset by:
1. setting `XENIUM_PROJECT_ROOT` and `XENIUM_DATASET` env vars (or rely on auto-detection),
2. providing `pipeline/config/markers/<dataset>.yaml` (templates: `lung_TBDs.yaml`, `liver_TBDs.yaml`),
3. editing `pipeline/config/config_<dataset>.yaml` for QC + analytical thresholds,
4. (optional) setting `SPACIA_PATH` to a Spacia checkout, or run `spacia:audited` Docker image.

Remaining for TBDs cohort: P1 (multi-sample integration M1–M5) — ~21 h of focused dev work; recommend two weeks part-time before sample 1 arrives.
