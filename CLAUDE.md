# CLAUDE.md — proyecto_demo_xenium

**Leer SENDA_DORADA.md al inicio de cada sesión.**
**Dataset:** 268,034 células × 289 genes | Human lung cancer, Xenium v1
**Última actualización:** 2026-05-04

---

## Estado del pipeline

| Fase | Herramientas | Estado |
|------|-------------|--------|
| A — Anotación v3 híbrida | immune_granular + Leiden + score_genes | ✅ COMPLETO |
| B — Dominios espaciales | Novae (F8) + Banksy (F1) consensus | ✅ COMPLETO |
| C — SVG Consensus | Hotspot + nnSVG + Moran's I | ✅ COMPLETO |
| D — Niche-DE | NicheDE | ❌ DESCARTADO (289 genes insuf.) |
| E — Spatial Pseudotime | Slingshot + tradeSeq + GAM | ✅ COMPLETO 2026-04-30 |
| F — CCC Validation | LIANA+ + Spacia Bayesian MIL | ✅ COMPLETO 2026-04-30 |
| G — Annotation x-val | CellTypist 5 modelos + ARI interno | ✅ COMPLETO |
| Master Figure | 5 paneles integrados | ✅ COMPLETO 2026-04-30 |
| F7 — Co-expresión | hdWGCNA / SCENIC | ❌ DESCARTADO (289 genes insuf.) |
| **Manuscrito** | Nature Communications Research Article | 🔄 EN PROGRESO 2026-05-04 |

**Pendiente:** Escribir manuscrito (spec aprobado). Reporte técnico (interno) — secundario.

---

## Datos canónicos

| Recurso | Ruta |
|---------|------|
| **sdata.zarr (v3, ACTIVO)** | `human_lung_cancer/results/sdata.zarr` |
| DGEA L2 v3 | `human_lung_cancer/results/02_biology/dgea_L2_v3/` |
| CCC LIANA+ L2 v3 | `human_lung_cancer/results/02_biology/ccc_liana_L2_v3/` |
| Spacia results | `human_lung_cancer/results/02_biology/phase_f_spacia/` |
| SVG Consensus | `human_lung_cancer/results/03_phase3_spatial/F2_svg_consensus/` |
| Spatial domains | `human_lung_cancer/results/03_phase3_spatial/F1_spatial_domains/` + `F8_novae/` |
| **Master Figure** | `human_lung_cancer/results/figures/master_figure/master_figure_v1.png` + `.pdf` |

**NO usar:** `sdata_DESACTUALIZADO_NO_USAR.zarr`, `dgea_L2_v2/`, `ccc_liana_L2_v2/`

---

## Anotación v3 (canónica)

14 L1 · 26 L2 (25 bio + Unassigned 4.8%) · 7 L3 | ARI = **0.742**
Columnas: `cell_type_L1`, `cell_type_L2`, `cell_type_L3`, `annotation_source_v3`
CellTypist descartado como validador — 289 genes << ~20k del entrenamiento.

---

## Hallazgos Spacia (Fase F, 2026-04-30)

**22/25 interacciones LIANA+ validadas espacialmente (88%)**

| Eje | Top sender→receiver | pval_adj |
|-----|---------------------|---------|
| **ADAM17→MUC1** | AT1→Tumor_resting | 1.10e-71 |
| **LTF→AGER** | Tumor_resting→AT1 | 9.37e-66 |
| **CDH1→EGFR** | Epithelial→Endothelial_lymphatic | 6.64e-42 |
| **S100B→AGER** | cDC1→AT1 | 1.89e-24 |

5 FAILs irrecuperables: tipos raros (<500 células sender: Plasma, pDC×1, Plasmablast, DC_mature).

---

## Scripts canónicos (`pipeline/scripts/analysis/`)

| Script | Función |
|--------|---------|
| `F0_reannotation_v3.py` | Anotación v3 |
| `F0b_dgea_granular_v2.py` | DGEA L2 |
| `F0c_ccc_liana_granular_v2.py` | LIANA+ L2 |
| `F1_*/F8_novae_*.py/R` | Spatial domains |
| `F2_*` | SVG Consensus |
| `F3_spacia_ccc_validation.py` | Spacia pipeline |
| `F3b_run_all_spacia_jobs.sh` | Runner 30 jobs |
| `F3c_aggregate_spacia_results.py` | Agrega + figuras |
| `F5_pseudotime_slingshot.R` | Slingshot + tradeSeq |
| `F_master_figure.py` | Master Figure 5 paneles |

---

## Entornos

```bash
conda activate xenium_pipeline   # Python: scanpy, spatialdata, liana, novae
conda activate xenium_R_analysis # R: Banksy, nnSVG, slingshot, tradeSeq
conda activate spacia            # Python 3.8 + R 4.5.3 (Spacia Bayesian MIL)
```

**Spacia output:** `Pathway_betas.csv` (pval_adj, criterio validación) + `B_and_FDR.csv`

---

## Convenciones

- SpatialData (Zarr) es el único formato. Nunca convertir.
- `adata.X` = log1p · `adata.layers["counts"]` = raw · `adata.layers["denoised"]` = ResolVI (solo viz)
- Salidas → `human_lung_cancer/results/` con subcarpetas por fase.
- Al modificar zarr: `sdata.delete_element_from_disk(key)` → `sdata.write_element(key)`.

---

## Manuscrito (Nature Communications)

**Spec:** `docs/superpowers/specs/2026-05-04-manuscript-design.md` ✅ aprobado  
**Directorio:** `manuscript/` (por crear)  
**Título:** *A multi-tool consensus framework for in situ spatial transcriptomics reveals spatially organized myeloid immunosuppression and validated cell–cell communication hubs in human lung adenocarcinoma*  
**Autores:** Miguel Ángel Díaz-Campos (1°) · Alfredo de Jesús Rodríguez Gómez (correspondencia)  
**Figuras listas:** `human_lung_cancer/results/figures/manuscript/` (9 PNGs) + pre-existentes  
**REGLA CRÍTICA:** Nunca usar códigos internos en el manuscrito (F1, F2, v3, L2, etc.)

## Reporte técnico (interno)

`Reporte_tecnico/Reporte_tecnico.tex` — borrador interno en español; secundario al manuscrito.
```bash
cd Reporte_tecnico && pdflatex -interaction=nonstopmode Reporte_tecnico.tex
```
