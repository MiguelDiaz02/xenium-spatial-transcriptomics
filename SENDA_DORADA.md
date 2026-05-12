# SENDA DORADA — Xenium Lung Cancer Analysis
**Dataset:** 268,034 células × 289 genes | Human lung cancer, Xenium v1  
**Objetivo:** Pipeline publicable en **Nature Communications** (multi-tool consensus por familia analítica)  
**Entornos:** `xenium_pipeline` (Python) | `xenium_R_analysis` (R) | `spacia` (Fase F)  
**Data:** `human_lung_cancer/results/sdata.zarr` (6.2 GB)  
**Hardware:** RTX 4500 Ada, 18 GB VRAM | Snakemake: sub-Snakefile dedicado (`Snakefile_analysis`)  
**Última actualización:** 2026-04-29 (anotación v3 + DGEA/CCC re-corridos)

---

## MAPA DE FASES

| Fase | Familia | Herramientas | Estado |
|------|---------|-------------|--------|
| **A** | Snakemake skeleton + Phase 0 re-labeling | Pipeline base + 26 cell types (v3 híbrida) | ✅ COMPLETO |
| **B** | **F8** Novae + **F1** Banksy (Spatial Domains) | Foundation model + Banksy consensus | ✅ COMPLETO |
| **C** | **F2** SVG Consensus | Hotspot + nnSVG + Moran's I | ✅ COMPLETO 2026-04-28 |
| **D** | **F4** Niche-DE / Niche-LR | Mason 2024 (Genome Biology) | ❌ DESCARTADO (289 genes insuf. + cell_type_fine = L1, sin L2 real) |
| **E** | **F5** Spatial Pseudotime | Slingshot + tradeSeq + GAM | ✅ COMPLETO 2026-04-30 |
| **F** | **F3** CCC Validation | LIANA+ bivariate + Spacia | ✅ COMPLETO 2026-04-30 |
| **G** | **F6** Annotation x-val + Master Figure | CellTypist (5 modelos) + validación interna ARI | ✅ COMPLETO 2026-04-29 (v3) |
| — | **F7** Co-expression | hdWGCNA / SCENIC | ❌ DESCARTADO (289 genes insuf.) |

**Análisis previos completados (preceden a la arquitectura de 8 familias — válidos como suplementarios):**
- ✅ **Phase 1A** — DGEA Wilcoxon + Moran's I: 288/289 genes con señal espacial
- ✅ **Phase 2B** — CCC LIANA+ rank_aggregate: 141 interacciones L-R significativas

---

## ENTORNOS CONDA

| Entorno | Activación | Contiene |
|---------|-----------|---------|
| `xenium_pipeline` | `conda activate xenium_pipeline` | scanpy, spatialdata, liana≥1.7.1, novae==1.0.5, hotspotsc==1.1.3 |
| `xenium_R_analysis` | `conda activate xenium_R_analysis` | Banksy (Bioc), nnSVG (Bioc), Slingshot (Bioc), tradeSeq (Bioc), SingleR+celldex (Bioc), NicheDE (devtools), rpy2, anndata2ri |
| `spacia` | `conda activate spacia` | Python 3.8 + R≥4.2 + Spacia (Bayesian MIL) — aislado por incompatibilidad de versiones |

**Yamls en:** `pipeline/envs/`

---

## DATOS & ANOTACIÓN CELULAR

- **268,034 células** | **289 genes** | 0% mitocondrial (tejido fresco, excelente calidad)
- **Re-labeling F0 completado** — sistema jerárquico L1/L2/L3:

### Anotación granular v3 en sdata.zarr (activa desde 2026-04-29)

**sdata.zarr ACTIVO** = `human_lung_cancer/results/sdata.zarr` (anotación **v3** — canónica)
**Backup obsoleto** = `human_lung_cancer/results/sdata_DESACTUALIZADO_NO_USAR.zarr` (NO usar — solo 7 tipos L1)

#### Por qué v3 (y no v2)

La v2 usaba únicamente `sc.tl.score_genes` con marcadores del panel → ARI=0.115 vs `cell_type_immune_granular` (subclustering computacional real). CD8_T_cell era clasificado como cDC1 (31%), Plasma→Unassigned (37%), M1 confundido con M2. Inaceptable para análisis downstream.

La v3 es un **enfoque híbrido anclado en los datos computacionales existentes:**

| Fuente | N células | % | Descripción |
|--------|----------|---|-------------|
| `cell_type_immune_granular` | 106,740 | 39.8% | Subclustering Leiden sobre células inmunes — fuente más confiable |
| Leiden directo (clusters puros ≥70%) | 57,650 | 21.5% | Leiden 2→Epithelial_general (92%), 6→Smooth_muscle (88%), 8→Epithelial_general |
| `sc.tl.score_genes` (no-inmunes restantes) | 103,644 | 38.7% | Marcadores INP-Pulmón 2025, solo para tipos no-inmunes |

**Validación post-v3:**

| Comparación | ARI v2 | ARI v3 | Interpretación |
|---|---|---|---|
| L2 vs immune_granular | 0.115 | **0.742** | +543% — concordancia excelente |
| L2 vs Leiden | 0.263 | 0.293 | Ligera mejora |
| L2 vs leiden_immune | 0.116 | 0.194 | Moderada mejora |

**Modelos de referencia utilizados (CellTypist — descartados por cobertura insuficiente):**
Se corrieron 5 modelos CellTypist para validación cruzada: `Human_Lung_Atlas`, `Human_PF_Lung`, `Human_IPF_Lung`, `Immune_All_High`, `Immune_All_Low`. Todos colapsaron a 4-6 tipos amplios (panel de 289 genes << los ~20,000 genes requeridos por los clasificadores). **La validación con CellTypist no es apropiada para paneles dirigidos Xenium v1.** La validación interna (ARI vs subclustering computacional) es la métrica relevante.

**Anotación jerárquica L1 / L2 / L3 (v3):**

| Columna | Categorías | Descripción |
|---------|------------|-------------|
| `cell_type_L1` | **14** | Broad: Epithelial_alveolar/airway/general, Stromal, T_cell, Macrophage, Endothelial, Dendritic_cell, Tumor, Plasma_cell, Monocyte, B_cell, NK_cell, Unassigned |
| `cell_type_L2` | **26** | Granular: AT1, B_memory/naive, CD4_T_helper, CD8_T_cytotoxic/exhausted, Ciliated, DC_mature, Endothelial_blood/lymphatic, Epithelial_general, Macrophage_M1/M2, Monocyte_classical, NK_cytotoxic, Pericyte, Plasma, Plasmablast, Smooth_muscle, Treg, Tumor_proliferating/resting, cDC1, cDC2, pDC |
| `cell_type_L3` | **7** | Estados: Antigen_presenting, Cytotoxic_effector, Exhaustion, Immune_regulatory, Pro_inflammatory, Proliferating, Steady_state |
| `cell_type` | 14 | Alias de L1 (compatibilidad) |
| `cell_type_fine` | 26 | Alias de L2 (compatibilidad) |
| `cell_type_L1_v3` / `cell_type_L2_v3` / `cell_type_L3_v3` | — | Columnas v3 originales (conservadas) |
| `annotation_source_v3` | categórica | immune_granular / leiden_X_direct / score_nonimmune |

**Script de anotación:** `pipeline/scripts/analysis/F0_reannotation_v3.py`
**Validación interna:** `pipeline/scripts/analysis/F6b_internal_annotation_validation.py`
**Validación CellTypist (documentación):** `pipeline/scripts/analysis/F6_celltypist_validation.py`

**Distribución L2 (% de 268,034 células — v3):**
- Epithelial_general 22.1% · CD4_T_helper 11.7% · Endothelial_blood 7.9% · Smooth_muscle 7.3%
- Macrophage_M2 7.0% · Treg 6.1% · AT1 5.1% · Unassigned 4.8% · Pericyte 4.7%
- Tumor_resting 3.3% · Monocyte_classical 3.1% · CD8_T_cytotoxic 2.3% · Ciliated 2.3%
- B_naive 2.2% · NK_cytotoxic 2.0% · Endothelial_lymphatic 1.7% · B_memory 1.6%
- Macrophage_M1 1.5% · Tumor_proliferating 1.1% · cDC1 0.7% · cDC2 0.5%
- CD8_T_exhausted 0.5% · Plasmablast 0.2% · pDC 0.1% · Plasma 0.1% · DC_mature 0.04%

---

## FASE C — SVG Consensus ✅ COMPLETO (2026-04-28)

**Objetivo:** Identificar Spatially Variable Genes robustos por consenso de ≥2 métodos independientes.  
**Benchmark:** Salas 2025 → Hotspot (99.2% TPR) + nnSVG (95.3% TPR) + Moran's I (92.1% TPR)  
**Criterio:** gene "robusto" = significativo en ≥2 de 3 métodos

### Resultados
- **288/289 genes (99.7%) = SVG consenso** (≥2 de 3 métodos)
- **268/289 genes (92.7%) = SVGs en los 3 métodos simultáneamente** (máxima robustez)
- **20 genes en exactamente 2 métodos** (SVGs débiles)
- **1 gene no-SVG** (1/289, excluido del consenso)

| Método | Genes significativos | Umbral |
|--------|---------------------|--------|
| Hotspot | 289/289 (100%) | FDR < 0.05 |
| nnSVG (10k subsample) | 268/289 (92.7%) | padj < 0.05 |
| Moran's I (Phase 1A) | 288/289 (99.7%) | FDR < 0.05 |

### Top SVGs consenso (top 5 por Hotspot C)
`CDH1` (C=0.637) · `EPCAM` (C=0.624) · `AGR3` (C=0.619) · `MUC1` (C=0.595) · `MYH11` (C=0.563)

### Scripts (4 total)
| Script | Estado | Entorno |
|--------|--------|---------|
| `F2_export_for_nnsvg.py` | ✅ Ejecutado | `xenium_pipeline` |
| `F2_hotspot_svg.py` | ✅ Ejecutado (29.5s) | `xenium_pipeline` |
| `F2_nnsvg.R` | ✅ Ejecutado (549.6s) | `xenium_R_analysis` |
| `F2_svg_consensus.py` | ✅ Ejecutado (0.5s) | `xenium_pipeline` |

### Outputs → `results/03_phase3_spatial/F2_svg_consensus/`
`hotspot_svg_scores.csv` · `nnsvg_svg_scores.csv` · `svg_consensus_table.csv` · `svg_all_genes.csv` · `svg_venn_data.json` · `svg_venn.png` · `svg_consensus_heatmap.png`

---

## FASE D — Niche-DE / Niche-LR ❌ DESCARTADO

**Razones de descarte (revisión sistemática 2026-04-29):**

### 1. 289 genes insuficientes — misma razón que F7
NicheDE necesita ≥1000 genes con counts adecuados para modelar expresión contexto-dependiente.
El panel Xenium v1 fue diseñado para identificación de tipos celulares, no para DE contextual.
El método mismo lo confirmó: *"Less than 1000 genes pass"* al ejecutarlo.

### 2. `cell_type_fine` = `cell_type` — no existe L2 en sdata.zarr
Tras inspección del objeto (`python -c "print(adata.obs['cell_type_fine'].unique())"`) se verificó que:
- `cell_type_fine` es **idéntica** a `cell_type` — ambas tienen los mismos 7 tipos
- Los 24 tipos L2 documentados (AT1, B_memory, Ciliated, etc.) **no fueron guardados en sdata.zarr**
- La columna detallada real es `cell_type_immune_granular` (11 tipos, solo células inmunes)

### 3. Escala incompatible
- Diseñado para Visium (~3-5k células), escalado a 268k requiere subsample al 19% (50k)
- Con subsample: niche × celltype combinations esparce demasiado los datos
- 4 bugs críticos para ejecutarlo confirman que la herramienta no fue diseñada para Xenium-scale

### Corrección a la documentación
El SENDA DORADA documentaba "24 tipos L2" que **nunca fueron ejecutados ni guardados**.
Lo que realmente existe: L1=7 tipos + immune_granular=11 tipos. El F0 re-labeling con 24 tipos fue diseño, no ejecución.

### Analogía con F7
| F7 hdWGCNA | F4 NicheDE |
|------------|------------|
| 289 genes → no co-expression network | 289 genes → no niche-DE signal |
| Requiere ≥3k genes | Requiere ≥1k genes |
| ❌ Descartado | ❌ Descartado |

**Aplicabilidad futura:** NicheDE es válido para Xenium Prime (5k genes, panel INP propuesto). Mantener scripts como referencia para esa aplicación.

---

## FASE E — Spatial Pseudotime ⏳

**Pregunta:** ¿Cómo cambia el microambiente tumoral a lo largo del gradiente tumor→estroma?  
**Herramientas:** Slingshot + tradeSeq + GAM  
**Template:** Vannan 2025, *Nature Genetics* (IPF lung) — adaptado de fibrosis → cáncer  
**Anchor confirmado:** distancia euclidiana al centroide tumoral (en lugar de % patología IPF)  
**Input requerido:** embeddings Novae 64d (Fase B ✅)  
**Entorno:** `xenium_R_analysis`  
**Novedad:** primera aplicación de este approach en Xenium lung cancer

---

## FASE F — CCC Validation ⏳

**Objetivo:** Validar top 30 hits de LIANA+ con Spacia (reduce FP 50% vs CellChat/CellPhoneDB)  
**Herramientas:** LIANA+ bivariate (cosine, n_perms=1000, radius=200µm) + Spacia (MCMC 5000 iter)  
**Entornos:** `xenium_pipeline` (LIANA+) + `spacia` (Spacia)  
**Input requerido:** Phase 2B LIANA+ results ✅  
**Output:** pares "doblemente validados" (aparecen significativos en AMBOS)  
**Cita Spacia:** Zhu J et al. Nat Methods 21, 1830–1842 (2024)

---

## FASE G — Integración Final + Master Figure ⏳

**Figura maestra (5 paneles, 300 DPI, Crameri colormaps):**

| Panel | Contenido | Fuente |
|-------|-----------|--------|
| A | Tissue domains (Novae + Banksy consensus overlay) | Fase B |
| B | SVG consensus heatmap (3 métodos) | Fase C |
| C | CCC cross-validated bubble plot (LIANA+ × Spacia) | Fase F |
| D | Pseudotime trajectory + curvas GAM por tipo celular | Fase E |
| E | Annotation x-val heatmap (SingleR vs manual) | Fase G |

**F6 SingleR (opcional):** cross-validar anotación manual L2 vs referencia HCA Lung (Sikkema 2023). Si concordancia >80%, refuerza credibilidad del re-labeling.

---

## RESULTADOS BIOLÓGICOS PRELIMINARES

### Señal Espacial Global
- **288/289 genes (99.7%)** con autocorrelación espacial significativa (Moran's I, p<0.05, FDR-corrected)
- Mediana Moran's I = 2.10 (clustering fuerte positivo, rango: −0.15 a 3.98)
- 100% de los top-50 marcadores DE por tipo celular tienen señal espacial

### Marcadores DE por tipo celular (top 3 por score Wilcoxon)
| Tipo | Top marcadores |
|------|--------------|
| Endothelial | VWF (322.8), PECAM1 (298.1), ADGRL4 (287.2) |
| Tumor | TOP2A (287.3), CCNB1 (276.5), MKI67 (245.1) |
| Macrophage | CD68 (234.8), AIF1 (195.2), TNF (162.1) |
| Epithelial | EPCAM (251.1), MUC1 (221.1), TSPAN8 (187.3) |
| T_cell | CD3E (203.5), IL7R (173.3), CD2 (170.0) |
| B_cell | CD19 (175.2), MS4A1 (142.5), CD79B (139.1) |
| NK_cell | KLRD1 (156.2), NKG7 (128.9), GNLY (125.4) |

### Comunicación Celular — LIANA+ Phase 2B
- **141 interacciones L-R significativas** (rank_aggregate, 5+ métodos, n=1000 permutaciones)
- Top hits: ADAM17/MUC1, CDH1/EGFR, CD86/CTLA4, **CD274/CD80 (PD-L1/PD-1)**, CD34/SELL
- **Macrophages como hub central:** 7 señales entrantes + 7 salientes
- **Eje PD-1/PD-L1 presente pero débil** → supresión dominada por mecanismos mieloides, no por exhaustion de checkpoint
- **Implicación terapéutica:** targeting M2 + bloqueo diferenciación monocito→M2 probablemente más efectivo que monoterapia anti-PD1

### Dominios Espaciales — Fase B
- **F8 Novae** (`MICS-Lab/novae-human-0`): 51.7s, embeddings 64d, dominios L5/L10/L20 (multi-resolución), 3 Hallmark pathways activadas
- **F1 Banksy**: 531.9s, λ=0.2, k_geom=6 → 14 dominios Leiden
- **Consenso F1+F8:** 94,799 células robustas (35.4% acuerdo), ARI=0.1513
- Biología: nichos inmunes robustos (ambos métodos coinciden) + zonas de transición biológicamente plásticas (desacuerdo = ambigüedad real)

### Nichos espaciales identificados (análisis preliminar)
1. **Macrophage ↔ T_cell** — nicho supresivo (inmunomodulación)
2. **Epithelial ↔ Endothelial** — frontera tisular (angiogénesis)
3. **T_cell ↔ B_cell** — clúster linfoide (activación adaptativa)
4. **Tumor ↔ Macrophage** — microambiente pro-tumoral
5. **NK_cell ↔ T_cell** — nicho citotóxico

---

## NARRATIVA DEL PAPER (Nature Comm)

**Hilo central:** El microambiente tumoral del cáncer de pulmón está organizado espacialmente a múltiples escalas. Usando un pipeline multi-herramienta de consenso aplicamos, por primera vez en datos Xenium de pulmón, el modelo fundacional Novae (Nat Methods 2025) junto con 5 familias analíticas independientes, triangulando cada hallazgo antes de reportarlo.

| Pregunta | Respuesta | Familias |
|----------|-----------|---------|
| ¿Qué hay en el tejido? | 7 tipos celulares + dominios espaciales jerárquicos | F1, F6, F8 |
| ¿Qué genes varían espacialmente? | SVG consensus (Hotspot × nnSVG × Moran's I) | F2 |
| ¿Cómo se comunican las células? | LIANA+ validado independientemente por Spacia | F3 |
| ¿Cómo evoluciona el microambiente hacia el tumor? | Pseudotime espacial con GAMs a lo largo del gradiente | F5 |

**Hallazgo diferencial clave:** Eje PD-1/PD-L1 débil + hub M2 macrophage dominante → el microambiente de este tumor está controlado mieloide, no por exhaustion. Implicación directa para estrategia terapéutica.

**Diferenciadores para revisores:**
1. Novae (Nat Methods Dec 2025) — primero en dataset clínico, primer lab México
2. Triangulación ≥2 métodos independientes por cada familia (filosofía LIANA+ generalizada)
3. CCC dual-validado: LIANA+ rank_aggregate (5 métodos) × Spacia (Bayesian MIL) — reduce FP 50%
4. Pseudotime adaptado de IPF (Vannan 2025, Nat Genetics) → cáncer de pulmón
5. Reproducibilidad: Snakemake DAG completo, conda locked, YAML configurable

---

## ARQUITECTURA TÉCNICA

```
Snakefile                  # pipeline preprocesamiento 01-12 (NO tocar)
Snakefile_analysis         # pipeline análisis F1-F8 (nuevo, separado)
pipeline/
├── analysis_rules/
│   ├── F1_spatial_domains.smk   ✅
│   ├── F2_svg_consensus.smk     ✅
│   ├── F3_ccc_validation.smk    ⏳
│   ├── F4_niche_de.smk          ⏳
│   ├── F5_pseudotime.smk        ⏳
│   ├── F6_annotation_xval.smk   ⏳
│   ├── F8_novae.smk             ✅
│   └── integration_viz.smk      ⏳
├── envs/
│   ├── xenium_pipeline.yaml     ✅
│   ├── xenium_R_analysis.yaml   ✅
│   └── spacia.yaml              ✅
└── scripts/analysis/
    ├── [Phase 1A] week3_01_dgea_benchmarked.py + _viz.py
    ├── [Phase 2B] phase2b_ccc_liana.py + _viz.py
    ├── [Fase B]   F8_novae_embeddings.py, F1_novae_domains.py,
    │              F1_banksy_domains.R, F1_export_for_banksy.py, F1_domains_consensus.py
    ├── [Fase C]   F2_hotspot_svg.py ✅, F2_nnsvg.R ✅, F2_svg_consensus.py ⚠️
    └── [D-G]      pendientes
```

**Outputs válidos actuales:**
```
results/
├── 02_biology/
│   ├── immune_DE_benchmarked/   Phase 1A (17 CSVs)
│   ├── ccc_liana/               Phase 2B (5 CSVs)
│   └── reannotation/            F0 (reannotation_summary.md)
└── 03_phase3_spatial/
    ├── F1_spatial_domains/      Fase B (6 CSVs + 1 JSON)
    ├── F8_novae/                Fase B (5 CSVs + 1 JSON)
    └── F2_svg_consensus/        Fase C ⚠️ pendiente ejecución
```

> **Nota:** Existen outputs en `phase3_tissue_region_classification/`, `phase3_tissue_region_refinement/` y otros directorios de análisis previos (pipeline de 4 semanas original). Son válidos como referencia preliminar pero **no forman parte del pipeline F1-F8 oficial**.

---

## SPEC Y DOCUMENTACIÓN

- `docs/superpowers/specs/2026-04-27-multi-tool-pipeline-spec.md` — spec completo 8 familias ★
- `docs/superpowers/specs/2026-04-27-phase3-spatial-mapping-design.md` — diseño Phase 3 detallado
- `docs/superpowers/specs/2026-04-28-phase0-validation-and-installation.md` — decisiones Phase 0

---

*Cargar al inicio de sesión: `proyecto_demo_xenium/SENDA_DORADA.md`*
