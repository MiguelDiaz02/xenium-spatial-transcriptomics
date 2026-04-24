"""
Rule: Step 06b — Immune Cell Subclustering & Granular Annotation
==================================================================
Performs high-resolution subclustering and granular annotation of immune cells.

Dependencies:
  - Input: SpatialData zarr from step 06 (annotation)
  - Output: Enhanced SpatialData with cell_type_immune_granular + validation report

Configuration (from config.yaml):
  - immune_annotation.enabled: true/false
  - immune_annotation.immune_cell_types: list of cell types to subset
  - immune_annotation.leiden_resolution_immune: ~1.2 (high resolution)
  - immune_annotation.immune_markers: dict of subtype → marker genes
  - immune_annotation.validation_threshold: 0.7 (purity cutoff for PASS)

This rule runs ONLY if immune_annotation.enabled == true in config.
"""

rule immune_subclustering:
    input:
        sdata = f"{SDATA}",
        done = f"{OUTDIR}/06_annotation.done"
    output:
        done = f"{OUTDIR}/06b_immune_subclustering.done",
        report = f"{OUTDIR}/06b_immune_subclustering_report.html"
    log:
        "logs/06b_immune_subclustering.log"
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 300  # 5 minutes (subclustering is fast for immune subset)
    conda:
        "envs/xenium_pipeline.yaml"
    params:
        immune_types = config.get("immune_annotation", {}).get("immune_cell_types", [
            "T_cell", "B_cell", "Macrophage", "NK_cell", "Mast_cell"
        ]),
        leiden_res_immune = config.get("immune_annotation", {}).get("leiden_resolution_immune", 1.2),
        n_pcs_immune = config.get("immune_annotation", {}).get("n_pcs_immune", 25),
        min_cluster_size = config.get("immune_annotation", {}).get("min_immune_cluster_size", 50),
        markers = config.get("immune_annotation", {}).get("immune_markers", {}),
        validation_threshold = config.get("immune_annotation", {}).get("validation_threshold", 0.7)
    script:
        "../scripts/06b_immune_subclustering.py"
