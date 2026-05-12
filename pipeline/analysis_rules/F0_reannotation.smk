"""F0 — Hierarchical cell type re-annotation (L1/L2/L3 + functional states)."""

rule F0_reannotation:
    input:
        sdata = SDATA,
        markers = "config/cell_markers_lung.yaml"
    output:
        done = ADIR / "F0_reannotation.done",
        l2_scores = ADIR / "F0_reannotation/L2_scores.csv",
        l3_scores = ADIR / "F0_reannotation/L3_scores.csv",
        assignments = ADIR / "F0_reannotation/cell_assignments.csv",
        crosstab = ADIR / "F0_reannotation/L1_vs_L2_crosstab.csv",
        summary = ADIR / "F0_reannotation/reannotation_summary.md"
    log: "logs/F0_reannotation.log"
    threads: 4
    resources:
        mem_mb = 16000
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F0_reannotation.py \
            --sdata {input.sdata} \
            --config {input.markers} \
            --outdir {ADIR}/F0_reannotation \
            --write-back \
            > {log} 2>&1
        touch {output.done}
        """
