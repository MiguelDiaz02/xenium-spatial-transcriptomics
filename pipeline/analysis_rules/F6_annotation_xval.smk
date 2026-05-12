"""F6 — Cell type annotation cross-validation: SingleR vs manual marker-based.

Optional module (config flag f6_annotation_xval.enabled). Validates the
manual + score_genes approach (F0) against reference-based SingleR.

Reference (decided 2026-04-28): HCA Lung Atlas Sikkema 2023 (most recent).

Citations:
- SingleR: Aran D et al. Nat Immunol 20, 163-172 (2019).
- HCA Lung Atlas: Sikkema L et al. Nat Med 29, 1563-1577 (2023).
- Cheng et al. BMC Bioinformatics 26, 22 (2025) — SingleR best for Xenium.
"""

rule F6_singler:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F6_singler.done",
        singler_calls = ADIR / "F6_annotation_xval/singler_predictions.csv",
        crosstab = ADIR / "F6_annotation_xval/manual_vs_singler_crosstab.csv",
        concordance = ADIR / "F6_annotation_xval/concordance_metrics.json"
    log: "logs/F6_singler.log"
    threads: 8
    resources:
        mem_mb = 64000
    conda: "../envs/xenium_R_analysis.yaml"
    params:
        reference = lambda w: ANALYSIS_CFG.get("f6_annotation_xval", {}).get("reference", "Sikkema2023_HCA_Lung")
    shell:
        r"""
        Rscript scripts/analysis/F6_singler.R \
            --sdata {input.sdata} \
            --reference "{params.reference}" \
            --celltype-col cell_type_L2 \
            --outdir {ADIR}/F6_annotation_xval \
            > {log} 2>&1
        touch {output.done}
        """
