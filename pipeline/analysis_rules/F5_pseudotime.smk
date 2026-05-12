"""F5 — Spatial pseudotime: rank tissue regions by tumor proximity.

Adapts the Vannan 2025 (Nat Genetics) IPF template to cancer:
  - Anchor: distance to tumor centroid (configured by user 2026-04-28)
  - Slingshot trajectory on Novae embeddings
  - tradeSeq GAM fits gene expression along pseudotime per cell type

Citations:
- Slingshot: Street K et al. BMC Genomics 19, 477 (2018).
- tradeSeq: Van den Berge K et al. Nat Commun 11, 1201 (2020).
- Vannan 2025: Vannan A et al. Nat Genet 57, 647-658 (2025).
"""

rule F5_compute_tumor_distance:
    """Step 1: compute distance from each cell to tumor centroid (anchor)."""
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        distances = ADIR / "F5_pseudotime/tumor_distance.csv"
    log: "logs/F5_distance.log"
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F5_tumor_distance.py \
            --sdata {input.sdata} \
            --tumor-celltype "Tumor_proliferating,Tumor_resting" \
            --celltype-col cell_type_L2 \
            --outdir {ADIR}/F5_pseudotime \
            > {log} 2>&1
        """


rule F5_slingshot:
    """Step 2: Slingshot lineage inference on Novae embeddings, anchored by tumor distance."""
    input:
        sdata = SDATA,
        embeddings = ADIR / "F1_spatial_domains/novae_latent.parquet",
        distances = ADIR / "F5_pseudotime/tumor_distance.csv",
        f8_done = ADIR / "F8_novae.done"
    output:
        pseudotime = ADIR / "F5_pseudotime/slingshot_pseudotime.csv",
        lineages = ADIR / "F5_pseudotime/slingshot_lineages.json"
    log: "logs/F5_slingshot.log"
    threads: 4
    conda: "../envs/xenium_R_analysis.yaml"
    shell:
        r"""
        Rscript scripts/analysis/F5_slingshot.R \
            --embeddings {input.embeddings} \
            --distances {input.distances} \
            --start-cluster "stroma_far" \
            --outdir {ADIR}/F5_pseudotime \
            > {log} 2>&1
        """


rule F5_tradeseq_gam:
    """Step 3: tradeSeq GAM fits gene expression ~ pseudotime per cell type."""
    input:
        sdata = SDATA,
        pseudotime = ADIR / "F5_pseudotime/slingshot_pseudotime.csv"
    output:
        done = ADIR / "F5_pseudotime.done",
        gam_results = ADIR / "F5_pseudotime/tradeseq_gam_results.csv",
        plots_dir = directory(ADIR / "F5_pseudotime/gam_plots")
    log: "logs/F5_tradeseq.log"
    threads: 8
    resources:
        mem_mb = 64000
    conda: "../envs/xenium_R_analysis.yaml"
    shell:
        r"""
        Rscript scripts/analysis/F5_tradeseq_gam.R \
            --sdata {input.sdata} \
            --pseudotime {input.pseudotime} \
            --celltype-col cell_type_L2 \
            --n-knots 6 \
            --pval-threshold 0.05 \
            --outdir {ADIR}/F5_pseudotime \
            > {log} 2>&1
        touch {output.done}
        """
