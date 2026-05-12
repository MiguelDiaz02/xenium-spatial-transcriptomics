"""F4 — NicheDE: Context-dependent differential expression per spatial niche.

Identifies genes differentially expressed in cell type X specifically when
X is located within spatial domain Y (from Phase B Banksy domains).

Multi-tool philosophy: NicheDE output is validated by cross-referencing with
Phase 2B LIANA+ CCC results (do niche-DE genes participate in L-R interactions?).

Citations:
- Mason JC et al. Genome Biol 25, 107 (2024)  — NicheDE method.
- F1 Banksy domains from Phase B provide the niche context.
"""

rule F4_export_for_nichede:
    """Export counts, coordinates, cell types, and domain labels for R."""
    input:
        sdata   = SDATA,
        domains = ADIR / "F1_spatial_domains/banksy_domains.csv",
        f0_done = ADIR / "F0_reannotation.done",
        fb_done = ADIR / "F1_spatial_domains.done",
    output:
        done         = ADIR / "F4_nichede_export.done",
        counts_mtx   = ADIR / "F4_niche_de/export/counts_matrix.mtx",
        cell_types   = ADIR / "F4_niche_de/export/cell_types.csv",
        cell_domains = ADIR / "F4_niche_de/export/cell_domains.csv",
        coords       = ADIR / "F4_niche_de/export/spatial_coords.csv",
    log: "logs/F4_export.log"
    threads: 4
    resources:
        mem_mb  = 24000,
        runtime = 30,
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F4_export_for_nichede.py \
            --sdata    {input.sdata} \
            --domains  {input.domains} \
            --outdir   {ADIR}/F4_niche_de/export \
            > {log} 2>&1
        touch {output.done}
        """


rule F4_niche_de:
    """Run NicheDE: niche-specific differential expression (R, xenium_R_analysis env)."""
    input:
        export_done  = ADIR / "F4_nichede_export.done",
        counts_mtx   = ADIR / "F4_niche_de/export/counts_matrix.mtx",
        cell_types   = ADIR / "F4_niche_de/export/cell_types.csv",
        cell_domains = ADIR / "F4_niche_de/export/cell_domains.csv",
        coords       = ADIR / "F4_niche_de/export/spatial_coords.csv",
    output:
        done      = ADIR / "F4_niche_de.done",
        sig_genes = ADIR / "F4_niche_de/nichede_significant_genes.csv",
        summary   = ADIR / "F4_niche_de/nichede_summary.json",
    log: "logs/F4_niche_de.log"
    threads: 4
    resources:
        mem_mb  = 48000,
        runtime = 120,
    conda: "../envs/xenium_R_analysis.yaml"
    shell:
        r"""
        Rscript scripts/analysis/F4_niche_de.R \
            --datadir      {ADIR}/F4_niche_de/export \
            --outdir       {ADIR}/F4_niche_de \
            --cell_type_col cell_type \
            --sigma        100 \
            --num_cores    {threads} \
            --C            150 \
            --M            10 \
            > {log} 2>&1
        touch {output.done}
        """


rule F4_niche_de_viz:
    """Visualize NicheDE results: bubble plot, heatmap, spatial map."""
    input:
        nde_done  = ADIR / "F4_niche_de.done",
        sig_genes = ADIR / "F4_niche_de/nichede_significant_genes.csv",
        sdata     = SDATA,
    output:
        done    = ADIR / "F4_niche_de_viz.done",
        bubble  = ADIR / "F4_niche_de/figures/nichede_bubble_plot.png",
        heatmap = ADIR / "F4_niche_de/figures/nichede_gene_heatmap.png",
    log: "logs/F4_viz.log"
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 30,
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F4_niche_de_viz.py \
            --results_dir {ADIR}/F4_niche_de \
            --sdata       {input.sdata} \
            --outdir      {ADIR}/F4_niche_de/figures \
            > {log} 2>&1
        touch {output.done}
        """
