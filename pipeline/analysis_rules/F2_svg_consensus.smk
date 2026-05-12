"""F2 — Spatially Variable Gene consensus: Wilcoxon+Moran's I (existing) + Hotspot + nnSVG.

Multi-tool philosophy: a gene is "robust SVG" if significant in >=2 of 3 methods.

Citations:
- Hotspot: DeTomaso D, Yosef N. Cell Syst 12, 446-456 (2021) — lowest FP rate in Salas 2025 benchmark.
- nnSVG: Weber LM et al. Nat Commun 14, 4059 (2023) — scalable Gaussian process.
- Salas et al. Nat Methods 22, 813-823 (2025) — recommends Hotspot for FP control.
"""

rule F2_hotspot_svg:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F2_svg_consensus/hotspot.done",
        scores = ADIR / "F2_svg_consensus/hotspot_svg_scores.csv"
    log: "logs/F2_hotspot.log"
    threads: 8
    resources:
        mem_mb = 32000
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F2_hotspot_svg.py \
            --sdata {input.sdata} \
            --n_neighbors 30 \
            --model danb \
            --fdr_threshold 0.05 \
            --outdir {ADIR}/F2_svg_consensus \
            > {log} 2>&1
        touch {output.done}
        """


rule F2_nnsvg:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F2_svg_consensus/nnsvg.done",
        scores = ADIR / "F2_svg_consensus/nnsvg_svg_scores.csv"
    log: "logs/F2_nnsvg.log"
    threads: 8
    resources:
        mem_mb = 64000
    conda: "../envs/xenium_R_analysis.yaml"
    shell:
        r"""
        Rscript scripts/analysis/F2_nnsvg.R \
            --sdata {input.sdata} \
            --n_neighbors 10 \
            --n_threads {threads} \
            --outdir {ADIR}/F2_svg_consensus \
            > {log} 2>&1
        touch {output.done}
        """


rule F2_svg_consensus:
    """Cross-tabulate SVG calls across Wilcoxon+Moran's (existing Phase 1A), Hotspot, nnSVG.
    Output: gene, n_methods_significant, robust_flag, top_per_method."""
    input:
        hotspot = ADIR / "F2_svg_consensus/hotspot_svg_scores.csv",
        nnsvg = ADIR / "F2_svg_consensus/nnsvg_svg_scores.csv",
        existing_morans = OUTDIR / "02_biology" / "immune_DE_benchmarked" / "morans_i_all_genes.csv"
    output:
        done = ADIR / "F2_svg_consensus.done",
        consensus = ADIR / "F2_svg_consensus/svg_consensus_table.csv",
        venn_data = ADIR / "F2_svg_consensus/svg_venn_data.csv",
        summary = ADIR / "F2_svg_consensus/svg_consensus_summary.md"
    log: "logs/F2_consensus.log"
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F2_svg_consensus.py \
            --hotspot {input.hotspot} \
            --nnsvg {input.nnsvg} \
            --morans {input.existing_morans} \
            --min-methods 2 \
            --outdir {ADIR}/F2_svg_consensus \
            > {log} 2>&1
        touch {output.done}
        """
