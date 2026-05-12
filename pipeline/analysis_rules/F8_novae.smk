"""F8 — Novae foundation model: embeddings + spatial pathway analysis.

Drives F1 (domains via Novae assignment) and provides substrate for F5 (Slingshot
on Novae embeddings). Also enables spatial pathway analysis per domain.

Citation: Blampey Q, Benkirane H, Bercovici N et al. Novae: a graph-based
foundation model for spatial transcriptomics data. Nat Methods 22, 2539-2550 (2025).
doi:10.1038/s41592-025-02899-6
Pretrained: MICS-Lab/novae-human-0 (HuggingFace, 1.6k+ downloads as of Dec 2025).
"""

rule F8_novae:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F8_novae.done",
        embeddings = ADIR / "F8_novae/novae_latent.parquet",
        domains_multilevel = ADIR / "F8_novae/novae_domains_multilevel.csv",
        pathway_scores = ADIR / "F8_novae/novae_pathway_scores.csv",
        svg_scores = ADIR / "F8_novae/novae_svg_scores.csv"
    log: "logs/F8_novae.log"
    threads: 4
    resources:
        mem_mb = 32000,
        gpu = 1
    conda: "../envs/xenium_pipeline.yaml"
    params:
        pretrained = "MICS-Lab/novae-human-0",
        n_levels = "5,10,20",
        radius = 200,
        pathway_db = "msigdb_hallmark"
    shell:
        r"""
        python scripts/analysis/F8_novae_embeddings.py \
            --sdata {input.sdata} \
            --pretrained {params.pretrained} \
            --n-levels {params.n_levels} \
            --radius {params.radius} \
            --pathway-db {params.pathway_db} \
            --batch-correct \
            --outdir {ADIR}/F8_novae \
            > {log} 2>&1
        touch {output.done}
        """
