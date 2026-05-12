"""F1 — Spatial domain identification: Novae (foundation model) + Banksy (validation)."""

rule F1_novae_domains:
    """Novae zero-shot domain inference using MICS-Lab/novae-human-0 pretrained.
    Citation: Blampey Q et al. Nat Methods 22, 2539-2550 (2025).
    """
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"   # require F0 to ensure cell_type_L2 is in zarr
    output:
        done = ADIR / "F1_novae_domains.done",
        embeddings = ADIR / "F1_spatial_domains/novae_latent.parquet",
        domains = ADIR / "F1_spatial_domains/novae_domains.csv"
    log: "logs/F1_novae.log"
    threads: 4
    resources:
        mem_mb = 32000,
        gpu = 1
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F1_novae_domains.py \
            --sdata {input.sdata} \
            --pretrained MICS-Lab/novae-human-0 \
            --n-domains-levels 5 10 20 \
            --radius 200 \
            --outdir {ADIR}/F1_spatial_domains \
            > {log} 2>&1
        touch {output.done}
        """


rule F1_export_for_banksy:
    """Export counts.h5 + coords.csv from sdata.zarr for Banksy R script."""
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        counts = ADIR / "F1_spatial_domains/counts.h5",
        coords = ADIR / "F1_spatial_domains/coords.csv"
    log: "logs/F1_export.log"
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F1_export_for_banksy.py \
            --sdata {input.sdata} \
            --outdir {ADIR}/F1_spatial_domains \
            > {log} 2>&1
        """


rule F1_banksy_domains:
    """Banksy spatial clustering (Nat Genet 2024 best-in-Salas-benchmark).
    Citation: Singhal V et al. Nat Genet 56, 431-441 (2024).
    """
    input:
        counts = ADIR / "F1_spatial_domains/counts.h5",
        coords = ADIR / "F1_spatial_domains/coords.csv"
    output:
        done = ADIR / "F1_banksy_domains.done",
        domains = ADIR / "F1_spatial_domains/banksy_domains.csv"
    log: "logs/F1_banksy.log"
    threads: 8
    resources:
        mem_mb = 32000
    conda: "../envs/xenium_R_analysis.yaml"
    shell:
        r"""
        Rscript scripts/analysis/F1_banksy_domains.R \
            --counts_h5 {input.counts} \
            --coords_csv {input.coords} \
            --lambda 0.2 \
            --k_geom 6 \
            --n_clusters 10 \
            --outdir {ADIR}/F1_spatial_domains \
            > {log} 2>&1
        touch {output.done}
        """


rule F1_consensus:
    """Cross-validate Novae and Banksy domains -> consensus assignment."""
    input:
        novae = ADIR / "F1_spatial_domains/novae_domains.csv",
        banksy = ADIR / "F1_spatial_domains/banksy_domains.csv"
    output:
        consensus = ADIR / "F1_spatial_domains/consensus_domains.csv",
        ari = ADIR / "F1_spatial_domains/consensus_metrics.json"
    log: "logs/F1_consensus.log"
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F1_domains_consensus.py \
            --novae {input.novae} \
            --banksy {input.banksy} \
            --outdir {ADIR}/F1_spatial_domains \
            > {log} 2>&1
        """
