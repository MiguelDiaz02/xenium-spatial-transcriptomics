"""Integration figure — master 5-panel figure combining all families."""

rule integration_master_figure:
    input:
        # F1 spatial domains
        consensus_domains = ADIR / "F1_spatial_domains/consensus_domains.csv",
        # F2 SVG consensus
        svg = ADIR / "F2_svg_consensus/svg_consensus_table.csv",
        # F3 CCC (Phase 2B + bivariate + spacia)
        ccc = ADIR / "F3_ccc_validation/liana_x_spacia_consensus.csv",
        # F4 Niche-DE
        niche_de = ADIR / "F4_niche_de/niche_de_genes.csv",
        # F5 Pseudotime
        pseudotime = ADIR / "F5_pseudotime/slingshot_pseudotime.csv",
        gam = ADIR / "F5_pseudotime/tradeseq_gam_results.csv",
        # F5b TLS
        tls = ADIR / "F5b_tls/tls_objects.csv",
        # F8 Novae embeddings
        novae_embeddings = ADIR / "F8_novae/novae_latent.parquet",
        sdata = SDATA
    output:
        master_fig = FIGDIR / "phase3_integration_master_figure.png",
        master_pdf = FIGDIR / "phase3_integration_master_figure.pdf",
        panel_data = ADIR / "integration/panel_data.json"
    log: "logs/integration_master.log"
    threads: 4
    resources:
        mem_mb = 32000
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/integration_master_figure.py \
            --sdata {input.sdata} \
            --domains {input.consensus_domains} \
            --svg {input.svg} \
            --ccc {input.ccc} \
            --niche-de {input.niche_de} \
            --pseudotime {input.pseudotime} \
            --gam {input.gam} \
            --tls {input.tls} \
            --novae-embeddings {input.novae_embeddings} \
            --outfig {output.master_fig} \
            --outpdf {output.master_pdf} \
            --outdata {output.panel_data} \
            > {log} 2>&1
        """
