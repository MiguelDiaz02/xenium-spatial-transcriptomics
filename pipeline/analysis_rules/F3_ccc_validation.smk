"""F3 — Cell-cell communication validation:
   LIANA+ bivariate (spatial L/R co-expression) + Spacia (Bayesian MIL).

Phase 2B already completed LIANA+ rank_aggregate. This module adds:
  1) LIANA+ bivariate: spatial L/R cosine local score
  2) Spacia validation of top 30 LIANA+ hits (independent Bayesian framework)

Citations:
- LIANA+: Dimitrov D et al. Nat Cell Biol 26, 1613-1623 (2024).
- Spacia: Zhu J et al. Nat Methods 21, 1830-1842 (2024).
"""

PHASE2B_RESULTS = OUTDIR / "02_biology" / "ccc_liana"

rule F3_liana_bivariate:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F3_liana_bivariate.done",
        bivariate = ADIR / "F3_ccc_validation/spatial_bivariate_lr.csv"
    log: "logs/F3_liana_bivariate.log"
    threads: 8
    resources:
        mem_mb = 32000
    conda: "../envs/xenium_pipeline.yaml"
    shell:
        r"""
        python scripts/analysis/F3_liana_bivariate.py \
            --sdata {input.sdata} \
            --local-name cosine \
            --resource consensus \
            --n-perms 1000 \
            --n-neighbors 15 \
            --radius 200 \
            --nz-prop 0.05 \
            --outdir {ADIR}/F3_ccc_validation \
            > {log} 2>&1
        touch {output.done}
        """


rule F3_spacia_validation:
    """Validate top-30 LIANA+ hits using Spacia (independent Bayesian MIL framework).
    Cross-validation: a pair is "doubly validated" if significant in BOTH methods."""
    input:
        sdata = SDATA,
        liana_significant = PHASE2B_RESULTS / "liana_significant_interactions.csv",
        f0_done = ADIR / "F0_reannotation.done"
    output:
        done = ADIR / "F3_spacia_validation.done",
        spacia_results = ADIR / "F3_ccc_validation/spacia_validation.csv",
        crossval = ADIR / "F3_ccc_validation/liana_x_spacia_consensus.csv"
    log: "logs/F3_spacia.log"
    threads: 8
    resources:
        mem_mb = 64000
    conda: "../envs/spacia.yaml"
    shell:
        r"""
        python scripts/analysis/F3_spacia_validation.py \
            --sdata {input.sdata} \
            --liana-pairs {input.liana_significant} \
            --top-n 30 \
            --spacia-repo /home/mdiaz/tools/Spacia \
            --n-mcmc-iters 5000 \
            --outdir {ADIR}/F3_ccc_validation \
            > {log} 2>&1
        touch {output.done}
        """
