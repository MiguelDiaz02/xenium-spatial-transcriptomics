"""F5b — Tertiary Lymphoid Structures (TLS) detection.

User-requested addition (2026-04-28). TLS are immune aggregates of B + T + DC
+ lymphatic chemokines (CCL19, CCL21, CXCL13). Their density is a known
prognostic marker in lung cancer immunotherapy response.

Strategy:
  1. F0 already computed TLS_signature score (CCL19, CCL21, CXCL13, MS4A1, CD3D, CCR7).
  2. Spatially cluster TLS+ cells (radius-based DBSCAN) -> discrete TLS objects.
  3. Quantify size, cellular composition, distance-to-tumor per TLS.
  4. Correlate TLS density with immune infiltration patterns.
"""

rule F5b_tls:
    input:
        sdata = SDATA,
        f0_done = ADIR / "F0_reannotation.done",
        tumor_dist = ADIR / "F5_pseudotime/tumor_distance.csv"
    output:
        done = ADIR / "F5b_tls.done",
        tls_objects = ADIR / "F5b_tls/tls_objects.csv",
        tls_composition = ADIR / "F5b_tls/tls_composition.csv",
        tls_summary = ADIR / "F5b_tls/tls_summary.md"
    log: "logs/F5b_tls.log"
    threads: 4
    resources:
        mem_mb = 16000
    conda: "../envs/xenium_pipeline.yaml"
    params:
        eps_um = 50,                  # DBSCAN radius
        min_samples = 30,             # minimum cells to form a TLS
        tls_score_thr = 1.5
    shell:
        r"""
        python scripts/analysis/F5b_tls_detection.py \
            --sdata {input.sdata} \
            --tumor-distance {input.tumor_dist} \
            --eps-um {params.eps_um} \
            --min-samples {params.min_samples} \
            --tls-score-threshold {params.tls_score_thr} \
            --outdir {ADIR}/F5b_tls \
            > {log} 2>&1
        touch {output.done}
        """
