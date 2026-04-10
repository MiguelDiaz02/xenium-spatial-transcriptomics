"""
Rule 09 — Downstream Analyses
Doublet detection (Scrublet) and pseudobulk DE (PyDESeq2).
Pseudobulk requires ≥3 biological samples; disabled for single-sample demos.
"""

if config["run"]["downstream"]:
    rule downstream:
        input:
            done  = str(OUTDIR / "08_spatial.done"),
            sdata = SDATA,
        output:
            done = str(OUTDIR / "09_downstream.done"),
        params:
            doublet_detection      = config["downstream"]["doublet_detection"],
            pseudobulk             = config["downstream"]["pseudobulk"],
            pseudobulk_group_key   = config["downstream"]["pseudobulk_group_key"],
            pseudobulk_sample_key  = config["downstream"]["pseudobulk_sample_key"],
        log:
            str(OUTDIR / "logs" / "09_downstream.log"),
        threads: 4
        script:
            "../scripts/09_downstream.py"
else:
    rule downstream_skip:
        input:
            done = str(OUTDIR / "08_spatial.done"),
        output:
            done = str(OUTDIR / "09_downstream.done"),
        shell:
            "echo '[09_downstream] Skipped (run.downstream: false)' && touch {output.done}"
