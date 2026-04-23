"""
Rule 11 — Downstream Analysis
Cell type-specific spatial patterns and statistical analysis.
"""

rule analysis:
    input:
        done  = str(OUTDIR / "10_visualize.done"),
        sdata = SDATA,
    output:
        done    = str(OUTDIR / "11_analysis.done"),
        results = directory(str(OUTDIR / "11_analysis")),
    log:
        str(OUTDIR / "logs" / "11_analysis.log"),
    threads: 4
    resources:
        mem_mb = 16000
    script:
        "../scripts/11_analysis.py"
