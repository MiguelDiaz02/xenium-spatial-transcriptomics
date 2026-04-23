"""
Rule 10 — Visualization
Publication-quality figures from spatial analysis results.
"""

rule visualize:
    input:
        done  = str(OUTDIR / "09_downstream.done"),
        sdata = SDATA,
    output:
        done    = str(OUTDIR / "10_visualize.done"),
        figures = directory(str(OUTDIR / "10_visualize_figures")),
    log:
        str(OUTDIR / "logs" / "10_visualize.log"),
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../scripts/10_visualize.py"
