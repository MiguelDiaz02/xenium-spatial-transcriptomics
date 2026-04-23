"""
Rule 12 — Export & Package
Format conversion and data packaging for external use.
"""

rule export:
    input:
        done  = str(OUTDIR / "11_analysis.done"),
        sdata = SDATA,
    output:
        done    = str(OUTDIR / "12_export.done"),
        exports = directory(str(OUTDIR / "12_exports")),
    log:
        str(OUTDIR / "logs" / "12_export.log"),
    threads: 2
    resources:
        mem_mb = 32000
    script:
        "../scripts/12_export.py"
