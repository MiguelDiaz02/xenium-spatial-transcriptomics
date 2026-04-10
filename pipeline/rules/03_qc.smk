"""
Rule 03 — Quality Control
Computes per-cell QC metrics and filters low-quality cells.
"""

rule qc:
    input:
        done  = str(OUTDIR / "02_segmentation.done"),
        sdata = SDATA,
    output:
        done   = str(OUTDIR / "03_qc.done"),
        report = str(OUTDIR / "03_qc_report.html"),
    params:
        min_transcripts  = config["qc"]["min_transcripts"],
        max_transcripts  = config["qc"]["max_transcripts"],
        min_genes        = config["qc"]["min_genes"],
        max_pct_neg_ctrl = config["qc"]["max_pct_neg_ctrl"],
    log:
        str(OUTDIR / "logs" / "03_qc.log"),
    threads: 4
    script:
        "../scripts/03_qc.py"
