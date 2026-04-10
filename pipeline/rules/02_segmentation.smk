"""
Rule 02 — Cell Re-segmentation
Runs the auto-selected method (Cellpose or Baysor) on the ingested SpatialData.
Optionally generates a comparison report against XOA and the alternative method.
SEG_METHOD is resolved at DAG-build time in the Snakefile.
"""

rule segmentation:
    input:
        done  = str(OUTDIR / "01_ingest.done"),
        sdata = SDATA,
    output:
        done  = str(OUTDIR / "02_segmentation.done"),
    params:
        method    = SEG_METHOD,
        cellpose  = config["segmentation"]["cellpose"],
        baysor    = config["segmentation"]["baysor"],
        he_image  = config["data"]["he_image"],
    log:
        str(OUTDIR / "logs" / "02_segmentation.log"),
    threads: 8
    resources:
        gpu = 1,
    script:
        "../scripts/02_segmentation.py"


rule compare_segmentation:
    input:
        done  = str(OUTDIR / "02_segmentation.done"),
        sdata = SDATA,
    output:
        report = str(OUTDIR / "02_segmentation_comparison.html"),
    params:
        seg_method = SEG_METHOD,
    log:
        str(OUTDIR / "logs" / "02_compare_segmentation.log"),
    threads: 4
    script:
        "../scripts/02b_compare_segmentation.py"
