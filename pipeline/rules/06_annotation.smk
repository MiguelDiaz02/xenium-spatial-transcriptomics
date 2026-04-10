"""
Rule 06 — Cell Type Annotation
CellTypist (default), manual marker scoring, or scANVI.
The celltypist_model key in config is tissue-specific — swap for liver.
"""

rule annotation:
    input:
        done  = str(OUTDIR / "05_reduction.done"),
        sdata = SDATA,
    output:
        done  = str(OUTDIR / "06_annotation.done"),
    params:
        method            = config["annotation"]["method"],
        celltypist_model  = config["annotation"]["celltypist_model"],
        majority_voting   = config["annotation"]["majority_voting"],
        markers           = config["annotation"]["markers"],
    log:
        str(OUTDIR / "logs" / "06_annotation.log"),
    threads: 4
    script:
        "../scripts/06_annotation.py"
