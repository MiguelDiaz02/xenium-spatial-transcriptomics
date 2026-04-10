"""
Rule 04 — Preprocessing
Normalization, log-transform, and optional HVG selection.
Raw counts are always preserved in adata.layers["counts"].
"""

rule preprocess:
    input:
        done  = str(OUTDIR / "03_qc.done"),
        sdata = SDATA,
    output:
        done  = str(OUTDIR / "04_preprocess.done"),
    params:
        target_sum = config["preprocess"]["normalization_target_sum"],
        log1p      = config["preprocess"]["log1p"],
        select_hvg = config["preprocess"]["select_hvg"],
        n_hvg      = config["preprocess"]["n_hvg"],
    log:
        str(OUTDIR / "logs" / "04_preprocess.log"),
    threads: 4
    script:
        "../scripts/04_preprocess.py"
