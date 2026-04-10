"""
Rule 05 — Dimensionality Reduction & Clustering
PCA → k-NN neighbors → UMAP → Leiden clustering.
"""

rule reduction:
    input:
        done  = str(OUTDIR / "04_preprocess.done"),
        sdata = SDATA,
    output:
        done  = str(OUTDIR / "05_reduction.done"),
    params:
        n_pcs         = config["reduction"]["n_pcs"],
        n_neighbors   = config["reduction"]["n_neighbors"],
        umap_min_dist = config["reduction"]["umap_min_dist"],
        resolution    = config["reduction"]["leiden_resolution"],
    log:
        str(OUTDIR / "logs" / "05_reduction.log"),
    threads: 8
    script:
        "../scripts/05_reduction.py"
