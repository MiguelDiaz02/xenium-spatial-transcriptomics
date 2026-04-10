"""
Rule 08 — Spatial Statistics
Squidpy-based spatial analysis: neighborhood enrichment, Moran's I,
co-occurrence, and Ripley's statistics.
"""

rule spatial:
    input:
        done  = str(OUTDIR / "07_denoising.done"),
        sdata = SDATA,
    output:
        done    = str(OUTDIR / "08_spatial.done"),
        figures = directory(str(OUTDIR / "08_spatial_figures")),
    params:
        n_neighs            = config["spatial"]["n_neighs"],
        coord_type          = config["spatial"]["coord_type"],
        cluster_key         = config["spatial"]["cluster_key"],
        n_perms_enrichment  = config["spatial"]["n_perms_enrichment"],
        autocorr            = config["spatial"]["spatial_autocorr"],
    log:
        str(OUTDIR / "logs" / "08_spatial.log"),
    threads: config["spatial"]["spatial_autocorr"]["n_jobs"]
    script:
        "../scripts/08_spatial.py"
