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
    resources:
        mem_mb = 32000,       # 32 GB for neighborhood enrichment permutations
        runtime = 3600        # 1 hour timeout (neighborhood enrichment is compute-intensive)
    script:
        "../scripts/08_spatial.py"
