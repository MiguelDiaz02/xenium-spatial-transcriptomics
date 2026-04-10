"""
Rule 01 — Ingest
Extracts Xenium ZIP bundle and converts to SpatialData Zarr store.
Registers H&E image using the alignment CSV if available.
"""

rule ingest:
    input:
        zip      = config["data"]["xenium_zip"],
        he_image = config["data"]["he_image"],
        he_align = config["data"]["he_alignment"],
    output:
        done  = str(OUTDIR / "01_ingest.done"),
        sdata = directory(SDATA),
    params:
        extracted_dir = config["data"]["xenium_extracted_dir"],
    log:
        str(OUTDIR / "logs" / "01_ingest.log"),
    threads: 4
    script:
        "../scripts/01_ingest.py"
