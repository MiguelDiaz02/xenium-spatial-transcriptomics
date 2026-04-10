"""
Rule 07 — ResolVI Denoising (optional, GPU)
Separates true expression from spatial contamination and ambient RNA.
Output stored in adata.layers["denoised"]; raw counts preserved in layers["counts"].
Skipped entirely if run.denoising_resolvi: false.
"""

if config["run"]["denoising_resolvi"]:
    rule denoising:
        input:
            done  = str(OUTDIR / "06_annotation.done"),
            sdata = SDATA,
        output:
            done  = str(OUTDIR / "07_denoising.done"),
        params:
            max_epochs   = config["denoising"]["max_epochs"],
            batch_size   = config["denoising"]["batch_size"],
            accelerator  = config["denoising"]["accelerator"],
        log:
            str(OUTDIR / "logs" / "07_denoising.log"),
        resources:
            gpu = 1,
        threads: 4
        script:
            "../scripts/07_denoising.py"
else:
    rule denoising_skip:
        input:
            done = str(OUTDIR / "06_annotation.done"),
        output:
            done = str(OUTDIR / "07_denoising.done"),
        shell:
            "echo '[07_denoising] Skipped (denoising_resolvi: false)' && touch {output.done}"
