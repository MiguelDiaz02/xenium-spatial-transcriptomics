#!/usr/bin/env python
"""End-to-end smoke test for P1 multi-sample infrastructure.

Builds a synthetic TBDs mini-cohort (2 TMAs, 10 donors, 20 cores, 5000 cells)
and runs the M1+M2 pipeline as a real subprocess chain — same CLI you'd run
on production data, just pointed at the smoke fixture.

Pipeline stages exercised:
    [1] concat_samples       (REAL)
    [2] scvi_integrate       (OPTIONAL — skipped by default; use --with-scvi)
    [3] pseudobulk_decoupler (REAL — PyDESeq2)
    [4] pseudobulk_consensus (REAL on Python-only side; muscat R-side skipped)
    [5] LIANA+ multisample   (REAL — n_perms=100 for speed)
    [6] M3-M5 dry-runs       (REAL — verifies plan.json output)

Each step validates expected output files and prints a one-line summary.

Run with:
    cd proyecto_demo_xenium
    conda run -n xenium_pipeline python tests/smoke_test_p1.py
    conda run -n xenium_pipeline python tests/smoke_test_p1.py --with-scvi
    conda run -n xenium_pipeline python tests/smoke_test_p1.py --keep-outputs
"""
from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests"))
sys.path.insert(0, str(REPO_ROOT / "pipeline" / "scripts"))

from fixtures.synthetic_cohort import build_smoke_cohort  # noqa: E402

log = logging.getLogger("smoke_test_p1")


# ─── helpers ────────────────────────────────────────────────────────────
class SmokeFail(RuntimeError):
    pass


def _run(label: str, cmd: List[str], cwd: Path,
         expect_files: Optional[List[Path]] = None,
         allow_fail: bool = False) -> subprocess.CompletedProcess:
    log.info("─" * 70)
    log.info("[%s] %s", label, " ".join(str(x) for x in cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    dt = time.time() - t0
    if proc.returncode != 0:
        log.error("[%s] ✗ failed in %.1fs (exit %d)", label, dt, proc.returncode)
        log.error("STDOUT:\n%s", proc.stdout[-2000:])
        log.error("STDERR:\n%s", proc.stderr[-2000:])
        if not allow_fail:
            raise SmokeFail(label)
        return proc
    log.info("[%s] ✓ ok in %.1fs", label, dt)
    if expect_files:
        for f in expect_files:
            if not f.exists():
                raise SmokeFail(f"{label}: expected output missing: {f}")
            log.info("    output: %s (%.1f KB)", f, f.stat().st_size / 1024)
    return proc


def _validate_cohort_h5ad(h5ad: Path) -> dict:
    import anndata as ad
    a = ad.read_h5ad(h5ad)
    summary = {
        "n_obs": int(a.n_obs), "n_vars": int(a.n_vars),
        "n_samples": int(a.obs["sample_id"].nunique()),
        "n_donors": int(a.obs["subject_id"].nunique()),
        "conditions": sorted(a.obs["condition"].astype(str).unique()),
        "organs": sorted(a.obs["organ"].astype(str).unique()),
    }
    log.info("  cohort summary: %s", summary)
    return summary


def _count_tsvs(root: Path) -> int:
    return sum(1 for _ in root.rglob("*.tsv"))


# ─── pipeline ───────────────────────────────────────────────────────────
def smoke_test(with_scvi: bool, keep_outputs: bool) -> bool:
    fixture_root = REPO_ROOT / "tests" / "fixtures" / "_smoke"

    # ─── [0] Build fixture ──────────────────────────────────────────────
    log.info("═" * 70)
    log.info("[0] Building synthetic cohort under %s", fixture_root)
    log.info("═" * 70)
    paths = build_smoke_cohort(fixture_root, project_root=REPO_ROOT)
    cohort_yaml = paths["cohort_yaml"]
    out_root = fixture_root / "cohort_results"
    out_root.mkdir(parents=True, exist_ok=True)

    python = sys.executable

    # ─── [1] concat_samples ─────────────────────────────────────────────
    cohort_h5ad = out_root / "cohort.h5ad"
    _run(
        "concat_samples",
        [
            python, "pipeline/scripts/integration/concat_samples.py",
            "--cohort", str(cohort_yaml),
            "--outdir", str(out_root),
            "--require-all",
        ],
        cwd=REPO_ROOT,
        expect_files=[cohort_h5ad, out_root / "cohort_meta.json"],
    )
    summary = _validate_cohort_h5ad(cohort_h5ad)
    assert summary["n_obs"] == 5000
    assert summary["n_donors"] == 10
    assert summary["n_samples"] == 20  # 10 donors × 2 cores
    assert set(summary["organs"]) == {"lung", "liver"}
    assert set(summary["conditions"]) >= {"control", "fibrotic_TBD"}

    # ─── [2] scvi_integrate (optional) ─────────────────────────────────
    cohort_integrated = out_root / "cohort_integrated.h5ad"
    if with_scvi:
        _run(
            "scvi_integrate",
            [
                python, "pipeline/scripts/integration/scvi_integrate.py",
                "--input", str(cohort_h5ad),
                "--outdir", str(out_root),
                "--cohort", str(cohort_yaml),
                "--no-gpu",
                "--max-epochs", "20",
            ],
            cwd=REPO_ROOT,
            expect_files=[cohort_integrated],
        )
    else:
        log.info("[scvi_integrate] skipped (re-run with --with-scvi to exercise)")
        # For downstream stages, use cohort.h5ad directly
        cohort_integrated = cohort_h5ad

    # ─── [3] pseudobulk_decoupler ──────────────────────────────────────
    pb_out = out_root / "pseudobulk_decoupler"
    _run(
        "pseudobulk_decoupler",
        [
            python, "pipeline/scripts/integration/pseudobulk_decoupler.py",
            "--input", str(cohort_integrated),
            "--cohort", str(cohort_yaml),
            "--outdir", str(pb_out),
            "--celltype-col", "cell_type_L2",
            "--min-cells", "5",
            "--min-donors", "2",
        ],
        cwd=REPO_ROOT,
    )
    n_de = _count_tsvs(pb_out)
    if n_de == 0:
        raise SmokeFail("pseudobulk_decoupler produced 0 DE TSVs")
    log.info("    %d DE TSVs written", n_de)

    # ─── [4] pseudobulk_consensus ──────────────────────────────────────
    # We don't run muscat here (requires R env); consensus runs on the
    # decoupler outputs alone — exercises the join logic and Fisher math.
    muscat_stub = out_root / "pseudobulk_muscat_stub"
    muscat_stub.mkdir(exist_ok=True)
    cons_out = out_root / "pseudobulk_consensus"
    _run(
        "pseudobulk_consensus",
        [
            python, "pipeline/scripts/integration/pseudobulk_consensus.py",
            "--decoupler-dir", str(pb_out),
            "--muscat-dir", str(muscat_stub),
            "--outdir", str(cons_out),
            "--skip-venn",
        ],
        cwd=REPO_ROOT,
    )
    n_cons = _count_tsvs(cons_out)
    if n_cons == 0:
        raise SmokeFail("pseudobulk_consensus produced 0 TSVs")
    log.info("    %d consensus TSVs", n_cons)

    # ─── [5] LIANA+ multisample ────────────────────────────────────────
    liana_out = out_root / "ccc_liana_multisample"
    _run(
        "LIANA+ multisample",
        [
            python, "pipeline/scripts/integration/F0c_ccc_liana_multisample.py",
            "--input", str(cohort_integrated),
            "--cohort", str(cohort_yaml),
            "--outdir", str(liana_out),
            "--celltype-col", "cell_type_L2",
            "--n-perms", "100",
            "--min-cells", "30",
        ],
        cwd=REPO_ROOT,
    )
    n_liana = _count_tsvs(liana_out)
    if n_liana == 0:
        raise SmokeFail("LIANA produced 0 tables")
    log.info("    %d LIANA TSVs", n_liana)

    # ─── [6] M3-M5 dry-runs ────────────────────────────────────────────
    # spacia_meta dry-run requires a pairs.tsv; build a minimal one
    pairs_tsv = out_root / "spacia_pairs.tsv"
    pairs_tsv.write_text("sender\treceiver\nMacrophage_M1\tFibroblast\n")
    _run(
        "spacia_meta dry-run",
        [
            python, "pipeline/scripts/integration/spacia_meta.py",
            "--cohort", str(cohort_yaml),
            "--input", str(cohort_h5ad),
            "--pairs", str(pairs_tsv),
            "--outdir", str(out_root / "spacia_meta"),
            "--organ", "lung",
        ],
        cwd=REPO_ROOT,
        expect_files=[out_root / "spacia_meta" / "plan.json"],
    )
    _run(
        "novae_crosssample dry-run",
        [
            python, "pipeline/scripts/integration/F1_novae_crosssample.py",
            "--cohort", str(cohort_yaml),
            "--outdir", str(out_root / "spatial_domains"),
            "--mode", "novae",
        ],
        cwd=REPO_ROOT,
        expect_files=[out_root / "spatial_domains" / "plan.json"],
    )

    # ─── Done ───────────────────────────────────────────────────────────
    log.info("═" * 70)
    log.info("SMOKE TEST PASSED ✓")
    log.info("    cohort.h5ad:        %d cells, %d donors, %d cores",
             summary["n_obs"], summary["n_donors"], summary["n_samples"])
    log.info("    DE tables:          %d (decoupler) + %d (consensus)",
             n_de, n_cons)
    log.info("    CCC tables:         %d (LIANA+)", n_liana)
    log.info("    plans (M3/M5):      OK")
    log.info("    outputs under:      %s", out_root)
    log.info("═" * 70)

    if not keep_outputs:
        log.info("cleaning up %s", fixture_root)
        shutil.rmtree(fixture_root, ignore_errors=True)
    else:
        log.info("outputs preserved under %s (--keep-outputs)", fixture_root)

    return True


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--with-scvi", action="store_true",
                   help="run scVI integration (CPU; slow on this fixture)")
    p.add_argument("--keep-outputs", action="store_true",
                   help="don't delete tests/fixtures/_smoke at the end")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    try:
        smoke_test(with_scvi=args.with_scvi, keep_outputs=args.keep_outputs)
    except SmokeFail as e:
        log.error("SMOKE TEST FAILED at step: %s", e)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
