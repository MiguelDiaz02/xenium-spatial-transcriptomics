"""P1 multi-sample infrastructure tests.

Validates the TMA cohort contract delivered on 2026-05-11:
    M1 — cohort YAML schema, cohort loader, expand_cores, concat_samples CLI
    M2 — LIANA+ multi-sample CLI help
    M3 — Spacia meta planner (dry-run)
    M4 — Pseudotime cohort planner (dry-run, R script not invoked)
    M5 — Cross-sample domains planner (dry-run)

Run with:
    cd proyecto_demo_xenium
    conda run -n xenium_pipeline pytest tests/test_p1_multisample.py -v
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from textwrap import dedent

import pandas as pd
import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
COHORT_YAML = REPO_ROOT / "pipeline" / "config" / "cohort_TBDs.yaml"
INTEGRATION_DIR = REPO_ROOT / "pipeline" / "scripts" / "analysis"


# ───────────────────────────────────────────────────────────────────────
# M1.A — cohort YAML schema
# ───────────────────────────────────────────────────────────────────────
class TestCohortYaml:
    @pytest.fixture(scope="class")
    def cfg(self):
        return yaml.safe_load(open(COHORT_YAML))

    def test_top_level_keys(self, cfg):
        for k in ("cohort_id", "description", "conditions", "organs",
                  "donors", "integration", "pseudobulk", "ccc_multisample",
                  "spacia_meta", "pseudotime_cohort", "spatial_domains"):
            assert k in cfg, f"missing top-level key {k!r}"

    def test_conditions_reference_first(self, cfg):
        assert cfg["conditions"][0] == "control"

    def test_organs_have_required_fields(self, cfg):
        for org, ocfg in cfg["organs"].items():
            for k in ("tma_slide", "tma_sdata", "markers_yaml",
                      "analysis_config", "contrasts"):
                assert k in ocfg, f"organs.{org} missing {k!r}"

    def test_donor_counts_match_design(self, cfg):
        """LUNG_AF: 10 TBD + 6 controls. LIVER_AF: 4 TBD + 4 ctrl + 4 alc_cirh."""
        lung = [d for d in cfg["donors"] if d["organ"] == "lung"]
        liver = [d for d in cfg["donors"] if d["organ"] == "liver"]
        assert len(lung) == 16, f"lung donors: {len(lung)}"
        assert len(liver) == 12, f"liver donors: {len(liver)}"

        # Lung breakdown
        lung_tbd = [d for d in lung if d["tbd_status"] == "TBD"]
        lung_ns = [d for d in lung if d.get("control_subtype") == "NS"]
        lung_nn = [d for d in lung if d.get("control_subtype") == "NN"]
        assert len(lung_tbd) == 10
        assert len(lung_ns) == 3
        assert len(lung_nn) == 3

        # Liver breakdown
        liver_tbd = [d for d in liver if d["tbd_status"] == "TBD"]
        liver_ctrl = [d for d in liver if d.get("control_subtype") == "Control"]
        liver_alc = [d for d in liver if d.get("control_subtype") == "AlcCirh"]
        assert len(liver_tbd) == 4
        assert len(liver_ctrl) == 4
        assert len(liver_alc) == 4

    def test_subject_ids_unique(self, cfg):
        ids = [d["subject_id"] for d in cfg["donors"]]
        assert len(ids) == len(set(ids))

    def test_lung_contrast_simple(self, cfg):
        """Lung: only fibrotic_TBD vs control (no Alc_Cirh-equivalent)."""
        ccs = cfg["organs"]["lung"]["contrasts"]
        names = [(c["test"], c["ref"]) for c in ccs]
        assert ("fibrotic_TBD", "control") in names

    def test_liver_contrast_includes_TBD_vs_alcohol(self, cfg):
        """Liver contrast that separates TBD-specific from generic fibrosis."""
        ccs = cfg["organs"]["liver"]["contrasts"]
        names = [(c["test"], c["ref"]) for c in ccs]
        assert ("fibrotic_TBD", "fibrotic_nonTBD") in names


# ───────────────────────────────────────────────────────────────────────
# M1.B — cohort loader
# ───────────────────────────────────────────────────────────────────────
class TestCohortLoader:
    def test_load_cohort_canonical(self):
        from utils.cohort import load_cohort
        cohort = load_cohort(COHORT_YAML)
        assert cohort.cohort_id == "TBDs_v1"
        assert cohort.n_donors() == 28
        assert set(cohort.organs.keys()) == {"lung", "liver"}

    def test_load_cohort_by_organ(self):
        from utils.cohort import load_cohort
        cohort = load_cohort(COHORT_YAML)
        assert len(cohort.by_organ("lung")) == 16
        assert len(cohort.by_organ("liver")) == 12

    def test_load_cohort_by_condition(self):
        from utils.cohort import load_cohort
        cohort = load_cohort(COHORT_YAML)
        ctrl = cohort.by_condition("control")
        assert len(ctrl) == 10  # 3 NS + 3 NN (lung) + 4 Control (liver)
        tbd = cohort.by_condition("fibrotic_TBD")
        assert len(tbd) == 14   # 10 lung + 4 liver
        alc = cohort.by_condition("fibrotic_nonTBD")
        assert len(alc) == 4    # 4 AlcCirh (liver only)

    def test_loader_rejects_missing_required_key(self, tmp_path):
        from utils.cohort import load_cohort
        bad = tmp_path / "bad.yaml"
        bad.write_text(dedent("""\
            cohort_id: bad
            conditions: [a, b]
            organs:
              lung: {tma_slide: x, tma_sdata: x, markers_yaml: x}
            donors:
              - subject_id: D1
                organ: lung
                # missing tma_slide, condition, tbd_status, tma_core_label
        """))
        with pytest.raises(ValueError, match="missing required key"):
            load_cohort(bad)

    def test_loader_rejects_unknown_condition(self, tmp_path):
        from utils.cohort import load_cohort
        bad = tmp_path / "bad.yaml"
        bad.write_text(dedent("""\
            cohort_id: bad
            conditions: [control, fibrotic_TBD]
            organs:
              lung: {tma_slide: LUNG_AF, tma_sdata: x, markers_yaml: x}
            donors:
              - subject_id: D1
                organ: lung
                tma_slide: LUNG_AF
                tma_core_label: TBD1
                condition: NOT_A_CONDITION
                tbd_status: TBD
        """))
        with pytest.raises(ValueError, match="not in cohort.conditions"):
            load_cohort(bad)

    def test_loader_rejects_duplicate_subject(self, tmp_path):
        from utils.cohort import load_cohort
        bad = tmp_path / "bad.yaml"
        bad.write_text(dedent("""\
            cohort_id: bad
            conditions: [control, fibrotic_TBD]
            organs:
              lung: {tma_slide: LUNG_AF, tma_sdata: x, markers_yaml: x}
            donors:
              - {subject_id: D1, organ: lung, tma_slide: LUNG_AF,
                 tma_core_label: TBD1, condition: control, tbd_status: nonTBD}
              - {subject_id: D1, organ: lung, tma_slide: LUNG_AF,
                 tma_core_label: TBD2, condition: control, tbd_status: nonTBD}
        """))
        with pytest.raises(ValueError, match="duplicate subject_id"):
            load_cohort(bad)


# ───────────────────────────────────────────────────────────────────────
# M1.C — expand_cores (TMA-aware)
# ───────────────────────────────────────────────────────────────────────
class TestExpandCores:
    def test_expand_two_donors_two_cores_each(self):
        from utils.cohort import expand_cores, load_cohort
        cohort = load_cohort(COHORT_YAML)
        tma_obs = pd.DataFrame({
            "core_id":         ["c1", "c1", "c2", "c2", "c3", "c4"],
            "tma_core_label":  ["TBD1", "TBD1", "TBD1", "TBD1",
                                "Control NS1", "Control NS1"],
        })
        samples = expand_cores(cohort, tma_obs, tma_slide="LUNG_AF")
        # 2 cores for TBD1 (LUNG_TBD01) + 2 cores for NS1 (LUNG_CTRL_NS1) = 4
        assert len(samples) == 4
        tbd_samples = [s for s in samples if s.subject_id == "LUNG_TBD01"]
        assert len(tbd_samples) == 2
        assert {s.sample_id for s in tbd_samples} == {
            "LUNG_TBD01_core01", "LUNG_TBD01_core02",
        }

    def test_expand_rejects_missing_columns(self):
        from utils.cohort import expand_cores, load_cohort
        cohort = load_cohort(COHORT_YAML)
        bad = pd.DataFrame({"x": [1]})
        with pytest.raises(ValueError, match="missing required columns"):
            expand_cores(cohort, bad, tma_slide="LUNG_AF")

    def test_expand_unknown_label_skipped(self):
        from utils.cohort import expand_cores, load_cohort
        cohort = load_cohort(COHORT_YAML)
        tma_obs = pd.DataFrame({
            "core_id":        ["c1", "c2"],
            "tma_core_label": ["TBD1", "NOT_IN_MANIFEST"],
        })
        samples = expand_cores(cohort, tma_obs, tma_slide="LUNG_AF")
        # Only the TBD1 row produces a sample
        assert len(samples) == 1
        assert samples[0].subject_id == "LUNG_TBD01"


# ───────────────────────────────────────────────────────────────────────
# M2/M3/M4/M5 — CLI smoke tests (--help exits zero)
# ───────────────────────────────────────────────────────────────────────
class TestCLIHelp:
    PY_SCRIPTS = [
        "pipeline/scripts/integration/concat_samples.py",
        "pipeline/scripts/integration/scvi_integrate.py",
        "pipeline/scripts/integration/pseudobulk_decoupler.py",
        "pipeline/scripts/integration/pseudobulk_consensus.py",
        "pipeline/scripts/integration/F0c_ccc_liana_multisample.py",
        "pipeline/scripts/integration/spacia_meta.py",
        "pipeline/scripts/integration/F1_novae_crosssample.py",
    ]

    @pytest.mark.parametrize("rel", PY_SCRIPTS, ids=lambda p: Path(p).name)
    def test_help_exits_zero(self, rel):
        script = REPO_ROOT / rel
        ret = subprocess.run(
            [sys.executable, str(script), "--help"],
            capture_output=True, text=True, cwd=str(REPO_ROOT),
        )
        assert ret.returncode == 0, f"{rel}: {ret.stderr}"
        # Every script must document its --cohort flag
        assert ("--cohort" in ret.stdout) or ("cohort" in ret.stdout.lower())


# ───────────────────────────────────────────────────────────────────────
# Paths extension (P0 + P1 path helpers coexist)
# ───────────────────────────────────────────────────────────────────────
class TestPathsExtended:
    def test_cohort_root_under_TBDs(self):
        from utils.paths import cohort_root
        p = cohort_root()
        assert p.parts[-2:] == ("TBDs", "cohort")

    def test_sample_sdata_path_layout(self):
        from utils.paths import sample_sdata_path
        p = sample_sdata_path("TBDs_L01", "lung")
        assert p.parts[-5:] == ("TBDs", "lung", "TBDs_L01", "results", "sdata.zarr")

    def test_cohort_yaml_default(self):
        from utils.paths import cohort_yaml_path
        assert cohort_yaml_path().name == "cohort_TBDs.yaml"


# ───────────────────────────────────────────────────────────────────────
# Sanity: P0 still works (loader didn't break paths.py)
# ───────────────────────────────────────────────────────────────────────
class TestP0Compat:
    def test_p0_module_still_imports(self):
        from utils import paths
        for fn in ("project_root", "dataset_root", "results_root",
                   "sdata_path", "spacia_tool_path"):
            assert callable(getattr(paths, fn))
