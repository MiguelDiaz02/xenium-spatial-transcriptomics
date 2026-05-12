"""P0 pipeline portability tests.

Verifies the resolution machinery delivered on 2026-05-08:
    B1 — path lifting (utils.paths + env-var overrides)
    B2 — marker YAML externalization + F0 loader
    B3 — threshold parameterization (config_lung.yaml + F3 loader)
    B4 — gene-panel-constrained method exclusions documented
    B5 — Spacia binary resolution via SPACIA_PATH env var

Run with:
    cd proyecto_demo_xenium
    conda run -n xenium_pipeline pytest tests/test_p0_portability.py -v
"""
from __future__ import annotations

import importlib
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "pipeline" / "config" / "config_lung.yaml"
MARKERS_DIR = REPO_ROOT / "pipeline" / "config" / "markers"
MANUSCRIPT_TEX = REPO_ROOT / "manuscript" / "manuscript.tex"
SCRIPTS_DIR = REPO_ROOT / "pipeline" / "scripts" / "analysis"
SHELL_SCRIPTS = [
    SCRIPTS_DIR / "F3b_run_all_spacia_jobs.sh",
    SCRIPTS_DIR / "F3d_retry_rare_celltypes.sh",
]


# ───────────────────────────────────────────────────────────────────────
# B1 — path lifting
# ───────────────────────────────────────────────────────────────────────
class TestB1Paths:
    def test_paths_module_imports(self):
        from utils import paths
        for fn in ("project_root", "dataset_root", "results_root",
                   "sdata_path", "spacia_tool_path"):
            assert callable(getattr(paths, fn)), f"missing {fn}"

    def test_project_root_walks_up_to_repo(self):
        from utils.paths import project_root
        assert project_root() == REPO_ROOT

    def test_project_root_respects_env_override(self, tmp_path, monkeypatch):
        monkeypatch.setenv("XENIUM_PROJECT_ROOT", str(tmp_path))
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        assert paths_mod.project_root() == tmp_path.resolve()

    def test_dataset_default_is_human_lung_cancer(self, monkeypatch):
        monkeypatch.delenv("XENIUM_DATASET", raising=False)
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        assert paths_mod.dataset_root().name == "human_lung_cancer"

    def test_dataset_respects_env_override(self, monkeypatch):
        monkeypatch.setenv("XENIUM_DATASET", "TBDs_lung_sample01")
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        assert paths_mod.dataset_root().name == "TBDs_lung_sample01"

    def test_canonical_sdata_zarr_exists(self):
        from utils.paths import sdata_path
        # Restore default env first by reloading
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        assert paths_mod.sdata_path().exists(), (
            f"canonical sdata.zarr missing at {paths_mod.sdata_path()}"
        )

    def test_no_hardcoded_mdiaz_paths_in_scripts(self):
        """B1 regression: no /home/mdiaz in pipeline/scripts/ or pipeline/config/."""
        offenders = []
        scan_dirs = [REPO_ROOT / "pipeline" / "scripts",
                     REPO_ROOT / "pipeline" / "config"]
        for d in scan_dirs:
            for p in d.rglob("*"):
                if not p.is_file():
                    continue
                if p.suffix in (".py", ".R", ".sh", ".yaml", ".yml"):
                    if "__pycache__" in p.parts:
                        continue
                    try:
                        txt = p.read_text(errors="ignore")
                    except Exception:
                        continue
                    if "/home/mdiaz" in txt:
                        offenders.append(str(p.relative_to(REPO_ROOT)))
        assert not offenders, f"hardcoded /home/mdiaz in: {offenders}"


# ───────────────────────────────────────────────────────────────────────
# B2 — markers externalization
# ───────────────────────────────────────────────────────────────────────
REQUIRED_KEYS = ["markers", "hierarchy_l1", "hierarchy_l3",
                 "immune_granular_to_l2", "nonimmune_l2",
                 "annotation_params"]

MARKER_FILES = ["lung_pilot.yaml", "lung_TBDs.yaml", "liver_TBDs.yaml"]


class TestB2Markers:
    @pytest.mark.parametrize("yaml_name", MARKER_FILES)
    def test_marker_yaml_exists(self, yaml_name):
        assert (MARKERS_DIR / yaml_name).exists()

    @pytest.mark.parametrize("yaml_name", MARKER_FILES)
    def test_marker_yaml_parses_with_required_keys(self, yaml_name):
        cfg = yaml.safe_load(open(MARKERS_DIR / yaml_name))
        for k in REQUIRED_KEYS:
            assert k in cfg, f"{yaml_name}: missing key {k}"

    @pytest.mark.parametrize("yaml_name", MARKER_FILES)
    def test_marker_signatures_have_min_2_markers(self, yaml_name):
        cfg = yaml.safe_load(open(MARKERS_DIR / yaml_name))
        short = [k for k, v in cfg["markers"].items() if len(v) < 2]
        assert not short, f"{yaml_name}: signatures with <2 markers: {short}"

    @pytest.mark.parametrize("yaml_name", MARKER_FILES)
    def test_every_marker_l2_has_l1_l3_mapping(self, yaml_name):
        cfg = yaml.safe_load(open(MARKERS_DIR / yaml_name))
        l2_keys = set(cfg["markers"].keys())
        l1_keys = set(cfg["hierarchy_l1"].keys()) - {"Unassigned"}
        l3_keys = set(cfg["hierarchy_l3"].keys()) - {"Unassigned"}
        missing_l1 = l2_keys - l1_keys
        missing_l3 = l2_keys - l3_keys
        assert not missing_l1, f"{yaml_name}: L2 without L1 mapping: {missing_l1}"
        assert not missing_l3, f"{yaml_name}: L2 without L3 mapping: {missing_l3}"

    def test_lung_pilot_matches_canonical_n_types(self):
        """Pilot YAML must produce the same 25 L2 types as the canonical run."""
        cfg = yaml.safe_load(open(MARKERS_DIR / "lung_pilot.yaml"))
        assert len(cfg["markers"]) == 25

    def test_lung_pilot_leiden_direct_panel_specific(self):
        cfg = yaml.safe_load(open(MARKERS_DIR / "lung_pilot.yaml"))
        # Pilot used clusters 6/8/2 for Smooth_muscle / Epithelial_general.
        assert cfg["leiden_direct"]["6"] == "Smooth_muscle"
        assert cfg["leiden_direct"]["8"] == "Epithelial_general"
        assert cfg["annotation_params"]["use_leiden_direct"] is True

    @pytest.mark.parametrize("yaml_name", ["lung_TBDs.yaml", "liver_TBDs.yaml"])
    def test_tbds_templates_disable_leiden_direct(self, yaml_name):
        """Leiden cluster IDs are pilot-specific; TBDs templates must not reuse them."""
        cfg = yaml.safe_load(open(MARKERS_DIR / yaml_name))
        ld = cfg.get("leiden_direct") or {}
        assert ld == {}, f"{yaml_name} reuses pilot leiden_direct: {ld}"
        assert cfg["annotation_params"]["use_leiden_direct"] is False

    @pytest.mark.parametrize("yaml_name", ["lung_TBDs.yaml", "liver_TBDs.yaml"])
    def test_tbds_templates_have_no_tumor_markers(self, yaml_name):
        """TBDs is fibrotic, not malignant — tumor signatures must be absent."""
        cfg = yaml.safe_load(open(MARKERS_DIR / yaml_name))
        tumor = [k for k in cfg["markers"] if "tumor" in k.lower()]
        assert not tumor, f"{yaml_name}: unexpected tumor markers {tumor}"

    def test_lung_tbds_includes_fibrotic_markers(self):
        """TBDs lung must include myofibroblast + AT2 dysfunction biology."""
        cfg = yaml.safe_load(open(MARKERS_DIR / "lung_TBDs.yaml"))
        assert "Myofibroblast" in cfg["markers"]
        assert "AT2_dysfunctional" in cfg["markers"]
        # Canonical fibrotic axis must be in Myofibroblast signature
        myof = cfg["markers"]["Myofibroblast"]
        for required in ("ACTA2", "COL1A1"):
            assert required in myof, f"Myofibroblast missing {required}"

    def test_liver_tbds_includes_stellate_zonation(self):
        cfg = yaml.safe_load(open(MARKERS_DIR / "liver_TBDs.yaml"))
        assert "HSC_quiescent" in cfg["markers"]
        assert "HSC_activated" in cfg["markers"]
        assert "LRAT" in cfg["markers"]["HSC_quiescent"]   # quiescent hallmark
        # Zonal hepatocytes
        for k in ("Hepatocyte_pericentral", "Hepatocyte_periportal"):
            assert k in cfg["markers"]

    def test_f0_loader_reads_canonical_yaml(self):
        """B2 contract: F0_reannotation_v3.load_markers_config returns required keys."""
        F0 = importlib.import_module("analysis.F0_reannotation_v3")
        cfg = F0.load_markers_config(MARKERS_DIR / "lung_pilot.yaml")
        for k in REQUIRED_KEYS:
            assert k in cfg
        assert len(cfg["markers"]) == 25

    def test_f0_loader_rejects_incomplete_yaml(self, tmp_path):
        """Loader must fail loudly on a YAML missing required sections."""
        F0 = importlib.import_module("analysis.F0_reannotation_v3")
        bad = tmp_path / "bad.yaml"
        bad.write_text("markers: {AT1: [AGER]}\n")  # missing every other key
        with pytest.raises(ValueError, match="missing required keys"):
            F0.load_markers_config(bad)


# ───────────────────────────────────────────────────────────────────────
# B3 — threshold parameterization
# ───────────────────────────────────────────────────────────────────────
class TestB3Thresholds:
    @pytest.fixture(scope="class")
    def config(self):
        return yaml.safe_load(open(CONFIG_PATH))

    def test_qc_max_cell_area(self, config):
        """B3: cell-area threshold must be externalized (was hardcoded)."""
        assert "max_cell_area_um2" in config["qc"]
        assert config["qc"]["max_cell_area_um2"] == 1500

    def test_spacia_mcmc_audit_correct(self, config):
        """Audit fix: MCMC must be 50000,20000,100,2 (not stale 5000)."""
        sv = config["analysis"]["f3_ccc"]["spacia_validation"]
        assert sv["mcmc_params"] == "50000,20000,100,2"

    def test_spacia_threshold_block_complete(self, config):
        sv = config["analysis"]["f3_ccc"]["spacia_validation"]
        for k in ("mcmc_params", "n_cells_subsample",
                  "dist_cutoff_um", "top_n_pairs", "multiple_testing"):
            assert k in sv, f"spacia_validation missing {k}"
        assert sv["multiple_testing"] == "bonferroni"   # audit correction
        # Hardcoded path must be gone
        assert "spacia_repo" not in sv

    def test_liana_section_externalized(self, config):
        li = config["analysis"]["f3_ccc"]["liana"]
        for k in ("method", "n_perms", "min_cells", "expr_prop",
                  "fdr_method", "fdr_threshold"):
            assert k in li, f"liana missing {k}"
        assert li["method"] == "rank_aggregate"     # audit correction
        assert li["fdr_method"] == "bh"             # audit correction

    def test_f3_load_spacia_params_from_config(self):
        F3 = importlib.import_module("analysis.F3_spacia_ccc_validation")
        params = F3.load_spacia_params(CONFIG_PATH)
        assert params["mcmc_params"] == "50000,20000,100,2"
        assert params["n_cells"] == 5000
        assert params["dist_cutoff"] == 30
        assert params["top_n"] == 30

    def test_f3_load_spacia_params_falls_back_to_defaults(self, tmp_path):
        F3 = importlib.import_module("analysis.F3_spacia_ccc_validation")
        # Non-existent path → defaults
        params = F3.load_spacia_params(tmp_path / "nope.yaml")
        assert params["mcmc_params"] == F3.MCMC_PARAMS
        assert params["n_cells"] == F3.N_CELLS


# ───────────────────────────────────────────────────────────────────────
# B4 — methods excluded due to gene-panel size, documented in manuscript
# ───────────────────────────────────────────────────────────────────────
class TestB4MethodExclusions:
    def test_manuscript_documents_panel_constraints(self):
        assert MANUSCRIPT_TEX.exists()
        tex = MANUSCRIPT_TEX.read_text()
        # Note S1.1 must enumerate the four excluded tools by name
        for tool in ("NicheDE", "hdWGCNA", "SCENIC", "CellTypist"):
            assert tool in tex, f"manuscript does not mention {tool}"
        assert "289" in tex  # panel size context
        # The audit-memo subsection header
        assert "S1.1" in tex or "Methods discarded" in tex


# ───────────────────────────────────────────────────────────────────────
# B5 — Spacia binary resolution
# ───────────────────────────────────────────────────────────────────────
class TestB5Spacia:
    def test_spacia_path_env_override(self, tmp_path, monkeypatch):
        fake = tmp_path / "spacia.py"
        fake.touch()
        monkeypatch.setenv("SPACIA_PATH", str(fake))
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        assert paths_mod.spacia_tool_path() == fake.resolve()

    def test_spacia_path_default_is_in_repo(self, monkeypatch):
        monkeypatch.delenv("SPACIA_PATH", raising=False)
        import utils.paths as paths_mod
        importlib.reload(paths_mod)
        p = paths_mod.spacia_tool_path()
        assert p.parts[-3:] == ("external", "Spacia", "spacia.py")

    @pytest.mark.parametrize("script", SHELL_SCRIPTS, ids=lambda p: p.name)
    def test_shell_scripts_pass_bash_syntax(self, script):
        ret = subprocess.run(["bash", "-n", str(script)],
                             capture_output=True, text=True)
        assert ret.returncode == 0, ret.stderr

    @pytest.mark.parametrize("script", SHELL_SCRIPTS, ids=lambda p: p.name)
    def test_shell_scripts_honor_env_resolution(self, script):
        """Shell scripts must derive COUNTS/META/OUTBASE from env vars (no /home/mdiaz)."""
        txt = script.read_text()
        assert "/home/mdiaz" not in txt
        assert "XENIUM_PROJECT_ROOT" in txt
        assert "SPACIA_PATH" in txt


# ───────────────────────────────────────────────────────────────────────
# Integration smoke test — CLI --help on the F0/F3 entry points
# ───────────────────────────────────────────────────────────────────────
class TestCLIEntryPoints:
    @pytest.mark.parametrize("rel", [
        "pipeline/scripts/analysis/F0_reannotation_v3.py",
        "pipeline/scripts/analysis/F3_spacia_ccc_validation.py",
    ])
    def test_script_help_exits_zero(self, rel):
        script = REPO_ROOT / rel
        ret = subprocess.run(
            [sys.executable, str(script), "--help"],
            capture_output=True, text=True, cwd=str(REPO_ROOT),
        )
        assert ret.returncode == 0, ret.stderr
        assert "--markers" in ret.stdout or "--config" in ret.stdout
