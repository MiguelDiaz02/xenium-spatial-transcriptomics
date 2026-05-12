"""Pytest fixtures for the Xenium pipeline test suite.

Adds `pipeline/scripts/` to sys.path so tests can `import utils.paths`,
`import analysis.F3_spacia_ccc_validation`, etc.
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "pipeline" / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
