"""
Phase F — Spacia CCC Validation
Validates top LIANA+ L-R interactions using Spacia (Bayesian MIL).
Spacia confirms spatial proximity of signaling pairs, reducing FP ~50%.

Usage (xenium_pipeline env):
    python F3_spacia_ccc_validation.py

Steps:
    1. Export Spacia inputs from sdata.zarr
    2. Filter top LIANA+ interactions to panel genes
    3. Run Spacia for each unique (source, target, ligand, receptor)
    4. Aggregate results and generate validation figure
"""

import os
import sys
import subprocess
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root, results_root, sdata_path, spacia_tool_path

# Defaults derived from environment (overridable via CLI / config)
DEFAULT_CONFIG_PATH = project_root() / "pipeline" / "config" / "config_lung.yaml"
ZARR = sdata_path()
LIANA_TOP = results_root() / "02_biology/ccc_liana_L2_v3/liana_L2_top50.csv"
OUT_DIR = results_root() / "02_biology/phase_f_spacia"
FIG_DIR = results_root() / "figures/phase_f_spacia"
SPACIA_PY = spacia_tool_path()
SPACIA_ENV = os.environ.get("SPACIA_CONDA_ENV", "spacia")

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUT_DIR / "F3_spacia.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ── Default Spacia parameters (audit-corrected; see manuscript §S1.5) ────────
# These match what produced the canonical results (50,000 iter / 20,000 burn-in).
# Override via CLI flags or pipeline/config/config_lung.yaml::analysis.f3_ccc.spacia_validation
MCMC_PARAMS = "50000,20000,100,2"   # Ntotal,Nwarm,Nthin,Nchains
N_CELLS = 5000                       # -nc: cells per type to subsample
DIST_CUTOFF = 30                     # -d:  juxtacrine range in µm
TOP_N = 30                           # LIANA interactions to validate


def load_spacia_params(config_path: Path | None) -> dict:
    """Load Spacia validation params from config YAML if present; otherwise return module defaults."""
    out = {"mcmc_params": MCMC_PARAMS, "n_cells": N_CELLS,
           "dist_cutoff": DIST_CUTOFF, "top_n": TOP_N}
    if config_path is None or not config_path.exists():
        return out
    import yaml
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}
    sv = (cfg.get("analysis", {})
             .get("f3_ccc", {})
             .get("spacia_validation", {}))
    if "mcmc_params" in sv:
        out["mcmc_params"] = str(sv["mcmc_params"])
    if "n_cells_subsample" in sv:
        out["n_cells"] = int(sv["n_cells_subsample"])
    if "dist_cutoff_um" in sv:
        out["dist_cutoff"] = int(sv["dist_cutoff_um"])
    if "top_n_pairs" in sv:
        out["top_n"] = int(sv["top_n_pairs"])
    return out


def export_spacia_inputs(sdata_path, out_dir):
    """Load sdata.zarr and export counts + metadata for Spacia."""
    import spatialdata as sd

    log.info("Loading sdata.zarr …")
    sdata = sd.read_zarr(str(sdata_path))
    adata = list(sdata.tables.values())[0]

    log.info(f"  {adata.n_obs} cells × {adata.n_vars} genes")

    # Raw counts (Spacia normalizes internally)
    import scipy.sparse as sp
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    counts_df = pd.DataFrame(
        X.astype(np.float32),
        index=adata.obs_names,
        columns=adata.var_names,
    )

    counts_path = out_dir / "spacia_counts.txt"
    log.info(f"  Writing counts ({counts_df.shape}) …")
    counts_df.to_csv(counts_path, sep="\t")

    # Metadata: X, Y, cell_type
    coords = adata.obsm["spatial"]
    meta_df = pd.DataFrame(
        {
            "X": coords[:, 0],
            "Y": coords[:, 1],
            "cell_type": adata.obs["cell_type_L2"].values,
        },
        index=adata.obs_names,
    )
    meta_path = out_dir / "spacia_metadata.txt"
    log.info(f"  Writing metadata …")
    meta_df.to_csv(meta_path, sep="\t")

    panel_genes = set(adata.var_names)
    return counts_path, meta_path, panel_genes


def select_interactions(liana_path, panel_genes, top_n):
    """Pick top N unique (source, target, ligand, receptor) pairs with both genes in panel."""
    df = pd.read_csv(liana_path)

    # Filter: both ligand and receptor must be in panel (handle complexes with '|')
    def genes_in_panel(gene_str, panel):
        parts = str(gene_str).split("|")
        return all(g in panel for g in parts)

    mask = df.apply(
        lambda r: genes_in_panel(r["ligand_complex"], panel_genes)
        and genes_in_panel(r["receptor_complex"], panel_genes),
        axis=1,
    )
    df_valid = df[mask].copy()
    log.info(f"LIANA interactions with panel genes: {mask.sum()}/{len(df)}")

    # Deduplicate: unique (source, target, ligand, receptor)
    df_unique = df_valid.drop_duplicates(
        subset=["source", "target", "ligand_complex", "receptor_complex"]
    ).head(top_n)

    interactions_path = OUT_DIR / "F3_selected_interactions.csv"
    df_unique.to_csv(interactions_path, index=False)
    log.info(f"Selected {len(df_unique)} unique interactions for Spacia")
    return df_unique


def run_spacia_job(counts_path, meta_path, source, target, ligand, receptor, out_dir,
                   mcmc_params=MCMC_PARAMS, n_cells=N_CELLS, dist_cutoff=DIST_CUTOFF):
    """Run a single Spacia job; returns output directory path."""
    job_name = f"{source}__{target}__{ligand}__{receptor}".replace("/", "-")
    job_out = out_dir / "jobs" / job_name
    job_out.mkdir(parents=True, exist_ok=True)

    result_marker = job_out / "B_and_FDR.csv"  # success marker
    if result_marker.exists():
        log.info(f"  [cached] {job_name}")
        return str(job_out), job_name, True

    # Handle complex ligands/receptors: use first gene if '|'
    sf_gene = ligand.split("|")[0]
    rf_gene = receptor.split("|")[0]

    cmd = (
        f"conda run -n {SPACIA_ENV} python {SPACIA_PY} "
        f"{counts_path} {meta_path} "
        f"-sc {source} -rc {target} "
        f"-sf {sf_gene} -rf {rf_gene} "
        f"-d {dist_cutoff} "
        f"-m {mcmc_params} "
        f"-nc {n_cells} "
        f"-o {job_out}"
    )

    log.info(f"  Running: {source} → {target} | {ligand} → {receptor}")
    ret = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)

    if ret.returncode != 0:
        log.warning(f"  FAILED: {job_name}\n{ret.stderr[-500:]}")
        return str(job_out), job_name, False

    return str(job_out), job_name, True


def parse_spacia_result(job_out_dir, job_name, ligand, receptor):
    """Parse Spacia Pathway_betas.csv + B_and_FDR.csv.

    Validation criterion: Pathway_betas pval_adj < 0.05 (sender effect on receiver).
    """
    job_out = Path(job_out_dir)

    result = {
        "job_name": job_name,
        "ligand": ligand,
        "receptor": receptor,
        "beta": np.nan,
        "pathway_pval": np.nan,
        "pathway_pval_adj": np.nan,
        "fdr": np.nan,
        "b_coef": np.nan,
        "spatially_validated": False,
    }

    pb_file = job_out / "Pathway_betas.csv"
    bfdr_file = job_out / "B_and_FDR.csv"

    if pb_file.exists():
        try:
            pb = pd.read_csv(pb_file, index_col=0)
            result["beta"] = float(pb["Beta"].iloc[0])
            result["pathway_pval"] = float(pb["pval"].iloc[0])
            result["pathway_pval_adj"] = float(pb["pval_adj"].iloc[0])
            result["spatially_validated"] = result["pathway_pval_adj"] < 0.05
        except Exception as e:
            log.warning(f"Could not parse Pathway_betas for {job_name}: {e}")

    if bfdr_file.exists():
        try:
            bfdr = pd.read_csv(bfdr_file, index_col=0)
            result["fdr"] = float(bfdr["FDR"].iloc[0])
            result["b_coef"] = float(bfdr["b"].iloc[0])
        except Exception:
            pass

    return result


def make_validation_figure(merged_df, fig_dir):
    """Bubble plot: LIANA scaled_weight vs Spacia PIP, colored by validation status."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # ── Panel A: LIANA weight vs Spacia PIP ──────────────────────────────────
    ax = axes[0]
    colors = ["#d73027" if v else "#4575b4" for v in merged_df["spatially_validated"]]
    sizes = (merged_df["scaled_weight"] * 200).clip(20, 500)

    sc = ax.scatter(
        merged_df["scaled_weight"],
        merged_df["pip"].fillna(0),
        c=colors,
        s=sizes,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.5,
    )

    ax.axhline(0.5, color="gray", linestyle="--", lw=1.2, label="PIP = 0.5 threshold")
    ax.set_xlabel("LIANA+ scaled_weight", fontsize=12)
    ax.set_ylabel("Spacia PIP", fontsize=12)
    ax.set_title("LIANA+ × Spacia Cross-validation", fontsize=13, fontweight="bold")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d73027", label="Spatially validated (PIP ≥ 0.5)"),
        Patch(facecolor="#4575b4", label="Not spatially validated"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9)

    # Label top interactions
    top_rows = merged_df.nlargest(8, "pip")
    for _, row in top_rows.iterrows():
        ax.annotate(
            f"{row['ligand_complex']}→{row['receptor_complex']}\n({row['source']}→{row['target']})",
            (row["scaled_weight"], row["pip"] if not np.isnan(row["pip"]) else 0),
            fontsize=6,
            xytext=(5, 5),
            textcoords="offset points",
        )

    # ── Panel B: Bar chart of PIP by interaction ─────────────────────────────
    ax2 = axes[1]
    df_sorted = merged_df.sort_values("pip", ascending=True).tail(25)
    y_labels = [
        f"{r['ligand_complex']}→{r['receptor_complex']}\n{r['source']}→{r['target']}"
        for _, r in df_sorted.iterrows()
    ]
    bar_colors = ["#d73027" if v else "#4575b4" for v in df_sorted["spatially_validated"]]
    bars = ax2.barh(range(len(df_sorted)), df_sorted["pip"].fillna(0), color=bar_colors, alpha=0.8)
    ax2.set_yticks(range(len(df_sorted)))
    ax2.set_yticklabels(y_labels, fontsize=7)
    ax2.axvline(0.5, color="gray", linestyle="--", lw=1.2)
    ax2.set_xlabel("Spacia PIP", fontsize=12)
    ax2.set_title("Spacia PIP — Top Interactions", fontsize=13, fontweight="bold")
    ax2.set_xlim(0, 1.05)

    plt.tight_layout()
    out_path = fig_dir / "F3_spacia_validation.png"
    plt.savefig(str(out_path), dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"Saved validation figure → {out_path}")


def make_summary_table_figure(merged_df, fig_dir):
    """Summary stats figure."""
    n_total = len(merged_df)
    n_validated = merged_df["spatially_validated"].sum()
    n_not_validated = n_total - n_validated
    pct = 100 * n_validated / n_total if n_total > 0 else 0

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")

    summary = merged_df[["source", "target", "ligand_complex", "receptor_complex",
                          "scaled_weight", "pip", "spatially_validated"]].copy()
    summary["pip"] = summary["pip"].round(3)
    summary["scaled_weight"] = summary["scaled_weight"].round(3)
    summary = summary.sort_values("pip", ascending=False)

    table_data = [summary.columns.tolist()] + summary.values.tolist()
    table = ax.table(
        cellText=summary.values.tolist(),
        colLabels=summary.columns.tolist(),
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.3)

    ax.set_title(
        f"Spacia Validation Summary: {n_validated}/{n_total} validated ({pct:.1f}%)",
        fontsize=12,
        fontweight="bold",
        pad=20,
    )

    out_path = fig_dir / "F3_spacia_summary_table.png"
    plt.savefig(str(out_path), dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"Saved summary table → {out_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase F — Spacia CCC Validation")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH,
                        help=f"Pipeline config YAML (default: {DEFAULT_CONFIG_PATH})")
    parser.add_argument("--mcmc_params", type=str, default=None,
                        help="Override config: 'Ntotal,Nwarm,Nthin,Nchains'")
    parser.add_argument("--n_cells", type=int, default=None, help="Override config: -nc subsample")
    parser.add_argument("--dist_cutoff", type=int, default=None, help="Override config: -d µm")
    parser.add_argument("--top_n", type=int, default=None, help="Override config: top-N LIANA pairs")
    args = parser.parse_args()

    params = load_spacia_params(args.config)
    if args.mcmc_params: params["mcmc_params"] = args.mcmc_params
    if args.n_cells:     params["n_cells"]     = args.n_cells
    if args.dist_cutoff: params["dist_cutoff"] = args.dist_cutoff
    if args.top_n:       params["top_n"]       = args.top_n

    log.info("=" * 60)
    log.info("Phase F — Spacia CCC Validation")
    log.info("=" * 60)
    log.info(f"  config:       {args.config}")
    log.info(f"  mcmc_params:  {params['mcmc_params']}")
    log.info(f"  n_cells (-nc):{params['n_cells']}")
    log.info(f"  dist_cutoff:  {params['dist_cutoff']} µm")
    log.info(f"  top_n_pairs:  {params['top_n']}")
    log.info(f"  spacia_py:    {SPACIA_PY}")

    # Step 1: Export inputs
    counts_path, meta_path, panel_genes = export_spacia_inputs(ZARR, OUT_DIR)

    # Step 2: Select interactions
    interactions = select_interactions(LIANA_TOP, panel_genes, params["top_n"])

    # Step 3: Run Spacia jobs
    results = []
    for i, (_, row) in enumerate(interactions.iterrows()):
        log.info(f"[{i+1}/{len(interactions)}] {row['source']} → {row['target']} | "
                 f"{row['ligand_complex']} → {row['receptor_complex']}")

        job_out_dir, job_name, success = run_spacia_job(
            counts_path, meta_path,
            row["source"], row["target"],
            row["ligand_complex"], row["receptor_complex"],
            OUT_DIR,
            mcmc_params=params["mcmc_params"],
            n_cells=params["n_cells"],
            dist_cutoff=params["dist_cutoff"],
        )

        if success:
            parsed = parse_spacia_result(job_out_dir, job_name,
                                          row["ligand_complex"], row["receptor_complex"])
        else:
            parsed = {
                "job_name": job_name,
                "ligand": row["ligand_complex"],
                "receptor": row["receptor_complex"],
                "pip": np.nan,
                "beta": np.nan,
                "spatially_validated": False,
            }

        parsed.update({
            "source": row["source"],
            "target": row["target"],
            "scaled_weight": row["scaled_weight"],
            "specificity_rank": row.get("specificity_rank", np.nan),
        })
        results.append(parsed)

    # Step 4: Aggregate and save
    results_df = pd.DataFrame(results)
    results_path = OUT_DIR / "F3_spacia_results.csv"
    results_df.to_csv(results_path, index=False)
    log.info(f"Saved results → {results_path}")

    # Merge with LIANA
    liana_df = pd.read_csv(LIANA_TOP)
    merged = interactions.merge(
        results_df[["source", "target", "ligand", "receptor", "pip", "beta", "spatially_validated"]],
        left_on=["source", "target", "ligand_complex", "receptor_complex"],
        right_on=["source", "target", "ligand", "receptor"],
        how="left",
    )
    merged_path = OUT_DIR / "F3_liana_spacia_merged.csv"
    merged.to_csv(merged_path, index=False)

    # Summary statistics
    n_total = len(results_df)
    n_validated = results_df["spatially_validated"].sum()
    n_run = results_df["pip"].notna().sum()
    log.info("=" * 60)
    log.info(f"RESULTS: {n_validated}/{n_run} interactions spatially validated (PIP ≥ 0.5)")
    log.info(f"         {n_total - n_run} jobs failed or incomplete")
    log.info("=" * 60)

    # Step 5: Figures
    if n_run > 0:
        make_validation_figure(merged, FIG_DIR)
        make_summary_table_figure(merged, FIG_DIR)

    log.info("Phase F — COMPLETE")
    return results_df


if __name__ == "__main__":
    main()
