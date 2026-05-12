"""
Phase F — Aggregate Spacia results + generate validation figures.
Run after F3b_run_all_spacia_jobs.sh completes.
Uses xenium_pipeline env (for matplotlib/pandas).
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import logging
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
SPACIA_DIR = BASE / "human_lung_cancer/results/02_biology/phase_f_spacia"
JOBS_DIR = SPACIA_DIR / "jobs"
INTERACTIONS_CSV = SPACIA_DIR / "F3_selected_interactions.csv"
FIG_DIR = BASE / "human_lung_cancer/results/figures/phase_f_spacia"
FIG_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

PVAL_THRESHOLD = 0.05  # Per-test pval_adj threshold (legacy — superseded by Bonferroni below)
FDR_THRESHOLD = 0.2    # B_and_FDR FDR threshold (secondary, spatial bag clustering)

# AUDIT FIX (2026-05-06): apply Bonferroni correction across ALL submitted tests.
# Spacia's internal pval_adj corrects only WITHIN a single pathway test, not
# across the N interactions submitted. The manuscript's 22/25 (88%) figure used
# uncorrected per-test p-values; the audited number with cross-test correction
# is reported below as `n_validated_bonferroni`.
APPLY_BONFERRONI_ACROSS_TESTS = True


def parse_spacia_result(job_dir):
    """Parse Pathway_betas.csv and B_and_FDR.csv from a Spacia output directory.

    Returns dict with: beta, pathway_pval, pathway_pval_adj, fdr, b_coef, has_result.
    Validation criterion: pathway_pval_adj < PVAL_THRESHOLD (sender effect on receiver).
    """
    job_dir = Path(job_dir)
    result = dict(
        beta=np.nan, pathway_pval=np.nan, pathway_pval_adj=np.nan,
        fdr=np.nan, b_coef=np.nan, has_result=False
    )

    pb_file = job_dir / "Pathway_betas.csv"
    bfdr_file = job_dir / "B_and_FDR.csv"

    if pb_file.exists():
        try:
            pb = pd.read_csv(pb_file, index_col=0)
            result["beta"] = float(pb["Beta"].iloc[0])
            result["pathway_pval"] = float(pb["pval"].iloc[0])
            result["pathway_pval_adj"] = float(pb["pval_adj"].iloc[0])
            result["has_result"] = True
        except Exception as e:
            pass

    if bfdr_file.exists():
        try:
            bfdr = pd.read_csv(bfdr_file, index_col=0)
            result["fdr"] = float(bfdr["FDR"].iloc[0])
            result["b_coef"] = float(bfdr["b"].iloc[0])
        except Exception:
            pass

    return result


def aggregate_results():
    interactions = pd.read_csv(INTERACTIONS_CSV)
    records = []

    for _, row in interactions.iterrows():
        job_name = f"{row['source']}__{row['target']}__{row['ligand_complex']}__{row['receptor_complex']}"
        job_dir = JOBS_DIR / job_name

        res = parse_spacia_result(job_dir)
        spatially_validated = (
            res["has_result"]
            and not np.isnan(res["pathway_pval_adj"])
            and res["pathway_pval_adj"] < PVAL_THRESHOLD
        )

        records.append({
            "source": row["source"],
            "target": row["target"],
            "ligand_complex": row["ligand_complex"],
            "receptor_complex": row["receptor_complex"],
            "scaled_weight": row["scaled_weight"],
            "beta": res["beta"],
            "pathway_pval": res["pathway_pval"],
            "pathway_pval_adj": res["pathway_pval_adj"],
            "fdr": res["fdr"],
            "b_coef": res["b_coef"],
            "has_result": res["has_result"],
            "spatially_validated": spatially_validated,
        })

    df = pd.DataFrame(records)

    # AUDIT FIX: Bonferroni across all submitted tests (the honest denominator
    # is the number of interactions we attempted, not the number that succeeded).
    if APPLY_BONFERRONI_ACROSS_TESTS:
        n_submitted = len(df)
        bonf_thresh = PVAL_THRESHOLD / n_submitted
        df["bonferroni_threshold"] = bonf_thresh
        df["validated_bonferroni"] = (
            df["has_result"]
            & df["pathway_pval"].notna()
            & (df["pathway_pval"] < bonf_thresh)
        )
        log.info(f"Bonferroni threshold across {n_submitted} submitted tests: p < {bonf_thresh:.2e}")

    results_path = SPACIA_DIR / "F3_spacia_results.csv"
    df.to_csv(results_path, index=False)

    n_run = df["has_result"].sum()
    n_val = df["spatially_validated"].sum()
    n_val_bonf = df["validated_bonferroni"].sum() if APPLY_BONFERRONI_ACROSS_TESTS else None
    log.info(f"Results: {n_run}/{len(df)} jobs completed")
    log.info(f"  Per-test pval_adj < {PVAL_THRESHOLD} (legacy):  {n_val}")
    if n_val_bonf is not None:
        log.info(f"  Bonferroni < 0.05 across {len(df)} attempted: {n_val_bonf}")
    return df


def figure_bubble_validation(df, fig_dir):
    """LIANA scaled_weight vs Spacia PIP — main validation figure."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor="white")

    # ── Panel A: scatter ──────────────────────────────────────────────────────
    ax = axes[0]
    df_plot = df[df["has_result"]].copy()
    # Use -log10(pval_adj) as y-axis
    df_plot["neg_log_pval"] = -np.log10(df_plot["pathway_pval_adj"].clip(1e-300))
    colors = ["#d73027" if v else "#4575b4" for v in df_plot["spatially_validated"]]
    sizes = (df_plot["scaled_weight"] * 200).clip(30, 500)

    ax.scatter(df_plot["scaled_weight"], df_plot["neg_log_pval"],
               c=colors, s=sizes, alpha=0.8, edgecolors="white", linewidths=0.6)
    ax.axhline(-np.log10(PVAL_THRESHOLD), color="#888", linestyle="--", lw=1.2, alpha=0.8,
               label=f"pval_adj = {PVAL_THRESHOLD}")
    ax.set_xlabel("LIANA+ scaled_weight", fontsize=12)
    ax.set_ylabel("-log10(Spacia pval_adj)", fontsize=12)
    ax.set_title("LIANA+ × Spacia Cross-Validation", fontsize=13, fontweight="bold")
    ax.set_ylim(-0.05, 1.1)

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#d73027", label=f"Spatially validated (pval_adj < {PVAL_THRESHOLD})"),
        Patch(color="#4575b4", label="Not validated"),
    ], fontsize=9, loc="upper left")

    # Annotate top 6 by significance
    for _, r in df_plot.nlargest(6, "neg_log_pval").iterrows():
        ax.annotate(
            f"{r['ligand_complex']}→{r['receptor_complex']}\n{r['source']}→{r['target']}",
            (r["scaled_weight"], r["neg_log_pval"]), fontsize=6,
            xytext=(5, 5), textcoords="offset points",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="#999"),
        )

    # ── Panel B: horizontal bar by PIP ───────────────────────────────────────
    ax2 = axes[1]
    df_bar = df_plot.sort_values("neg_log_pval", ascending=True)
    labels = [f"{r['ligand_complex']}→{r['receptor_complex']}  "
              f"({r['source']}→{r['target']})" for _, r in df_bar.iterrows()]
    bar_colors = ["#d73027" if v else "#4575b4" for v in df_bar["spatially_validated"]]
    ax2.barh(range(len(df_bar)), df_bar["neg_log_pval"], color=bar_colors, alpha=0.8, height=0.7)
    ax2.set_yticks(range(len(df_bar)))
    ax2.set_yticklabels(labels, fontsize=7)
    ax2.axvline(-np.log10(PVAL_THRESHOLD), color="#888", linestyle="--", lw=1.2)
    ax2.set_xlabel("-log10(Spacia pval_adj)", fontsize=12)
    ax2.set_title("Top Interactions — Spatial Validation", fontsize=13, fontweight="bold")

    plt.tight_layout()
    out = fig_dir / "F3_spacia_validation_bubble.png"
    plt.savefig(str(out), dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"Saved → {out}")


def figure_validation_summary(df, fig_dir):
    """Summary statistics and key findings figure."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")

    # ── Panel A: pie chart of validation status ───────────────────────────────
    ax = axes[0]
    n_val = df["spatially_validated"].sum()
    n_not = (df["has_result"] & ~df["spatially_validated"]).sum()
    n_no_result = (~df["has_result"]).sum()

    sizes = [n_val, n_not, n_no_result]
    labels = [f"Validated\n(n={n_val})", f"Not validated\n(n={n_not})", f"No result\n(n={n_no_result})"]
    colors_pie = ["#d73027", "#4575b4", "#aaaaaa"]
    wedge_props = dict(edgecolor="white", linewidth=2)
    ax.pie([s for s in sizes if s > 0],
           labels=[l for s, l in zip(sizes, labels) if s > 0],
           colors=[c for s, c in zip(sizes, colors_pie) if s > 0],
           autopct="%1.0f%%", startangle=90,
           wedgeprops=wedge_props, textprops=dict(fontsize=11))
    ax.set_title("Spacia Validation Status\n(Top 30 LIANA+ interactions)", fontsize=12, fontweight="bold")

    # ── Panel B: PIP distribution ─────────────────────────────────────────────
    ax2 = axes[1]
    df_has = df[df["has_result"]]
    ax2.hist(df_has["pathway_pval_adj"], bins=20, color="#4575b4", alpha=0.8, edgecolor="white")
    ax2.axvline(PVAL_THRESHOLD, color="#d73027", linestyle="--", lw=2, label=f"pval_adj = {PVAL_THRESHOLD}")
    ax2.set_xlabel("Spacia pval_adj", fontsize=12)
    ax2.set_ylabel("Number of interactions", fontsize=12)
    ax2.set_title("Distribution of Spacia PIP Values", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=10)

    plt.tight_layout()
    out = fig_dir / "F3_spacia_summary.png"
    plt.savefig(str(out), dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"Saved → {out}")


def main():
    log.info("=" * 60)
    log.info("Phase F — Aggregating Spacia results")
    log.info("=" * 60)

    df = aggregate_results()
    figure_bubble_validation(df, FIG_DIR)
    figure_validation_summary(df, FIG_DIR)

    # Print validated interactions
    val = df[df["spatially_validated"]].sort_values("pathway_pval_adj")
    log.info(f"\nSPATIALLY VALIDATED INTERACTIONS (pval_adj < {PVAL_THRESHOLD}):")
    for _, r in val.iterrows():
        log.info(f"  {r['source']} → {r['target']} | {r['ligand_complex']} → {r['receptor_complex']} "
                 f"| beta={r['beta']:.3f} | pval_adj={r['pathway_pval_adj']:.2e} | LIANA_w={r['scaled_weight']:.3f}")

    log.info("Phase F aggregation COMPLETE")
    return df


if __name__ == "__main__":
    main()
