#!/usr/bin/env python3
"""
F1 — Spatial Domain Consensus (Novae × BANKSY).

Compares domain assignments from two independent methods (Novae foundation model
and BANKSY) using:
  - Adjusted Rand Index (ARI) — agreement between partitions, chance-corrected
  - Adjusted Mutual Information (AMI) — information-theoretic alternative
  - Hungarian-matched relabeling — find optimal 1-to-1 domain mapping

Output: consensus assignment per cell + cross-tabulation table.

Citation context:
  - ARI: Hubert & Arabie (1985) Comparing partitions. J Classif 2, 193-218.
  - AMI: Vinh et al. (2010) Information theoretic measures for clustering comparison. JMLR 11, 2837-2854.
  - Salas et al. Nat Methods 22, 813-823 (2025) — recommends consensus across multiple domain methods.
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_rand_score,
    adjusted_mutual_info_score,
    normalized_mutual_info_score,
    confusion_matrix,
)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--novae", required=True, type=Path,
                   help="Path to novae_domains.csv from F1_novae_domains.py")
    p.add_argument("--banksy", required=True, type=Path,
                   help="Path to banksy_domains.csv from F1_banksy_domains.R")
    p.add_argument("--novae-col", default="novae_domain",
                   help="Column in novae csv with primary domain assignment")
    p.add_argument("--banksy-col", default="banksy_domain")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def hungarian_relabel(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, dict]:
    """Relabel `b` to maximize agreement with `a` using Hungarian algorithm on
    confusion matrix.
    Returns: (relabeled_b, mapping_dict)
    """
    from scipy.optimize import linear_sum_assignment
    a_uniq = np.array(sorted(set(a.tolist())))
    b_uniq = np.array(sorted(set(b.tolist())))
    a_to_idx = {lbl: i for i, lbl in enumerate(a_uniq.tolist())}
    b_to_idx = {lbl: i for i, lbl in enumerate(b_uniq.tolist())}
    a_idx = np.array([a_to_idx[x] for x in a])
    b_idx = np.array([b_to_idx[x] for x in b])
    n_a, n_b = len(a_uniq), len(b_uniq)
    # Build confusion matrix manually to avoid sklearn auto-handling of non-aligned labels
    cm = np.zeros((n_a, n_b), dtype=np.int64)
    np.add.at(cm, (a_idx, b_idx), 1)
    # Pad to square
    n = max(n_a, n_b)
    padded = np.zeros((n, n), dtype=cm.dtype)
    padded[:n_a, :n_b] = cm
    row_ind, col_ind = linear_sum_assignment(-padded)
    # row_ind[k]=i means assignment row->col, both span [0,n)
    mapping = {}
    for r, c in zip(row_ind.tolist(), col_ind.tolist()):
        if c >= n_b:
            continue                                # padded col, no real banksy class
        b_label = b_uniq[c]
        if r < n_a:
            mapping[str(b_label)] = str(a_uniq[r])
        else:
            mapping[str(b_label)] = f"unmatched_b{c}"
    relabeled = np.array([mapping.get(str(x), str(x)) for x in b], dtype=object)
    return relabeled, mapping


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log.info("=" * 70)
    log.info("F1 — Domain Consensus (Novae × BANKSY)")
    log.info("=" * 70)

    log.info(f"Loading Novae domains: {args.novae}")
    novae_df = pd.read_csv(args.novae, index_col=0)
    log.info(f"  {novae_df.shape}, columns: {list(novae_df.columns)[:5]}")

    log.info(f"Loading BANKSY domains: {args.banksy}")
    banksy_df = pd.read_csv(args.banksy)
    if "cell_id" in banksy_df.columns:
        banksy_df = banksy_df.set_index("cell_id")
    log.info(f"  {banksy_df.shape}, columns: {list(banksy_df.columns)}")

    # Align indices
    common = novae_df.index.intersection(banksy_df.index)
    log.info(f"  Common cells: {len(common):,} ({len(common)/max(len(novae_df), len(banksy_df)):.1%})")
    if len(common) < 100:
        raise RuntimeError("Too few common cells between the two methods")

    novae_labels = novae_df.loc[common, args.novae_col].astype(str).values
    banksy_labels = banksy_df.loc[common, args.banksy_col].astype(str).values

    # Metrics
    ari = adjusted_rand_score(novae_labels, banksy_labels)
    ami = adjusted_mutual_info_score(novae_labels, banksy_labels)
    nmi = normalized_mutual_info_score(novae_labels, banksy_labels)
    log.info(f"\n  Adjusted Rand Index (ARI)        : {ari:.4f}")
    log.info(f"  Adjusted Mutual Information (AMI): {ami:.4f}")
    log.info(f"  Normalized MI (NMI)              : {nmi:.4f}")

    # Hungarian relabeling
    relabeled_banksy, mapping = hungarian_relabel(novae_labels, banksy_labels)
    agreement = (novae_labels == relabeled_banksy).mean()
    log.info(f"  Post-Hungarian agreement         : {agreement:.4f}")

    # Confusion table
    ct = pd.crosstab(
        pd.Series(novae_labels, name="novae_domain", index=common),
        pd.Series(banksy_labels, name="banksy_domain", index=common),
    )
    ct_path = args.outdir / "novae_vs_banksy_crosstab.csv"
    ct.to_csv(ct_path)
    log.info(f"  Saved: {ct_path} (shape {ct.shape})")

    # Consensus column: agree-only labels (Novae domain when Hungarian-matched)
    consensus = pd.DataFrame(index=common)
    consensus["novae_domain"] = novae_labels
    consensus["banksy_domain"] = banksy_labels
    consensus["banksy_relabeled_to_novae"] = relabeled_banksy
    consensus["agreement"] = (consensus["novae_domain"] == consensus["banksy_relabeled_to_novae"])
    consensus["consensus_domain"] = np.where(
        consensus["agreement"],
        consensus["novae_domain"],
        "Disagreement",
    )
    cons_path = args.outdir / "consensus_domains.csv"
    consensus.to_csv(cons_path)
    log.info(f"  Saved: {cons_path}")

    # Distribution
    log.info(f"\n  Consensus distribution:")
    log.info(consensus["consensus_domain"].value_counts().head(15).to_string())

    # Metrics JSON
    metrics = {
        "n_common_cells": int(len(common)),
        "n_novae_domains": int(len(np.unique(novae_labels))),
        "n_banksy_domains": int(len(np.unique(banksy_labels))),
        "ari": round(ari, 4),
        "ami": round(ami, 4),
        "nmi": round(nmi, 4),
        "hungarian_agreement": round(float(agreement), 4),
        "n_cells_consensus": int((consensus["consensus_domain"] != "Disagreement").sum()),
        "execution_seconds": round(time.time() - t0, 1),
    }
    (args.outdir / "consensus_metrics.json").write_text(json.dumps(metrics, indent=2))
    log.info(f"\n  Metrics: {metrics}")

    log.info(f"\n✓ F1 consensus complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
