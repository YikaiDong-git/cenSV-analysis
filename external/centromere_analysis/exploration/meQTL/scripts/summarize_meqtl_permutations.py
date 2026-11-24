#!/usr/bin/env python3
"""Compare an observed retained meQTL scan with genotype-label permutations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--permuted", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def scan_metrics(path: Path, label: str) -> dict[str, object]:
    frame = pd.read_csv(path, sep="\t")
    threshold = frame[frame["retained_reason"].astype(str).str.contains("threshold")].copy()
    sv_cols = [column for column in ["sv_source", "sv_chrom", "sv_index"] if column in frame]
    per_sv_min = frame.groupby(sv_cols, dropna=False)["pvalue"].min() if sv_cols else pd.Series(dtype=float)
    return {
        "scan": label,
        "n_retained": len(frame),
        "n_p_le_1e_4": len(threshold),
        "n_q_lt_0_01": int((frame.get("qvalue_retained_bh", pd.Series(dtype=float)) < 0.01).sum()),
        "min_pvalue": float(frame["pvalue"].min()) if not frame.empty else np.nan,
        "median_per_sv_min_p": float(per_sv_min.median()) if not per_sv_min.empty else np.nan,
        "n_svs_with_min_p_le_1e_6": int((per_sv_min <= 1e-6).sum()),
    }


def main() -> int:
    args = parse_args()
    rows = [scan_metrics(args.observed, "observed")]
    rows.extend(scan_metrics(path, f"permutation_{i + 1}") for i, path in enumerate(args.permuted))
    table = pd.DataFrame(rows)
    perm = table[table["scan"].str.startswith("permutation_")]
    observed = table.iloc[0]
    mean_perm_p_tail = float(perm["n_p_le_1e_4"].mean())
    mean_perm_q_hits = float(perm["n_q_lt_0_01"].mean())
    summary = {
        "n_permutations": len(perm),
        "observed_n_p_le_1e_4": int(observed["n_p_le_1e_4"]),
        "permutation_mean_n_p_le_1e_4": mean_perm_p_tail,
        "retained_count_enrichment": float(observed["n_p_le_1e_4"] / max(mean_perm_p_tail, 1.0)),
        "empirical_p_for_n_p_le_1e_4": float((1 + (perm["n_p_le_1e_4"] >= observed["n_p_le_1e_4"]).sum()) / (len(perm) + 1)),
        "observed_n_q_lt_0_01": int(observed["n_q_lt_0_01"]),
        "permutation_mean_n_q_lt_0_01": mean_perm_q_hits,
        "permutation_min_n_q_lt_0_01": int(perm["n_q_lt_0_01"].min()),
        "permutation_max_n_q_lt_0_01": int(perm["n_q_lt_0_01"].max()),
        "q_lt_0_01_count_enrichment": float(
            observed["n_q_lt_0_01"] / max(mean_perm_q_hits, 1.0)
        ),
        "empirical_p_for_n_q_lt_0_01": float(
            (1 + (perm["n_q_lt_0_01"] >= observed["n_q_lt_0_01"]).sum())
            / (len(perm) + 1)
        ),
        "note": (
            f"{len(perm)} permutations calibrate gross inflation only; increase the seed set "
            "for tail probabilities."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
