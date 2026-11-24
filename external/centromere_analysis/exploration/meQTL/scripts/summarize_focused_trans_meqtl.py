#!/usr/bin/env python3
"""Merge focused trans-meQTL chunks and apply permutation max-statistic calibration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from summarize_meqtl_pilot import add_retained_bh, read_json_maybe_empty, read_tsv_maybe_empty


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-manifest", type=Path, required=True)
    parser.add_argument("--jobs-dir", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def empirical_p(null: np.ndarray, observed: np.ndarray) -> np.ndarray:
    return (1.0 + np.sum(null[:, None] >= observed[None, :], axis=0)) / (len(null) + 1.0)


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.chunk_manifest, sep="\t")
    candidates = pd.read_csv(args.candidates, sep="\t")
    observed_frames: list[pd.DataFrame] = []
    null_frames: list[pd.DataFrame] = []
    summaries: list[dict] = []
    for row in manifest.itertuples(index=False):
        job_dir = args.jobs_dir / row.job_id
        observed = read_tsv_maybe_empty(job_dir / "signals.tsv.gz")
        if not observed.empty:
            observed.insert(0, "job_id", row.job_id)
            observed.insert(1, "chunk_id", row.chunk_id)
            observed_frames.append(observed)
        null = read_tsv_maybe_empty(job_dir / "null_stats.tsv.gz")
        if not null.empty:
            null.insert(0, "job_id", row.job_id)
            null_frames.append(null)
        summary = read_json_maybe_empty(job_dir / "summary.json")
        if summary:
            summary.update({"job_id": row.job_id, "chunk_id": row.chunk_id})
            summaries.append(summary)

    observed = pd.concat(observed_frames, ignore_index=True) if observed_frames else pd.DataFrame()
    null_chunks = pd.concat(null_frames, ignore_index=True) if null_frames else pd.DataFrame()
    summary_df = pd.DataFrame(summaries)
    total_tests = int(summary_df["n_tests"].sum())
    observed = add_retained_bh(observed, total_tests)

    null_by_sv = (
        null_chunks.groupby(["permutation_seed", "sv_index", "sv_id"], as_index=False)
        .agg(max_abs_t=("max_abs_t", "max"), min_pvalue=("min_pvalue", "min"), n_p_le_threshold=("n_p_le_threshold", "sum"))
    )
    null_global = null_by_sv.groupby("permutation_seed", as_index=False).agg(
        max_abs_t=("max_abs_t", "max"), min_pvalue=("min_pvalue", "min")
    )

    if not observed.empty:
        observed_abs_t = observed["t_stat"].abs().to_numpy(dtype=float)
        global_null = null_global["max_abs_t"].to_numpy(dtype=float)
        observed["empirical_p_global_fwer"] = empirical_p(global_null, observed_abs_t)
        observed["empirical_p_sv_fwer"] = np.nan
        for sv_id, indices in observed.groupby("sv_id").groups.items():
            sv_null = null_by_sv.loc[null_by_sv["sv_id"] == sv_id, "max_abs_t"].to_numpy(dtype=float)
            observed.loc[indices, "empirical_p_sv_fwer"] = empirical_p(
                sv_null, observed.loc[indices, "t_stat"].abs().to_numpy(dtype=float)
            )
        observed = observed.merge(
            candidates[
                [
                    "sv_index",
                    "physical_locus_id",
                    "conditioning_consensus",
                    "n_robust_panels",
                    "n_observed_carriers",
                    "collapsed_sv_ids",
                    "n_collapsed_records",
                ]
            ],
            on="sv_index",
            how="left",
        )
        observed = observed.sort_values(["empirical_p_global_fwer", "pvalue"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    observed.to_csv(args.output_dir / "retained_trans_assoc.tsv.gz", sep="\t", index=False, compression="gzip")
    observed.head(1000).to_csv(args.output_dir / "top_trans_hits.tsv", sep="\t", index=False)
    null_by_sv.to_csv(args.output_dir / "permutation_maxima_by_sv.tsv.gz", sep="\t", index=False, compression="gzip")
    null_global.to_csv(args.output_dir / "permutation_global_maxima.tsv", sep="\t", index=False)
    summary_df.to_csv(args.output_dir / "job_summary.tsv", sep="\t", index=False)
    significant = observed[observed["empirical_p_global_fwer"] < 0.05] if not observed.empty else observed
    significant.to_csv(args.output_dir / "trans_hits_empirical_fwer_005.tsv", sep="\t", index=False)
    result = {
        "n_jobs": int(len(manifest)),
        "n_completed_jobs": int(len(summary_df)),
        "n_candidate_dosage_profiles": int(len(candidates)),
        "n_total_tests": total_tests,
        "n_permutations": int(null_global["permutation_seed"].nunique()),
        "n_retained": int(len(observed)),
        "n_bh_q_lt_0_01": int((observed["qvalue_retained_bh"] < 0.01).sum()) if not observed.empty else 0,
        "n_sv_fwer_lt_0_05": int((observed["empirical_p_sv_fwer"] < 0.05).sum()) if not observed.empty else 0,
        "n_global_fwer_lt_0_05": int(len(significant)),
        "minimum_empirical_p": 1.0 / (int(null_global["permutation_seed"].nunique()) + 1.0),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
