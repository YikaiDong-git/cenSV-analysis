#!/usr/bin/env python3
"""Merge and summarize meQTL pilot association job outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-manifest", type=Path, required=True)
    parser.add_argument("--jobs-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_tsv_maybe_empty(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except EmptyDataError:
        return pd.DataFrame()


def read_json_maybe_empty(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    return json.loads(path.read_text())


def add_retained_bh(df: pd.DataFrame, total_tests: int) -> pd.DataFrame:
    if df.empty or total_tests <= 0:
        df["qvalue_retained_bh"] = np.nan
        return df
    df = df.sort_values("pvalue").reset_index(drop=True)
    threshold_mask = df["retained_reason"].astype(str).str.contains("threshold")
    qvals = np.full(len(df), np.nan)
    threshold_idx = np.where(threshold_mask.to_numpy())[0]
    if len(threshold_idx):
        p = df.loc[threshold_idx, "pvalue"].to_numpy(dtype=float)
        ranks = np.arange(1, len(p) + 1, dtype=float)
        q = p * float(total_tests) / ranks
        q = np.minimum.accumulate(q[::-1])[::-1]
        qvals[threshold_idx] = np.minimum(q, 1.0)
    df["qvalue_retained_bh"] = qvals
    return df


def plot_retained_qq(df: pd.DataFrame, output: Path) -> None:
    threshold = df[df["retained_reason"].astype(str).str.contains("threshold")]
    fig, ax = plt.subplots(figsize=(5, 5))
    if threshold.empty:
        ax.text(0.5, 0.5, "No threshold-retained p-values", ha="center", va="center")
        ax.set_axis_off()
    else:
        p = np.sort(threshold["pvalue"].clip(lower=np.nextafter(0, 1)).to_numpy(dtype=float))
        observed = -np.log10(p)
        expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        ax.scatter(expected, observed, s=8, alpha=0.6)
        lim = max(float(expected.max()), float(observed.max()))
        ax.plot([0, lim], [0, lim], color="black", linewidth=1)
        ax.set_xlabel("Expected -log10(p), retained set")
        ax.set_ylabel("Observed -log10(p)")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def plot_manhattan(df: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4))
    if df.empty:
        ax.text(0.5, 0.5, "No retained associations", ha="center", va="center")
        ax.set_axis_off()
    else:
        plot_df = df.nsmallest(min(len(df), 200000), "pvalue").copy()
        plot_df["neglog10p"] = -np.log10(plot_df["pvalue"].clip(lower=np.nextafter(0, 1)))
        ax.scatter(plot_df["methyl_start"], plot_df["neglog10p"], s=5, alpha=0.5)
        ax.set_xlabel("Methylation position")
        ax.set_ylabel("-log10(p)")
        ax.set_title("Retained SV-methylation associations")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.job_manifest, sep="\t")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    summaries = []
    for row in manifest.itertuples(index=False):
        job_dir = args.jobs_dir / row.job_id
        df = read_tsv_maybe_empty(job_dir / "signals.tsv.gz")
        if not df.empty:
            df.insert(0, "job_id", row.job_id)
            df.insert(1, "chunk_id", row.chunk_id)
            df.insert(2, "sv_batch_id", row.sv_batch_id)
            if hasattr(row, "sv_source"):
                df.insert(3, "sv_source", row.sv_source)
            frames.append(df)
        summary = read_json_maybe_empty(job_dir / "summary.json")
        if summary:
            summary.update({"job_id": row.job_id, "chunk_id": row.chunk_id, "sv_batch_id": row.sv_batch_id})
            if hasattr(row, "sv_source"):
                summary["sv_source"] = row.sv_source
            summaries.append(summary)

    assoc = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    summary_df = pd.DataFrame(summaries)
    total_tests = int(summary_df["n_tests"].sum()) if "n_tests" in summary_df else 0
    assoc = add_retained_bh(assoc, total_tests)

    assoc_path = args.output_dir / "retained_assoc.tsv.gz"
    assoc.to_csv(assoc_path, sep="\t", index=False, compression="gzip")
    top_hits = assoc.sort_values("pvalue").head(1000) if not assoc.empty else pd.DataFrame()
    top_hits.to_csv(args.output_dir / "top_hits.tsv", sep="\t", index=False)
    summary_df.to_csv(args.output_dir / "job_summary.tsv", sep="\t", index=False)

    qvalue_summary = {
        "n_jobs": int(len(manifest)),
        "n_completed_summaries": int(len(summary_df)),
        "n_total_tests": total_tests,
        "n_retained_rows": int(len(assoc)),
        "n_threshold_rows": int(assoc["retained_reason"].astype(str).str.contains("threshold").sum()) if not assoc.empty else 0,
        "min_pvalue": float(assoc["pvalue"].min()) if not assoc.empty else None,
        "n_retained_bh_q_lt_0_1": int((assoc["qvalue_retained_bh"] < 0.1).sum()) if not assoc.empty else 0,
    }
    (args.output_dir / "qvalue_summary.json").write_text(json.dumps(qvalue_summary, indent=2) + "\n")
    pd.DataFrame([qvalue_summary]).to_csv(args.output_dir / "qvalue_summary.tsv", sep="\t", index=False)

    plot_retained_qq(assoc, args.output_dir / "qq_retained_pvalues.png")
    plot_manhattan(assoc, args.output_dir / "methylation_manhattan_retained.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
