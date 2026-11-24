#!/usr/bin/env python3
import argparse
import csv
import json
import os
import socket
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_pheno_diag(pheno_dir: Path, pheno_name: str) -> pd.DataFrame:
    rows = []
    diag_files = sorted(pheno_dir.rglob("*__diag.tsv"))
    sys.stderr.write(
        f"[load] pheno={pheno_name} dir={pheno_dir}  diag_files={len(diag_files)}\n"
    )
    for p in diag_files:
        with open(p, newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            recs = list(rdr)
        if len(recs) != 1:
            raise ValueError(f"expected 1 row in {p}, got {len(recs)}")
        rows.append(recs[0])
    full = pd.DataFrame(rows)
    sys.stderr.write(f"[load] pheno={pheno_name} n_rows={len(full)}\n")
    return full


def load_manifest(manifest_path: Path) -> pd.DataFrame:
    df = pd.read_csv(manifest_path, sep="\t")
    expected_cols = {"pair_id", "sv_locus", "arm", "sv_chrom"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"manifest missing required cols: {missing}")
    sys.stderr.write(
        f"[manifest] {manifest_path}  n_rows={len(df)}  "
        f"cols={len(df.columns)}\n"
    )
    return df


def atomic_write_tsv(df: pd.DataFrame, out_path: Path) -> None:
    tmp = out_path.with_suffix(out_path.suffix + f".tmp.{os.getpid()}")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, out_path)


def atomic_write_json(obj: dict, out_path: Path) -> None:
    tmp = out_path.with_suffix(out_path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=str)
    os.replace(tmp, out_path)


def derive_qc_pass(df: pd.DataFrame, suffix: str) -> pd.Series:
    conv = df.get(f"converged_{suffix}", pd.Series(["FALSE"] * len(df)))
    warn = df.get(f"warn_unconverged_{suffix}", pd.Series(["TRUE"] * len(df)))
    cap = df.get(f"near_cap_flag_{suffix}", pd.Series(["TRUE"] * len(df)))
    return (
        (conv.astype(str).str.upper() == "TRUE") &
        (warn.astype(str).str.upper() != "TRUE") &
        (cap.astype(str).str.upper() != "TRUE")
    )


def compute_partition_counts(wide: pd.DataFrame,
                             pheno_a: str, pheno_b: str) -> pd.DataFrame:
    pip_a = pd.to_numeric(wide.get(f"sv_pip_{pheno_a}"), errors="coerce")
    pip_b = pd.to_numeric(wide.get(f"sv_pip_{pheno_b}"), errors="coerce")
    status_a = wide.get(f"status_{pheno_a}", pd.Series([""] * len(wide))).fillna("")
    status_b = wide.get(f"status_{pheno_b}", pd.Series([""] * len(wide))).fillna("")
    ok_a = ~pip_a.isna()
    ok_b = ~pip_b.isna()
    has_a = ok_a | (status_a != "")
    has_b = ok_b | (status_b != "")
    stub_a = has_a & ~ok_a
    stub_b = has_b & ~ok_b
    return pd.DataFrame([
        {"partition": "n_both_ok",            "count": int((ok_a & ok_b).sum())},
        {"partition": "n_a_ok_b_stub",        "count": int((ok_a & stub_b).sum())},
        {"partition": "n_a_stub_b_ok",        "count": int((stub_a & ok_b).sum())},
        {"partition": "n_both_stub",          "count": int((stub_a & stub_b).sum())},
        {"partition": "n_a_only_in_manifest", "count": int((has_a & ~has_b).sum())},
        {"partition": "n_b_only_in_manifest", "count": int((~has_a & has_b).sum())},
        {"partition": "n_total_rows",         "count": int(len(wide))},
    ])


def scatter_jpg(wide: pd.DataFrame, pheno_x: str, pheno_y: str,
                out_path: Path, partition: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 8.0), dpi=200)
    x = pd.to_numeric(wide[f"sv_pip_{pheno_x}"], errors="coerce")
    y = pd.to_numeric(wide[f"sv_pip_{pheno_y}"], errors="coerce")
    mask = (~x.isna()) & (~y.isna())
    n_paired = int(mask.sum())
    x_p = x[mask].values
    y_p = y[mask].values
    high = np.maximum(x_p, y_p) >= 0.5
    ax.scatter(x_p[~high], y_p[~high], s=6, alpha=0.25, c="#888888",
               edgecolors="none", label=f"low-PIP (max<0.5) n={(~high).sum()}")
    if high.any():
        ax.scatter(x_p[high], y_p[high], s=20, alpha=1.0, c="#ff7700",
                   edgecolors="#cc4400", linewidths=0.6,
                   label=f"high-PIP (max≥0.5) n={high.sum()}")
    ax.plot([0, 1], [0, 1], "--", c="#bbbbbb", lw=0.8)
    for thr in (0.5, 0.9):
        ax.axhline(thr, c="#cc4400", lw=0.5, ls=":")
        ax.axvline(thr, c="#cc4400", lw=0.5, ls=":")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(f"sv_pip  {pheno_x}")
    ax.set_ylabel(f"sv_pip  {pheno_y}")
    ax.set_title(f"SV-PIP cross-pheno  n_paired={n_paired}")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    pmap = {row["partition"]: int(row["count"]) for _, row in partition.iterrows()}
    txt = (
        f"both_ok={pmap.get('n_both_ok',0)}\n"
        f"{pheno_x}_only_ok={pmap.get('n_a_ok_b_stub',0)}\n"
        f"{pheno_y}_only_ok={pmap.get('n_a_stub_b_ok',0)}\n"
        f"both_stub={pmap.get('n_both_stub',0)}\n"
        f"total={pmap.get('n_total_rows',0)}"
    )
    ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, family="monospace",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="#cccccc"))
    fig.tight_layout()
    tmp = out_path.with_suffix(out_path.suffix + f".tmp.{os.getpid()}")
    fig.savefig(tmp, format="jpg", pil_kwargs={"quality": 95})
    plt.close(fig)
    os.replace(tmp, out_path)


def summary_table(wide: pd.DataFrame, pheno_a: str, pheno_b: str) -> pd.DataFrame:
    a = pd.to_numeric(wide[f"sv_pip_{pheno_a}"], errors="coerce")
    b = pd.to_numeric(wide[f"sv_pip_{pheno_b}"], errors="coerce")
    paired_mask = (~a.isna()) & (~b.isna())
    rows = []
    for thr in (0.1, 0.3, 0.5, 0.7, 0.9):
        n_a = int(((a >= thr) & paired_mask).sum())
        n_b = int(((b >= thr) & paired_mask).sum())
        n_both = int(((a >= thr) & (b >= thr) & paired_mask).sum())
        rows.append({
            "pip_thr": thr,
            "n_paired": int(paired_mask.sum()),
            f"n_{pheno_a}_ge_thr": n_a,
            f"n_{pheno_b}_ge_thr": n_b,
            "n_both_ge_thr": n_both,
        })
    return pd.DataFrame(rows)


def top_50_per_pheno(wide: pd.DataFrame, pheno: str) -> pd.DataFrame:
    col_pip = f"sv_pip_{pheno}"
    col_top = f"top_snp_pip_{pheno}"
    col_status = f"status_{pheno}"
    df = wide.dropna(subset=[col_pip]).copy()
    df[col_pip] = pd.to_numeric(df[col_pip], errors="coerce")
    df[col_top] = pd.to_numeric(df.get(col_top, np.nan), errors="coerce")
    df["qc_pass"] = derive_qc_pass(df, pheno)
    sort_cols = ["qc_pass", col_pip, col_top, "sv_locus"]
    ascending = [False, False, True, True]
    df_sorted = df.sort_values(sort_cols, ascending=ascending)
    return df_sorted.head(50)


def git_sha(repo_dir: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode == 0:
        return result.stdout.strip()
    return ""


def write_run_manifest(out_dir: Path, args: argparse.Namespace,
                       df_a: pd.DataFrame, df_b: pd.DataFrame,
                       wide: pd.DataFrame, partition: pd.DataFrame,
                       coverage_ok: bool, n_manifest: int) -> None:
    repo_dir = Path(CENSV_ROOT + "")
    sha = git_sha(repo_dir)
    distinct = {
        col: sorted(set(df_a.get(col, pd.Series([], dtype=str)).dropna().astype(str).tolist() +
                        df_b.get(col, pd.Series([], dtype=str)).dropna().astype(str).tolist()))
        for col in ("wrapper_id", "susieR_version", "sessionInfo_hash", "R_version")
        if (col in df_a.columns) or (col in df_b.columns)
    }
    manifest = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "host": socket.gethostname(),
        "argv": sys.argv,
        "git_sha": sha,
        "pheno_a": args.pheno_a,
        "pheno_b": args.pheno_b,
        "sweep_root": str(args.sweep_root),
        "manifest_arg": str(args.manifest) if args.manifest else None,
        "n_manifest": n_manifest,
        "n_diag_a": len(df_a),
        "n_diag_b": len(df_b),
        "n_wide_rows": len(wide),
        "partial_flag": bool(args.partial),
        "coverage_ok": coverage_ok,
        "partition_counts": partition.to_dict(orient="records"),
        "distinct_provenance": distinct,
    }
    atomic_write_json(manifest, out_dir / "RUN_MANIFEST.json")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep-root", type=Path, required=True,
                    help="sweep_6519_locked_v3/diag/")
    ap.add_argument("--pheno-a", default="I9_HYPTENSESS")
    ap.add_argument("--pheno-b", default="30210_irnt")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, default=None,
                    help="Optional 6,519-row MASTER_MANIFEST_6519.tsv for "
                         "coverage cross-check. If provided, refuses to "
                         "write top-50 unless coverage >= 95%% per pheno "
                         "(override with --partial).")
    ap.add_argument("--partial", action="store_true",
                    help="Bypass manifest coverage gate (writes top-50 even "
                         "when partial). Required when --manifest is given "
                         "and coverage < 95%%.")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df_a = load_pheno_diag(args.sweep_root / args.pheno_a, args.pheno_a)
    df_b = load_pheno_diag(args.sweep_root / args.pheno_b, args.pheno_b)

    cols_keep = [
        "sv_locus", "arm", "chrom",
        "sv_pip", "top_snp_pip", "sum_pip",
        "Z_SV", "N_GWAS", "n_LD_variants",
        "s_rss", "s_rss_tier", "kriging_outlier_count",
        "converged", "warn_unconverged", "near_cap_flag",
        "wrapper_id", "susieR_version", "status",
    ]
    df_a_k = df_a[[c for c in cols_keep if c in df_a.columns]].copy()
    df_b_k = df_b[[c for c in cols_keep if c in df_b.columns]].copy()
    df_a_k = df_a_k.rename(columns={
        c: f"{c}_{args.pheno_a}" for c in df_a_k.columns
        if c not in ("sv_locus", "arm", "chrom")
    })
    df_b_k = df_b_k.rename(columns={
        c: f"{c}_{args.pheno_b}" for c in df_b_k.columns
        if c not in ("sv_locus", "arm", "chrom")
    })
    wide = df_a_k.merge(df_b_k, on=["sv_locus", "arm", "chrom"], how="outer")
    sys.stderr.write(f"[wide] n_rows={len(wide)}  n_cols={wide.shape[1]}\n")

    coverage_ok = True
    n_manifest = 0
    if args.manifest is not None:
        m = load_manifest(args.manifest)
        n_manifest = len(m)
        cov_a = len(df_a) / n_manifest if n_manifest else 0
        cov_b = len(df_b) / n_manifest if n_manifest else 0
        sys.stderr.write(
            f"[coverage] {args.pheno_a}={len(df_a)}/{n_manifest} "
            f"({100*cov_a:.1f}%)  "
            f"{args.pheno_b}={len(df_b)}/{n_manifest} ({100*cov_b:.1f}%)\n"
        )
        if cov_a < 0.95 or cov_b < 0.95:
            coverage_ok = False
            if not args.partial:
                sys.stderr.write(
                    "[FATAL] coverage < 95%% per pheno. Pass --partial "
                    "to dry-run on the in-flight tree.\n"
                )
                return 3

    atomic_write_tsv(wide, args.out_dir / "wide_pip_matrix.tsv")

    summary = summary_table(wide, args.pheno_a, args.pheno_b)
    atomic_write_tsv(summary, args.out_dir / "cross_pheno_summary.tsv")
    sys.stderr.write(f"[summary]\n{summary.to_string(index=False)}\n")

    partition = compute_partition_counts(wide, args.pheno_a, args.pheno_b)
    atomic_write_tsv(partition, args.out_dir / "cross_pheno_partition.tsv")
    sys.stderr.write(f"[partition]\n{partition.to_string(index=False)}\n")

    scatter_jpg(wide, args.pheno_a, args.pheno_b,
                args.out_dir / "cross_pheno_scatter.jpg", partition)

    for pheno in (args.pheno_a, args.pheno_b):
        col = f"sv_pip_{pheno}"
        if col not in wide.columns:
            continue
        top = top_50_per_pheno(wide, pheno)
        atomic_write_tsv(top, args.out_dir / f"top_50_pip_{pheno}.tsv")

    write_run_manifest(args.out_dir, args, df_a, df_b, wide, partition,
                       coverage_ok, n_manifest)

    sys.stderr.write(
        f"[done] wrote wide_pip_matrix.tsv + cross_pheno_summary.tsv "
        f"+ cross_pheno_partition.tsv + cross_pheno_scatter.jpg + "
        f"top_50 per pheno + RUN_MANIFEST.json to {args.out_dir}\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
