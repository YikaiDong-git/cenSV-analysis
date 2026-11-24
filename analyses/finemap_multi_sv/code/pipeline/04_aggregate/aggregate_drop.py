#!/usr/bin/env python3
"""Aggregate diag.tsv from output/full_restart/finemap_drop/ + compare vs Phase A (pad)."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]


BASE = Path(
    CENSV_ROOT + "/"
    "analyses/finemap_multi_sv"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--drop-dir", type=Path, default=BASE / "output/full_restart/finemap_drop")
    p.add_argument("--pad-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/full_restart_aggregate.tsv")
    p.add_argument("--out-dir", type=Path, default=BASE / "output/full_restart/aggregate")
    return p.parse_args()


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}{path.suffix}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    diag_files = sorted(args.drop_dir.glob("chr*/*__diag.tsv"))
    sys.stderr.write(f"[aggregate-drop] found {len(diag_files)} drop diag.tsv files\n")
    if not diag_files:
        sys.stderr.write("[aggregate-drop] no drop diag.tsv yet\n")
        return 0

    frames = []
    for f in diag_files:
        df = pd.read_csv(f, sep="\t")
        df["src_file"] = f.name
        frames.append(df)
    drop = pd.concat(frames, ignore_index=True)
    sys.stderr.write(f"[aggregate-drop] total drop rows: {len(drop)}\n")
    atomic_write_tsv(drop, args.out_dir / "drop_aggregate.tsv")

    drop_hits01 = drop[drop["sv_pip"] >= 0.1].copy().sort_values("sv_pip", ascending=False)
    drop_hits05 = drop[drop["sv_pip"] >= 0.5].copy().sort_values("sv_pip", ascending=False)
    atomic_write_tsv(drop_hits01, args.out_dir / "drop_hits_pip0p1.tsv")
    atomic_write_tsv(drop_hits05, args.out_dir / "drop_hits_pip0p5.tsv")

    if not args.pad_aggregate.exists():
        sys.stderr.write(f"[aggregate-drop] pad aggregate missing: {args.pad_aggregate}\n")
        return 0
    pad = pd.read_csv(args.pad_aggregate, sep="\t")
    sys.stderr.write(f"[aggregate-drop] pad rows: {len(pad)}\n")

    key_cols = ["region_id", "pheno", "sv_locus", "sv_index"]
    pad_s = pad[key_cols + ["sv_pip", "cv_r2", "converged", "n_iter"]].rename(
        columns={"sv_pip": "sv_pip_pad", "converged": "conv_pad", "n_iter": "niter_pad"}
    )
    drop_s = drop[key_cols + ["sv_pip", "cv_r2", "converged", "n_iter"]].rename(
        columns={"sv_pip": "sv_pip_drop", "converged": "conv_drop", "n_iter": "niter_drop"}
    )
    cmp = pad_s.merge(drop_s[key_cols + ["sv_pip_drop", "conv_drop", "niter_drop"]],
                      on=key_cols, how="inner")
    cmp["pip_diff"] = cmp["sv_pip_drop"] - cmp["sv_pip_pad"]
    cmp["pip_abs_diff"] = cmp["pip_diff"].abs()
    cmp = cmp.sort_values("pip_abs_diff", ascending=False).reset_index(drop=True)
    atomic_write_tsv(cmp, args.out_dir / "pad_vs_drop_full_compare.tsv")

    def to_bool(s):
        s = str(s).strip().lower()
        return s in ("true", "t", "1")
    cmp["conv_pad_b"] = cmp["conv_pad"].apply(to_bool)
    cmp["conv_drop_b"] = cmp["conv_drop"].apply(to_bool)

    summary = pd.DataFrame([
        {"metric": "n_pairs_compared", "value": len(cmp)},
        {"metric": "n_pad_hits_pip0p1", "value": int((cmp["sv_pip_pad"] >= 0.1).sum())},
        {"metric": "n_drop_hits_pip0p1", "value": int((cmp["sv_pip_drop"] >= 0.1).sum())},
        {"metric": "n_pad_hits_pip0p5", "value": int((cmp["sv_pip_pad"] >= 0.5).sum())},
        {"metric": "n_drop_hits_pip0p5", "value": int((cmp["sv_pip_drop"] >= 0.5).sum())},
        {"metric": "n_drop_survivors_pip0p5", "value": int(((cmp["sv_pip_pad"] >= 0.5) & (cmp["sv_pip_drop"] >= 0.5)).sum())},
        {"metric": "n_pad_only_pip0p5", "value": int(((cmp["sv_pip_pad"] >= 0.5) & (cmp["sv_pip_drop"] < 0.5)).sum())},
        {"metric": "n_drop_only_pip0p5", "value": int(((cmp["sv_pip_pad"] < 0.5) & (cmp["sv_pip_drop"] >= 0.5)).sum())},
        {"metric": "n_conv_pad", "value": int(cmp["conv_pad_b"].sum())},
        {"metric": "n_conv_drop", "value": int(cmp["conv_drop_b"].sum())},
        {"metric": "spearman_rho_sv_pip", "value": float(cmp[["sv_pip_pad", "sv_pip_drop"]].corr("spearman").iloc[0, 1])},
    ])
    atomic_write_tsv(summary, args.out_dir / "pad_vs_drop_summary.tsv")

    print("\n=== Phase D: DROP-variant full-restart sweep summary ===")
    print(summary.to_string(index=False))
    print(f"\nTop 25 PIP differences (drop − pad), sorted by |Δ|:")
    print(cmp[["region_id", "pheno", "sv_locus", "sv_pip_pad", "sv_pip_drop",
               "pip_diff", "conv_pad", "conv_drop"]]
          .head(25).to_string(index=False))
    print(f"\nTop 25 DROP HITS at sv_pip ≥ 0.1:")
    print(cmp[cmp["sv_pip_drop"] >= 0.1]
          .sort_values("sv_pip_drop", ascending=False)
          [["region_id", "pheno", "sv_locus", "sv_pip_pad", "sv_pip_drop", "conv_pad", "conv_drop"]]
          .head(25).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
