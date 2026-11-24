#!/usr/bin/env python3
"""Aggregate Phase E (cTWAS-faithful) diag.tsv + 3-way compare vs pad + drop."""
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
    p.add_argument("--ctwas-dir", type=Path, default=BASE / "output/full_restart/finemap_ctwas")
    p.add_argument("--pad-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/full_restart_aggregate.tsv")
    p.add_argument("--drop-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/drop_aggregate.tsv")
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

    diag_files = sorted(args.ctwas_dir.glob("chr*/*__diag.tsv"))
    sys.stderr.write(f"[aggregate-ctwas] found {len(diag_files)} ctwas diag.tsv\n")
    if not diag_files:
        sys.stderr.write("[aggregate-ctwas] no ctwas diag.tsv yet\n")
        return 0

    frames = []
    for f in diag_files:
        df = pd.read_csv(f, sep="\t")
        df["src_file"] = f.name
        frames.append(df)
    ctwas = pd.concat(frames, ignore_index=True)
    sys.stderr.write(f"[aggregate-ctwas] total ctwas rows: {len(ctwas)}\n")
    atomic_write_tsv(ctwas, args.out_dir / "ctwas_aggregate.tsv")

    cw_hits01 = ctwas[ctwas["sv_pip"] >= 0.1].copy().sort_values("sv_pip", ascending=False)
    cw_hits05 = ctwas[ctwas["sv_pip"] >= 0.5].copy().sort_values("sv_pip", ascending=False)
    atomic_write_tsv(cw_hits01, args.out_dir / "ctwas_hits_pip0p1.tsv")
    atomic_write_tsv(cw_hits05, args.out_dir / "ctwas_hits_pip0p5.tsv")

    if not args.pad_aggregate.exists() or not args.drop_aggregate.exists():
        sys.stderr.write("[aggregate-ctwas] pad or drop aggregate missing; skipping 3-way\n")
        return 0
    pad = pd.read_csv(args.pad_aggregate, sep="\t")
    drop = pd.read_csv(args.drop_aggregate, sep="\t")
    key_cols = ["region_id", "pheno", "sv_locus", "sv_index"]

    pad_s = pad[key_cols + ["sv_pip", "cv_r2", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_pad", "converged": "conv_pad", "n_iter": "niter_pad"}
    )
    drop_s = drop[key_cols + ["sv_pip", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_drop", "converged": "conv_drop", "n_iter": "niter_drop"}
    )
    ctwas_s = ctwas[key_cols + ["sv_pip", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_ctwas", "converged": "conv_ctwas", "n_iter": "niter_ctwas"}
    )

    cmp = (pad_s
           .merge(drop_s, on=key_cols, how="inner")
           .merge(ctwas_s, on=key_cols, how="inner"))
    cmp["abs_max_diff"] = (
        (cmp["pip_pad"] - cmp["pip_ctwas"]).abs().combine(
            (cmp["pip_drop"] - cmp["pip_ctwas"]).abs(), max
        )
    )
    cmp = cmp.sort_values("abs_max_diff", ascending=False).reset_index(drop=True)
    atomic_write_tsv(cmp, args.out_dir / "three_way_compare.tsv")

    def to_bool(s):
        s = str(s).strip().lower()
        return s in ("true", "t", "1")
    for col in ("conv_pad", "conv_drop", "conv_ctwas"):
        cmp[f"{col}_b"] = cmp[col].apply(to_bool)

    summary = pd.DataFrame([
        {"metric": "n_pairs_compared", "value": len(cmp)},
        {"metric": "n_pad_hits_pip0p5", "value": int((cmp["pip_pad"] >= 0.5).sum())},
        {"metric": "n_drop_hits_pip0p5", "value": int((cmp["pip_drop"] >= 0.5).sum())},
        {"metric": "n_ctwas_hits_pip0p5", "value": int((cmp["pip_ctwas"] >= 0.5).sum())},
        {"metric": "n_all3_hits_pip0p5", "value": int(((cmp["pip_pad"] >= 0.5) & (cmp["pip_drop"] >= 0.5) & (cmp["pip_ctwas"] >= 0.5)).sum())},
        {"metric": "n_pad_hits_pip0p1", "value": int((cmp["pip_pad"] >= 0.1).sum())},
        {"metric": "n_drop_hits_pip0p1", "value": int((cmp["pip_drop"] >= 0.1).sum())},
        {"metric": "n_ctwas_hits_pip0p1", "value": int((cmp["pip_ctwas"] >= 0.1).sum())},
        {"metric": "spearman_pad_drop", "value": float(cmp[["pip_pad", "pip_drop"]].corr("spearman").iloc[0, 1])},
        {"metric": "spearman_pad_ctwas", "value": float(cmp[["pip_pad", "pip_ctwas"]].corr("spearman").iloc[0, 1])},
        {"metric": "spearman_drop_ctwas", "value": float(cmp[["pip_drop", "pip_ctwas"]].corr("spearman").iloc[0, 1])},
        {"metric": "n_conv_pad", "value": int(cmp["conv_pad_b"].sum())},
        {"metric": "n_conv_drop", "value": int(cmp["conv_drop_b"].sum())},
        {"metric": "n_conv_ctwas", "value": int(cmp["conv_ctwas_b"].sum())},
    ])
    atomic_write_tsv(summary, args.out_dir / "three_way_summary.tsv")

    print("\n=== Phase E: cTWAS-faithful summary + 3-way compare ===")
    print(summary.to_string(index=False))

    print(f"\nTop 20 cTWAS-faithful hits (sv_pip ≥ 0.1):")
    print(cmp[cmp["pip_ctwas"] >= 0.1]
          .sort_values("pip_ctwas", ascending=False)
          [["region_id", "pheno", "sv_locus", "pip_pad", "pip_drop", "pip_ctwas",
            "conv_pad", "conv_drop", "conv_ctwas"]]
          .head(20).to_string(index=False))

    print(f"\nHits common to all 3 variants (≥ 0.5):")
    common = cmp[(cmp["pip_pad"] >= 0.5) & (cmp["pip_drop"] >= 0.5) & (cmp["pip_ctwas"] >= 0.5)]
    print(common[["region_id", "pheno", "sv_locus", "pip_pad", "pip_drop", "pip_ctwas"]]
          .to_string(index=False))

    print(f"\nCase 1 (chr1_sv616) head-to-head:")
    c1 = cmp[cmp["sv_locus"].str.contains("chr1_sv616", na=False)]
    print(c1[["region_id", "pheno", "sv_locus", "pip_pad", "pip_drop", "pip_ctwas",
              "conv_pad", "conv_drop", "conv_ctwas"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
