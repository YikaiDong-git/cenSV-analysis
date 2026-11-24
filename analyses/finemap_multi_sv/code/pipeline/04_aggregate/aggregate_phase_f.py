#!/usr/bin/env python3
"""Phase F aggregate + 4-way compare (pad, drop, ctwas-100, ctwas-500)."""
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
    p.add_argument("--f-dir", type=Path, default=BASE / "output/full_restart/finemap_ctwas_iter500")
    p.add_argument("--pad-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/full_restart_aggregate.tsv")
    p.add_argument("--drop-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/drop_aggregate.tsv")
    p.add_argument("--ctwas-aggregate", type=Path,
                   default=BASE / "output/full_restart/aggregate/ctwas_aggregate.tsv")
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

    diag_files = sorted(args.f_dir.glob("chr*/*__diag.tsv"))
    sys.stderr.write(f"[aggregate-F] found {len(diag_files)} F diag.tsv\n")
    if not diag_files:
        return 0

    frames = []
    for f in diag_files:
        df = pd.read_csv(f, sep="\t")
        df["src_file"] = f.name
        frames.append(df)
    f_df = pd.concat(frames, ignore_index=True)
    sys.stderr.write(f"[aggregate-F] total F rows: {len(f_df)}\n")
    atomic_write_tsv(f_df, args.out_dir / "ctwas_iter500_aggregate.tsv")

    fhits01 = f_df[f_df["sv_pip"] >= 0.1].copy().sort_values("sv_pip", ascending=False)
    fhits05 = f_df[f_df["sv_pip"] >= 0.5].copy().sort_values("sv_pip", ascending=False)
    atomic_write_tsv(fhits01, args.out_dir / "ctwas_iter500_hits_pip0p1.tsv")
    atomic_write_tsv(fhits05, args.out_dir / "ctwas_iter500_hits_pip0p5.tsv")

    if not all(p.exists() for p in (args.pad_aggregate, args.drop_aggregate, args.ctwas_aggregate)):
        sys.stderr.write("[aggregate-F] some aggregate missing\n")
        return 0
    pad = pd.read_csv(args.pad_aggregate, sep="\t")
    drop = pd.read_csv(args.drop_aggregate, sep="\t")
    cw = pd.read_csv(args.ctwas_aggregate, sep="\t")
    key_cols = ["region_id", "pheno", "sv_locus", "sv_index"]

    pad_s = pad[key_cols + ["sv_pip", "cv_r2", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_pad", "converged": "conv_pad", "n_iter": "niter_pad"})
    drop_s = drop[key_cols + ["sv_pip", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_drop", "converged": "conv_drop", "n_iter": "niter_drop"})
    cw_s = cw[key_cols + ["sv_pip", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_ctwas100", "converged": "conv_ctwas100", "n_iter": "niter_ctwas100"})
    f_s = f_df[key_cols + ["sv_pip", "converged", "n_iter"]].rename(
        columns={"sv_pip": "pip_ctwas500", "converged": "conv_ctwas500", "n_iter": "niter_ctwas500"})

    cmp = (pad_s.merge(drop_s, on=key_cols)
                .merge(cw_s, on=key_cols)
                .merge(f_s, on=key_cols))
    atomic_write_tsv(cmp, args.out_dir / "four_way_compare.tsv")

    def to_bool(s):
        s = str(s).strip().lower()
        return s in ("true", "t", "1")
    for col in ("conv_pad", "conv_drop", "conv_ctwas100", "conv_ctwas500"):
        cmp[f"{col}_b"] = cmp[col].apply(to_bool)

    summary = pd.DataFrame([
        {"metric": "n_pairs_compared", "value": len(cmp)},
        {"metric": "n_pad_hits_pip0p1", "value": int((cmp["pip_pad"] >= 0.1).sum())},
        {"metric": "n_drop_hits_pip0p1", "value": int((cmp["pip_drop"] >= 0.1).sum())},
        {"metric": "n_ctwas100_hits_pip0p1", "value": int((cmp["pip_ctwas100"] >= 0.1).sum())},
        {"metric": "n_ctwas500_hits_pip0p1", "value": int((cmp["pip_ctwas500"] >= 0.1).sum())},
        {"metric": "n_pad_hits_pip0p5", "value": int((cmp["pip_pad"] >= 0.5).sum())},
        {"metric": "n_drop_hits_pip0p5", "value": int((cmp["pip_drop"] >= 0.5).sum())},
        {"metric": "n_ctwas100_hits_pip0p5", "value": int((cmp["pip_ctwas100"] >= 0.5).sum())},
        {"metric": "n_ctwas500_hits_pip0p5", "value": int((cmp["pip_ctwas500"] >= 0.5).sum())},
        {"metric": "n_all4_hits_pip0p5", "value": int(((cmp["pip_pad"] >= 0.5) & (cmp["pip_drop"] >= 0.5) & (cmp["pip_ctwas100"] >= 0.5) & (cmp["pip_ctwas500"] >= 0.5)).sum())},
        {"metric": "spearman_pad_drop", "value": float(cmp[["pip_pad", "pip_drop"]].corr("spearman").iloc[0,1])},
        {"metric": "spearman_pad_ctwas500", "value": float(cmp[["pip_pad", "pip_ctwas500"]].corr("spearman").iloc[0,1])},
        {"metric": "spearman_drop_ctwas500", "value": float(cmp[["pip_drop", "pip_ctwas500"]].corr("spearman").iloc[0,1])},
        {"metric": "spearman_ctwas100_ctwas500", "value": float(cmp[["pip_ctwas100", "pip_ctwas500"]].corr("spearman").iloc[0,1])},
        {"metric": "n_conv_pad", "value": int(cmp["conv_pad_b"].sum())},
        {"metric": "n_conv_drop", "value": int(cmp["conv_drop_b"].sum())},
        {"metric": "n_conv_ctwas100", "value": int(cmp["conv_ctwas100_b"].sum())},
        {"metric": "n_conv_ctwas500", "value": int(cmp["conv_ctwas500_b"].sum())},
    ])
    atomic_write_tsv(summary, args.out_dir / "four_way_summary.tsv")

    print("\n=== Phase F + 4-way compare summary ===")
    print(summary.to_string(index=False))

    print(f"\nTop 25 cTWAS-iter500 hits (sv_pip ≥ 0.1):")
    print(cmp[cmp["pip_ctwas500"] >= 0.1]
          .sort_values("pip_ctwas500", ascending=False)
          [["region_id", "pheno", "sv_locus", "pip_pad", "pip_drop", "pip_ctwas100", "pip_ctwas500",
            "conv_ctwas100", "conv_ctwas500"]]
          .head(25).to_string(index=False))

    print(f"\nHits common to all 4 variants (≥ 0.5):")
    common = cmp[(cmp["pip_pad"] >= 0.5) & (cmp["pip_drop"] >= 0.5) & (cmp["pip_ctwas100"] >= 0.5) & (cmp["pip_ctwas500"] >= 0.5)]
    if len(common) > 0:
        print(common[["region_id", "pheno", "sv_locus", "pip_pad", "pip_drop", "pip_ctwas100", "pip_ctwas500"]].to_string(index=False))
    else:
        print("(none)")

    print(f"\nCase 1 chr1_sv616 × ER-neg BC head-to-head:")
    c1 = cmp[(cmp["sv_locus"].str.contains("chr1_sv616", na=False)) & (cmp["pheno"] == "C3_BREAST_ERNEG_EXALLC")]
    print(c1[["sv_locus", "pheno", "pip_pad", "pip_drop", "pip_ctwas100", "pip_ctwas500",
              "conv_pad", "conv_drop", "conv_ctwas100", "conv_ctwas500"]].to_string(index=False))

    print(f"\nchr11_sv253 × 5255_irnt (the only Phase E all-3 survivor) under 4-way:")
    s = cmp[(cmp["sv_locus"] == "chr11_sv253_INS_ori_best8_midpoint500") & (cmp["pheno"] == "5255_irnt")]
    print(s[["sv_locus", "pheno", "pip_pad", "pip_drop", "pip_ctwas100", "pip_ctwas500",
             "conv_pad", "conv_drop", "conv_ctwas100", "conv_ctwas500"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
