#!/usr/bin/env python3
"""Phase B aggregate: compare pad vs drop sv_pip head-to-head."""
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
    p.add_argument("--root", type=Path,
                   default=BASE / "output/full_restart/zpad_headhead")
    p.add_argument("--out-dir", type=Path,
                   default=BASE / "output/full_restart/zpad_headhead/aggregate")
    return p.parse_args()


def atomic_write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.parent / f"{path.stem}.tmp.{os.getpid()}{path.suffix}"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def load_variant(root: Path, variant: str) -> pd.DataFrame:
    files = sorted((root / variant).glob("chr*/*__diag.tsv"))
    if not files:
        return pd.DataFrame()
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        df["variant"] = variant
        df["src_file"] = f.name
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    pad = load_variant(args.root, "pad")
    drop = load_variant(args.root, "drop")
    sys.stderr.write(f"[zpad-compare] pad rows: {len(pad)}; drop rows: {len(drop)}\n")
    if pad.empty or drop.empty:
        sys.stderr.write("[zpad-compare] missing variant data; can't compare\n")
        return 1

    key_cols = ["region_id", "pheno", "sv_locus", "sv_index"]
    keep = key_cols + ["sv_pip", "cv_r2", "converged", "n_iter", "p_snps", "G_svs"]
    pad_s = pad[keep].rename(columns={"sv_pip": "sv_pip_pad", "converged": "converged_pad",
                                       "n_iter": "n_iter_pad", "p_snps": "p_snps_pad"})
    drop_s = drop[keep].rename(columns={"sv_pip": "sv_pip_drop", "converged": "converged_drop",
                                         "n_iter": "n_iter_drop", "p_snps": "p_snps_drop"})

    cmp = pad_s.merge(
        drop_s[key_cols + ["sv_pip_drop", "converged_drop", "n_iter_drop", "p_snps_drop"]],
        on=key_cols, how="inner",
    )

    cmp["pip_diff"] = cmp["sv_pip_drop"] - cmp["sv_pip_pad"]
    cmp["pip_abs_diff"] = cmp["pip_diff"].abs()
    cmp = cmp.sort_values("pip_abs_diff", ascending=False).reset_index(drop=True)
    atomic_write_tsv(cmp, args.out_dir / "zpad_headhead_compare.tsv")

    def to_bool(s):
        s = str(s).strip().lower()
        return s in ("true", "t", "1")
    cmp["converged_pad_b"] = cmp["converged_pad"].apply(to_bool)
    cmp["converged_drop_b"] = cmp["converged_drop"].apply(to_bool)

    summary = pd.DataFrame([
        {"metric": "n_pairs", "value": len(cmp)},
        {"metric": "n_pip_diff_gt_0p05", "value": int((cmp["pip_abs_diff"] > 0.05).sum())},
        {"metric": "n_pip_diff_gt_0p1",  "value": int((cmp["pip_abs_diff"] > 0.10).sum())},
        {"metric": "n_pip_diff_gt_0p5",  "value": int((cmp["pip_abs_diff"] > 0.50).sum())},
        {"metric": "n_hits_pad_only_pip0p5",  "value": int(((cmp["sv_pip_pad"] >= 0.5) & (cmp["sv_pip_drop"] < 0.5)).sum())},
        {"metric": "n_hits_drop_only_pip0p5", "value": int(((cmp["sv_pip_drop"] >= 0.5) & (cmp["sv_pip_pad"] < 0.5)).sum())},
        {"metric": "n_hits_both_pip0p5",      "value": int(((cmp["sv_pip_pad"] >= 0.5) & (cmp["sv_pip_drop"] >= 0.5)).sum())},
        {"metric": "n_converged_pad",  "value": int(cmp["converged_pad_b"].sum())},
        {"metric": "n_converged_drop", "value": int(cmp["converged_drop_b"].sum())},
        {"metric": "spearman_rho_sv_pip", "value": float(cmp[["sv_pip_pad", "sv_pip_drop"]].corr("spearman").iloc[0, 1])},
    ])
    atomic_write_tsv(summary, args.out_dir / "zpad_headhead_summary.tsv")

    print("\n=== Phase B: z=0 PAD vs DROP head-to-head summary ===")
    print(summary.to_string(index=False))
    print(f"\nTop 15 PIP differences (drop - pad):")
    print(cmp[["region_id", "pheno", "sv_locus", "sv_pip_pad", "sv_pip_drop",
               "pip_diff", "p_snps_pad", "p_snps_drop", "converged_pad", "converged_drop"]]
          .head(15).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
