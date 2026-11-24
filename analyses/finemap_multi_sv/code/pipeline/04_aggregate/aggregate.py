#!/usr/bin/env python3
"""Aggregate diag.tsv from output/full_restart/finemap/ into one master TSV."""
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
    p.add_argument("--finemap-dir", type=Path, default=BASE / "output/full_restart/finemap")
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

    diag_files = sorted(args.finemap_dir.glob("chr*/*__diag.tsv"))
    sys.stderr.write(f"[aggregate-v2] found {len(diag_files)} diag.tsv files\n")
    if not diag_files:
        sys.stderr.write("[aggregate-v2] no diag.tsv yet; nothing to aggregate\n")
        return 0

    frames = []
    for f in diag_files:
        df = pd.read_csv(f, sep="\t")
        df["src_dir"] = f.parent.name
        df["src_file"] = f.name
        frames.append(df)
    agg = pd.concat(frames, ignore_index=True)
    sys.stderr.write(f"[aggregate-v2] total rows: {len(agg)}\n")

    atomic_write_tsv(agg, args.out_dir / "full_restart_aggregate.tsv")

    hits01 = agg[agg["sv_pip"] >= 0.1].copy().sort_values("sv_pip", ascending=False)
    hits05 = agg[agg["sv_pip"] >= 0.5].copy().sort_values("sv_pip", ascending=False)
    atomic_write_tsv(hits01, args.out_dir / "full_restart_hits_pip0p1.tsv")
    atomic_write_tsv(hits05, args.out_dir / "full_restart_hits_pip0p5.tsv")

    summary_rows = []
    for (chrom, conv), g in agg.groupby([agg["region_id"].str.extract(r"^(chr\d+)")[0], "converged"]):
        summary_rows.append({
            "chrom": chrom,
            "converged": conv,
            "n_rows": len(g),
            "max_sv_pip": float(g["sv_pip"].max()),
            "n_hits_pip0p1": int((g["sv_pip"] >= 0.1).sum()),
            "n_hits_pip0p5": int((g["sv_pip"] >= 0.5).sum()),
        })
    summary = pd.DataFrame(summary_rows)
    atomic_write_tsv(summary, args.out_dir / "full_restart_summary.tsv")

    print(f"\n=== V2 FULL-RESTART summary ===")
    print(f"Tasks: {agg[['region_id','pheno']].drop_duplicates().shape[0]}")
    print(f"SV rows: {len(agg)}")
    print(f"Hits at sv_pip ≥ 0.1: {len(hits01)}")
    print(f"Hits at sv_pip ≥ 0.5: {len(hits05)}")
    print(f"Converged TRUE: {int((agg['converged']=='TRUE').sum())} / {len(agg)}  ({100*(agg['converged']=='TRUE').mean():.1f}%)")
    print(f"\nTop 20 hits (sv_pip):")
    cols_print = ["sv_locus", "pheno", "sv_pip", "cv_r2", "converged", "n_iter", "region_id"]
    cols_print = [c for c in cols_print if c in hits01.columns]
    print(hits01[cols_print].head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
