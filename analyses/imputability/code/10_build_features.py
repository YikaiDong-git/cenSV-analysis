#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
BEST6 = ROOT / "results" / "imputation" / "decompose" / "cosine_nj" / f"best6_{SUFFIX}"
OUT = ROOT / "analyses" / "imputability" / "output"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def load_best6():
    rows = []
    for chrom in CHROMS:
        f = BEST6 / chrom / f"{chrom}_imputation_summary.tsv"
        if not f.exists():
            print(f"[WARN] missing {f}", flush=True)
            continue
        df = pd.read_csv(f, sep="\t", low_memory=False)
        df["chrom"] = chrom
        rows.append(df[["chrom", "sv_index", "cv_r2"]].rename(
            columns={"cv_r2": "cv_r2_best6"}))
    return pd.concat(rows, ignore_index=True)


def main():
    v2_path = OUT / "per_sv_features_with_ld.tsv"
    if not v2_path.exists():
        print(f"[FAIL] {v2_path} missing", flush=True)
        sys.exit(1)
    v2 = pd.read_csv(v2_path, sep="\t", low_memory=False)
    print(f"[INFO] v2 base table: {len(v2)} SVs x {v2.shape[1]} cols", flush=True)

    best6 = load_best6()
    print(f"[INFO] best6 union: {len(best6)} SVs across {best6['chrom'].nunique()} chroms",
          flush=True)

    n_dup = best6.duplicated(subset=["chrom", "sv_index"]).sum()
    print(f"[INFO] best6 duplicate (chrom, sv_index) count: {n_dup}", flush=True)
    if n_dup > 0:
        print("[FAIL] sv_index not unique within chrom — stopping", flush=True)
        sys.exit(1)

    v2["cv_r2_best8"] = v2["cv_r2"]
    v2["sv_index"] = v2["sv_index"].astype(int)
    best6["sv_index"] = best6["sv_index"].astype(int)
    merged = v2.merge(best6, on=["chrom", "sv_index"], how="left")
    print(f"[INFO] merged: {len(merged)} SVs x {merged.shape[1]} cols", flush=True)

    cov = merged["cv_r2_best6"].notna().mean()
    n_paired = merged["cv_r2_best6"].notna().sum()
    print(f"[GATE] cv_r2_best6 coverage on best8 SV table: {cov*100:.2f}% "
          f"({n_paired}/{len(merged)})", flush=True)
    if n_paired >= 100:
        sub = merged.dropna(subset=["cv_r2_best6", "cv_r2_best8"])
        r = sub["cv_r2_best6"].corr(sub["cv_r2_best8"])
        print(f"[GATE] Pearson r(cv_r2_best6, cv_r2_best8) = {r:.4f}  "
              f"(threshold 0.95)", flush=True)
        if r < 0.95:
            print("[WARN] sanity gate failed: best6 vs best8 r below 0.95. "
                  "Continuing but flagging in log for manual review.", flush=True)
    print(f"[INFO] best6 mean cv_r2 = {merged['cv_r2_best6'].mean():.4f}, "
          f"best8 mean cv_r2 = {merged['cv_r2_best8'].mean():.4f}", flush=True)

    out_path = OUT / "per_sv_features_v3.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(merged)} x {merged.shape[1]})", flush=True)


if __name__ == "__main__":
    main()
