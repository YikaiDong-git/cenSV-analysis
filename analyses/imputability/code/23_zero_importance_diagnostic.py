#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

ZERO_FEATS = [
    "log_n_window_snps", "pericen_len_500000",
    "svtype_DUP", "svtype_INV",
    "best_snp_dist_mean",
    "sat_active",
]


def main():
    df = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                      sep="\t", low_memory=False)
    print(f"[INFO] {len(df)} SVs", flush=True)

    rows = []
    for f in ZERO_FEATS:
        if f not in df.columns:
            print(f"[SKIP] {f} not in columns", flush=True)
            continue
        col = df[f]
        n_nan = col.isna().sum()
        n_unique = col.nunique()
        cv = col.std() / col.mean() if col.mean() != 0 else np.nan
        rows.append({
            "feature": f,
            "n_nan": int(n_nan),
            "n_unique": int(n_unique),
            "min": float(col.min()) if n_nan < len(col) else np.nan,
            "max": float(col.max()) if n_nan < len(col) else np.nan,
            "mean": float(col.mean()) if n_nan < len(col) else np.nan,
            "std": float(col.std()) if n_nan < len(col) else np.nan,
            "value_counts_top3": col.value_counts().head(3).to_dict(),
        })
        print(f"\n=== {f} ===", flush=True)
        print(f"  n_nan: {n_nan}/{len(col)} ({n_nan/len(col)*100:.1f}%)", flush=True)
        print(f"  n_unique values: {n_unique}", flush=True)
        if n_nan < len(col):
            print(f"  range: [{col.min():.3f}, {col.max():.3f}]  "
                  f"mean={col.mean():.3f} std={col.std():.3f}", flush=True)
            print(f"  top-3 value counts:", flush=True)
            print(f"  {col.value_counts().head(5).to_string()}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "v3_zero_imp_diagnostic.tsv", sep="\t", index=False)
    print(f"\n[OK] wrote v3_zero_imp_diagnostic.tsv", flush=True)

    print("\n=== chrom-level features per-chrom ===", flush=True)
    for f in ["pericen_len_500000", "log_active_HOR_span", "n_HOR_intervals",
              "n_train_haps", "acrocentric"]:
        if f not in df.columns: continue
        per_chrom = df.groupby("chrom")[f].first()
        print(f"\n{f}:", flush=True)
        print(per_chrom.to_string(), flush=True)


if __name__ == "__main__":
    main()
