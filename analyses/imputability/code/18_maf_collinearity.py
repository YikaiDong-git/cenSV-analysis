#!/usr/bin/env python3
"""Phase D verification — quantify collinearity among MAF-related features."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"


def main():
    df = pd.read_csv(OUT / "per_sv_features_with_genome.tsv",
                      sep="\t", low_memory=False)
    cols = ["n_hom_alt", "n_carriers", "log_maf", "maf",
            "total_af_sv", "af_skew_sv", "max_FST_pair_sv"]
    cols = [c for c in cols if c in df.columns]
    sub = df[cols].dropna()
    print(f"[INFO] {len(sub)} SVs with all MAF-family features non-NA",
          flush=True)
    sp = sub.corr(method="spearman").round(3)
    pe = sub.corr(method="pearson").round(3)
    print("\nSpearman ρ matrix:")
    print(sp.to_string())
    print("\nPearson r matrix:")
    print(pe.to_string())

    out_rows = []
    for i in cols:
        for j in cols:
            if i < j:
                out_rows.append({
                    "feat_a": i, "feat_b": j,
                    "spearman_rho": float(sp.loc[i, j]),
                    "pearson_r": float(pe.loc[i, j]),
                })
    pd.DataFrame(out_rows).to_csv(OUT / "v3_maf_collinearity.tsv",
                                    sep="\t", index=False)
    print(f"\n[OK] wrote {OUT/'v3_maf_collinearity.tsv'}", flush=True)


if __name__ == "__main__":
    main()
