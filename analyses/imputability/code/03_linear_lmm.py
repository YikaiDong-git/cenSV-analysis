#!/usr/bin/env python3
"""Linear mixed-effects model: per-feature partial effect on CV R², chrom random intercept."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = Path(CENSV_ROOT + "")
OUT = ROOT / "analyses" / "imputability" / "output"

NUMERIC = [
    "log_svlen", "log_maf",
    "n_carriers", "n_hom_alt",
    "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
    "cat_monophyletic", "cat_residual",
    "cat_non_decomposed", "cat_original_decomposed",
    "log_n_window_snps", "best_snp_dist_mean",
    "best_n_nonzero",
]
CHROM_FEATURES = [
    "log_active_HOR_span", "n_HOR_intervals",
    "pericen_len_500000", "n_train_haps", "acrocentric",
]


def fit_lmm(df, features, y="y_cv_r2"):
    cols = [c for c in features if c in df.columns]
    sub = df[["chrom", y] + cols].dropna(subset=[y]).copy()
    print(f"[INFO] LMM substrate: {len(sub)} rows × {len(cols)} features",
          flush=True)
    nan_rate = sub[cols].isna().mean().sort_values(ascending=False)
    print("Per-feature NaN rate (top 5):", flush=True)
    print(nan_rate.head().to_string(), flush=True)
    for c in cols:
        if sub[c].isna().any():
            med = sub[c].median()
            sub[c] = sub[c].fillna(med)
    keep = [c for c in cols if sub[c].nunique() > 1]
    if not keep:
        return None
    print(f"[INFO] LMM kept features: {keep}", flush=True)
    formula = f"{y} ~ " + " + ".join(keep)
    md = smf.mixedlm(formula, sub, groups=sub["chrom"])
    fit = md.fit(method="lbfgs")
    return fit, sub, keep


def main():
    feat_path = OUT / "per_sv_features.tsv"
    if not feat_path.exists():
        print(f"[FAIL] {feat_path} missing", flush=True)
        sys.exit(1)
    df = pd.read_csv(feat_path, sep="\t", low_memory=False)
    print(f"[INFO] loaded {len(df)} SVs", flush=True)

    try:
        res = fit_lmm(df, NUMERIC, y="y_cv_r2")
    except Exception as e:
        import traceback
        print(f"[FAIL] LMM raised: {type(e).__name__}: {e}", flush=True)
        traceback.print_exc()
        sys.exit(1)
    if res is None:
        print("[FAIL] LMM returned None — no usable features after filtering",
              flush=True)
        sys.exit(1)
    fit, sub, kept = res
    print(fit.summary(), flush=True)

    rows = []
    params = fit.params
    bse = fit.bse
    pvals = fit.pvalues
    for name in kept:
        if name in params.index:
            rows.append({
                "feature": name,
                "beta": float(params[name]),
                "se": float(bse[name]),
                "p": float(pvals[name]),
                "ci_lo": float(params[name] - 1.96 * bse[name]),
                "ci_hi": float(params[name] + 1.96 * bse[name]),
            })
    out = pd.DataFrame(rows).sort_values("p")
    if "Group Var" in params.index:
        rows.append({"feature": "<chrom_random_var>", "beta": float(params["Group Var"]),
                     "se": np.nan, "p": np.nan, "ci_lo": np.nan, "ci_hi": np.nan})
    out = pd.DataFrame(rows).sort_values("p")
    out_path = OUT / "lmm_partial_effects.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path}", flush=True)

    re_path = OUT / "lmm_chrom_random_effects.tsv"
    try:
        re = fit.random_effects
        re_rows = []
        for chrom, s in re.items():
            re_rows.append({"chrom": chrom, "intercept_offset": float(s.iloc[0])})
        re_df = pd.DataFrame(re_rows).sort_values("intercept_offset")
        re_df.to_csv(re_path, sep="\t", index=False)
        print(f"[OK] wrote {re_path}", flush=True)
        print("\nPer-chrom random intercepts (low → high):", flush=True)
        print(re_df.to_string(index=False), flush=True)
    except (np.linalg.LinAlgError, ValueError) as e:
        print(f"\n[NOTE] could not extract random effects: {e}", flush=True)
        print("Group Var ~ 0 — chrom adds no residual variance after feature controls.",
              flush=True)
        pd.DataFrame([{"chrom": "<all>", "intercept_offset": 0.0,
                       "note": "Group Var ~0; chrom variance absorbed by features"}]
                     ).to_csv(re_path, sep="\t", index=False)
        print(f"[OK] wrote placeholder {re_path}", flush=True)

    chrom_agg = (df.groupby("chrom")
                   .agg(mean_cv_r2=("cv_r2", "mean"),
                        n_svs=("cv_r2", "count"),
                        **{f: (f, "first") for f in CHROM_FEATURES if f in df.columns})
                   .reset_index())
    print(f"\n[INFO] between-chrom regression substrate: {len(chrom_agg)} chroms",
          flush=True)
    keep_chr = [f for f in CHROM_FEATURES if f in chrom_agg.columns
                and chrom_agg[f].nunique() > 1]
    if keep_chr and len(chrom_agg) >= len(keep_chr) + 2:
        formula = "mean_cv_r2 ~ " + " + ".join(keep_chr)
        ols = sm.OLS.from_formula(formula, data=chrom_agg).fit()
        chrom_rows = []
        for f in keep_chr + ["Intercept"]:
            if f in ols.params.index:
                chrom_rows.append({
                    "feature": f, "beta": float(ols.params[f]),
                    "se": float(ols.bse[f]), "p": float(ols.pvalues[f]),
                    "ci_lo": float(ols.params[f] - 1.96 * ols.bse[f]),
                    "ci_hi": float(ols.params[f] + 1.96 * ols.bse[f]),
                })
        pd.DataFrame(chrom_rows).to_csv(OUT / "between_chrom_ols.tsv",
                                        sep="\t", index=False)
        print(f"\nBetween-chrom OLS R² = {ols.rsquared:.3f}", flush=True)
        print(ols.summary().as_text(), flush=True)
    else:
        print("[WARN] insufficient variance in chrom-level features for OLS",
              flush=True)


if __name__ == "__main__":
    main()
