#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, t as tdist
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

SUB = CENSV_ROOT + "/subproject_methylation"
DISC = CENSV_ROOT + "/analyses/cdr_qtl"
INP = f"{SUB}/results/meqtl_hprc_direct/inputs"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"

LEADS = [
    ("chr3", 2240, "INS339 f0.19 (manuscript cleanest)", "cdrs_per_mb", None),
    ("chr3", 2241, "INS339 f0.30", "cdrs_per_mb", None),
    ("chr3", 2238, "INS338 f0.86 (dominant)", "cdr_fraction", None),
    ("chr6", 964,  "INS337 f0.18 (pop-strat)", "cdr_pos_offset", 2.7347e-12),
    ("chr6", 2900, "INS506 f0.26 (pop-strat)", "cdr_pos_var", 3.1857e-10),
    ("chr8", 988,  "DEL342 f0.24", "cdr_fraction", 7.415e-10),
    ("chr8", 821,  "DEL1356 f0.22", "cdr_fraction", 1.0258e-09),
    ("chr17", 4944, "INS167 f0.32", "cdrs_per_mb", 2.0249e-09),
]
PHENOS = ["mean_5mC", "n_cdrs", "cdrs_per_mb", "cdr_fraction",
          "cdr_pos_offset", "cdr_pos_var", "cdr_pos_iqr"]


def rint(x):
    x = np.asarray(x, dtype=float)
    mask = ~np.isnan(x)
    r = np.full_like(x, np.nan, dtype=float)
    if mask.sum() == 0:
        return r
    r[mask] = norm.ppf((rankdata(x[mask]) - 0.5) / mask.sum())
    return r


def ols_test(y, dos, bc):
    m = ~(np.isnan(y) | np.isnan(dos) | np.isnan(bc))
    n = int(m.sum())
    if n < 30 or dos[m].std() == 0 or y[m].std() == 0:
        return None
    if bc[m].std() == 0:
        X = np.column_stack([np.ones(n), dos[m]])
        k = 2
    else:
        X = np.column_stack([np.ones(n), dos[m], bc[m]])
        k = 3
    Y = y[m]
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    rss = (resid ** 2).sum()
    sigma2 = rss / (n - k)
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    tval = beta[1] / se[1]
    p = 2 * (1 - tdist.cdf(abs(tval), n - k))
    r2 = 1 - rss / ((Y - Y.mean()) ** 2).sum()
    return dict(n=n, beta=float(beta[1]), t=float(tval), p=float(p), r2=float(r2))


def main():
    ph = pd.read_csv(PHENO, sep="\t")
    indiv = ph.groupby(["sample", "chm13_chrom"], as_index=False).agg(
        **{p: (p, "mean") for p in PHENOS},
        basecaller=("basecaller", lambda x: x.iloc[0]),
    )
    indiv["bc_R1041"] = (indiv["basecaller"] == "R1041").astype(float)
    print(f"[pheno] {len(indiv)} (sample,chrom) rows; "
          f"R1041 frac = {indiv['bc_R1041'].mean():.2f}")

    out_rows = []
    for chrom, sv_index, label, prior_ph, prior_p in LEADS:
        dos_df = pd.read_csv(f"{INP}/{chrom}/{chrom}_imputed_dosages.tsv", sep="\t")
        dos_df["sv_index"] = pd.to_numeric(dos_df["sv_index"], errors="coerce")
        row = dos_df[dos_df["sv_index"] == sv_index]
        if row.empty:
            print(f"!! {chrom} sv {sv_index} not in dosages")
            continue
        sample_cols = [c for c in dos_df.columns if c != "sv_index"]
        dser = row[sample_cols].iloc[0].astype(float)

        sub = indiv[indiv["chm13_chrom"] == chrom].set_index("sample")
        dos = np.array([dser[s] if s in dser.index else np.nan for s in sub.index])
        bc = sub["bc_R1041"].values
        print(f"\n=== {chrom} sv{sv_index} {label}  (prior {prior_ph}, prior_p="
              f"{('%.2e' % prior_p) if prior_p else 'NA'}) ===")
        for cohort, mask_bc in [("ALL", np.ones(len(sub), bool)),
                                ("R941only", (bc == 0))]:
            for pheno in PHENOS:
                y = rint(np.where(mask_bc, sub[pheno].values, np.nan))
                res = ols_test(y, np.where(mask_bc, dos, np.nan), bc)
                if res is None:
                    continue
                star = ""
                if pheno == prior_ph:
                    star = " <-prior"
                if pheno == "mean_5mC":
                    star += " [PRIMARY]"
                flag = "***" if res["p"] < 1e-4 else ("*" if res["p"] < 0.05 else "")
                out_rows.append(dict(chrom=chrom, sv_index=sv_index, label=label,
                                     cohort=cohort, phenotype=pheno, **res,
                                     prior_phenotype=prior_ph, prior_p=prior_p))
                print(f"  [{cohort:8s}] {pheno:16s} n={res['n']:3d} beta={res['beta']:+.3f} "
                      f"t={res['t']:+.2f} p={res['p']:.2e} r2={res['r2']:.3f} {flag}{star}")

    df = pd.DataFrame(out_rows)
    out = f"{DISC}/output/gates/lead_survival_check.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.4g")
    print(f"\n[written] {out}: {len(df)} tests")
    print("\n=== SURVIVAL SUMMARY (min p per lead, ALL cohort) ===")
    for chrom, sv_index, label, prior_ph, prior_p in LEADS:
        g = df[(df.sv_index == sv_index) & (df.chrom == chrom) & (df.cohort == "ALL")]
        if g.empty:
            continue
        best = g.loc[g["p"].idxmin()]
        prior_now = g[g.phenotype == prior_ph]
        pn = prior_now["p"].iloc[0] if len(prior_now) else np.nan
        print(f"  {chrom} sv{sv_index:<5d} prior_p={('%.1e'%prior_p) if prior_p else 'NA':>8s} | "
              f"prior_pheno_now p={pn:.1e} | best={best['phenotype']}({best['p']:.1e})")


if __name__ == "__main__":
    main()
