#!/usr/bin/env python3
"""Step 2: refit every MAC8 discovery with and without the mean-imputed genotype calls."""
import glob
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
TRA = f"{M}/results/meqtl_trans/chm13_observed_sv_all_family_pc10_mac8"
U = (CENSV_ROOT + "/analyses/"
     "meqtl_replication")
OUT = f"{U}/output/arms/refit_genotype_complete"
N_COV_PC10 = 16


def rint(a):
    a = np.asarray(a, float)
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


P = pd.read_csv(f"{OUT}/r1_pairs.tsv", sep="\t")
z = np.load(f"{TRA}/manifests/profile_dosages.npz", allow_pickle=True)
D, prof_ids, dos_samples = z["dosages"], list(z["profile_ids"]), list(z["sample_ids"])
prow = {p: i for i, p in enumerate(prof_ids)}

dc = np.load(f"{M}/results/meqtl_cis/chm13_observed_sv_robust_mac8_pc10/shared/"
             f"canonical_design.npz", allow_pickle=True)
dt = np.load(f"{M}/results/meqtl_hierarchical_pc_sensitivity/shared/canonical_pc_design.npz",
             allow_pickle=True)
COV_CIS = dc["covariates"][:, :N_COV_PC10]
COV_TRA = dt["covariates"][:, :N_COV_PC10]
for nm, dd in (("cis", dc), ("trans", dt)):
    assert list(dd["sample_ids"]) == dos_samples, f"{nm} design sample order must match dosages"
assert list(dc["covariate_names"])[:N_COV_PC10] == list(dt["covariate_names"])[:N_COV_PC10]
print(f"dosages {D.shape}  cis cov {COV_CIS.shape}  trans cov {COV_TRA.shape}")
print(f"designs identical? {np.allclose(COV_CIS, COV_TRA)}"
      f"  max abs diff {np.abs(COV_CIS - COV_TRA).max():.4g}")

meth = {}
files = sorted(glob.glob(f"{OUT}/extracted/*/*.tsv"))
print(f"reading {len(files)} extracted shards")
for f in files:
    t = pd.read_csv(f, sep="\t", low_memory=False)
    if not len(t):
        continue
    cols = [c for c in t.columns if c in set(dos_samples)]
    sub = t.set_index("feature_id")[cols].apply(pd.to_numeric, errors="coerce")
    for fid, row in sub.iterrows():
        meth[fid] = row.reindex(dos_samples).values.astype(float)
print(f"CpGs recovered: {len(meth)} of {P.feature_id.nunique()} wanted")


def fit(y, g, mask, COV):
    n = int(mask.sum())
    if n < 40:
        return None
    gg = g[mask]
    if np.std(gg) < 1e-9 or len(np.unique(gg)) < 2:
        return None
    X = sm.add_constant(np.column_stack([gg, COV[mask]]), has_constant="add")
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    r = sm.OLS(rint(y[mask]), X).fit(cov_type="HC3")
    return dict(n=n, beta=float(r.params[1]), se=float(r.bse[1]),
                t=float(r.tvalues[1]), p=float(r.pvalues[1]))


rows = []
for _, r in P.iterrows():
    y = meth.get(r.feature_id)
    if y is None:
        rows.append(dict(r, status="cpg_not_extracted"))
        continue
    g = D[prow[r.profile_id]]
    obs = np.isin(g, [0.0, 1.0, 2.0])
    mA = np.isfinite(y)
    mB = mA & obs
    CV = COV_TRA if r.arm == "trans" else COV_CIS
    a, b = fit(y, g, mA, CV), fit(y, g, mB, CV)
    rec = dict(r, status="ok", n_imputed_dropped=int((mA & ~obs).sum()))
    for tag, res in (("A", a), ("B", b)):
        if res is None:
            rec[f"{tag}_n"] = np.nan
            for k in ("beta", "se", "t", "p"):
                rec[f"{tag}_{k}"] = np.nan
        else:
            rec[f"{tag}_n"] = res["n"]
            for k in ("beta", "se", "t", "p"):
                rec[f"{tag}_{k}"] = res[k]
    rows.append(rec)

R = pd.DataFrame(rows)
R["retention"] = R.B_beta / R.A_beta
R["sign_kept"] = np.sign(R.B_beta) == np.sign(R.A_beta)
R.to_csv(f"{OUT}/r4_refit.tsv", sep="\t", index=False, float_format="%.6g")

ok = R[(R.status == "ok") & R.A_beta.notna()]
dev = (ok.A_beta - ok.pub_rint_beta).abs()
rho = np.corrcoef(ok.A_beta, ok.pub_rint_beta)[0, 1] if len(ok) > 2 else np.nan
print("\n=== INSTRUMENT CHECK: does fit A reproduce the published betas? ===")
print(f"  pairs compared      : {len(ok)}")
print(f"  Pearson r           : {rho:.6f}")
print(f"  median |deviation|  : {dev.median():.3e}")
print(f"  max |deviation|     : {dev.max():.3e}")
print(f"  sign agreement      : {(np.sign(ok.A_beta)==np.sign(ok.pub_rint_beta)).mean():.4f}")
print(f"  n agreement (A_n==pub_n): {(ok.A_n==ok.pub_n).mean():.4f}")
verdict = "PASS" if (rho > 0.999 and dev.median() < 1e-3) else "FAIL"
print(f"  VERDICT: {verdict}"
      + ("" if verdict == "PASS" else "  -- fit B is NOT interpretable; stop here."))

print("\n=== FIT B: genotype-complete refit, by arm ===")
S = []
for arm, g in R[R.status == "ok"].groupby("arm"):
    gg = g[g.A_beta.notna() & g.B_beta.notna()]
    aff = gg[gg.n_imputed_dropped > 0]
    S.append(dict(
        arm=arm, n_pairs=len(gg),
        n_affected=len(aff),
        median_dropped=float(aff.n_imputed_dropped.median()) if len(aff) else 0.0,
        median_retention=float(gg.retention.median()),
        median_retention_affected=float(aff.retention.median()) if len(aff) else np.nan,
        frac_sign_kept=float(gg.sign_kept.mean()),
        frac_p_lt_1e4=float((gg.B_p < 1e-4).mean()),
        frac_p_lt_005=float((gg.B_p < 0.05).mean()),
        median_absB_over_absA=float((gg.B_beta.abs() / gg.A_beta.abs()).median())))
SS = pd.DataFrame(S).sort_values("arm")
SS.to_csv(f"{OUT}/r5_summary.tsv", sep="\t", index=False, float_format="%.4g")
print(SS.to_string(index=False))

print("\n=== the single SV that carries the distal/trans layer ===")
top = R[(R.family_id == "chr17_26788109_26788109_INS_083a35b6c0") & (R.status == "ok")]
if len(top):
    print(top[["arm", "feature_id", "n_imputed_dropped", "A_n", "A_beta", "A_p",
               "B_n", "B_beta", "B_p", "retention"]].to_string(index=False))
print(f"\nwrote {OUT}/r4_refit.tsv and r5_summary.tsv")
