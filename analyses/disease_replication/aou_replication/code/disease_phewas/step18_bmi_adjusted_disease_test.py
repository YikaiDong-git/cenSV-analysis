#!/usr/bin/env python3
"""Step 18: FORMAL vertical-pleiotropy test for the two 'clean' disease hits."""
import os, gzip, csv, gc, numpy as np, pandas as pd
from scipy.stats import norm, rankdata

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds"
PHENO = f"{HOME}/workspace/working/PRS/results/phewas/phenotypes.tsv"
COV = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_eur.tsv"
MEAS = f"{RV}/step6_222/measurement_medians_222.tsv"
DOSFILES = [f"{RV}/step6_222/sv_dosages_222.tsv.gz", f"{RV}/step2/sv_dosages.tsv.gz"]
OUT = f"{RV}/step18_bmi_adjusted"; os.makedirs(OUT, exist_ok=True)

TESTS = [("chr7_sv114_DEL_non_best8_edge1000", "cerebrovascular_disease"),
         ("chr11_sv121_DEL_non_best8_fksize100k", "hypertension")]
BMI_LOINC, WEIGHT_LOINC = "39156-5", "29463-7"


def logit_coef(X, y, maxit=200, tol=1e-8):
    beta = np.zeros(X.shape[1]); conv = False
    for _ in range(maxit):
        mu = 1 / (1 + np.exp(-(X @ beta))); W = np.clip(mu * (1 - mu), 1e-9, None)
        try:
            d = np.linalg.solve((X * W[:, None]).T @ X, X.T @ (y - mu))
        except np.linalg.LinAlgError:
            return np.nan, np.nan, np.nan, False
        beta += d
        if np.max(np.abs(d)) < tol:
            conv = True; break
    mu = 1 / (1 + np.exp(-(X @ beta)))
    se = np.sqrt(np.diag(np.linalg.inv((X * (mu * (1 - mu))[:, None]).T @ X)))
    return beta[1], se[1], 2 * norm.sf(abs(beta[1] / se[1])), conv


def irnt(x):
    m = np.isfinite(x); out = np.full(x.shape, np.nan)
    r = rankdata(x[m]); out[m] = norm.ppf((r - 0.5) / m.sum()); return out


base = ['age_at_recruitment_spline_1', 'age_at_recruitment_spline_2', 'age_at_recruitment_spline_3',
        'sex', *[f'gPC{i}' for i in range(1, 17)]]
cov = pd.read_csv(COV, sep='\t', dtype={'IID': str}, usecols=lambda c: c in ['IID', *base])
cov['sex'] = cov['sex'].map({'Male': 1.0, 'Female': 0.0})
for c in base:
    cov[c] = pd.to_numeric(cov[c], errors='coerce')
cov = cov.dropna(subset=base).reset_index(drop=True)
iids = cov['IID'].tolist(); idx = {i: k for k, i in enumerate(iids)}; N = len(iids)
C = cov[base].to_numpy(np.float64); del cov; gc.collect()

ph = pd.read_csv(PHENO, sep='\t', dtype=str)
dz = {d: np.full(N, np.nan) for _, d in TESTS}
for d in dz:
    col = pd.to_numeric(ph[d], errors='coerce').to_numpy()
    for k, i in enumerate(ph['IID'].tolist()):
        j = idx.get(i)
        if j is not None:
            dz[d][j] = col[k]
del ph; gc.collect()

want = {sv for sv, _ in TESTS}
dos = {sv: np.full(N, np.nan, np.float64) for sv in want}
got = set()
for path in DOSFILES:
    with gzip.open(path, 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        pos = {sv: hdr.index(sv) for sv in want if sv in hdr and sv not in got}
        if not pos:
            continue
        got |= set(pos)
        for line in f:
            p = line.rstrip('\n').split('\t'); j = idx.get(p[0])
            if j is not None:
                for sv, c in pos.items():
                    dos[sv][j] = float(p[c])

raw = {BMI_LOINC: np.full(N, np.nan), WEIGHT_LOINC: np.full(N, np.nan)}
for r in csv.DictReader(open(MEAS), delimiter='\t'):
    j = idx.get(r['IID'])
    if j is not None and r['loinc'] in raw:
        raw[r['loinc']][j] = float(r['median_value'])
BMIi = irnt(raw[BMI_LOINC]); WTi = irnt(raw[WEIGHT_LOINC])
Cmask = np.isfinite(C).all(1)


def fit_fixed(sv, disease, adjusters, force_finite):
    x = dos[sv]; y = dz[disease]
    m = np.isfinite(x) & np.isfinite(y) & Cmask
    for _, v in force_finite:
        m &= np.isfinite(v)
    n = int(m.sum()); ncase = int(np.nansum(y[m]))
    xz = (x[m] - x[m].mean()) / x[m].std()
    cols = [np.ones(n), xz, C[m]] + [v[m] for _, v in adjusters]
    X = np.column_stack(cols).astype(np.float64)
    keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9; X = X[:, keep]
    b, se, pv, conv = logit_coef(X, y[m].astype(np.float64))
    return dict(n=n, ncase=ncase, beta=b, se=se, p=pv, OR=np.exp(b) if np.isfinite(b) else np.nan, conv=conv)


rows = []
for sv, disease in TESTS:
    both = [("BMI", BMIi), ("weight", WTi)]
    models = [("unadjusted", []), ("+BMI", [("BMI", BMIi)]),
              ("+weight", [("weight", WTi)]), ("+BMI+weight", both)]
    ref = None
    for label, adj in models:
        r = fit_fixed(sv, disease, adj, both)
        if label == "unadjusted":
            ref = r["beta"]
        att = (1 - r["beta"] / ref) * 100 if (ref and np.isfinite(r["beta"])) else float('nan')
        rows.append(['_'.join(sv.split('_')[:2]), disease, label, r["n"], r["ncase"],
                     f"{r['OR']:.4f}", f"{r['beta']:.5g}", f"{r['se']:.5g}", f"{r['p']:.3g}",
                     f"{att:+.1f}%", r["conv"]])


HDR = "sv\tdisease\tmodel\tn\tn_cases\tOR_perSD\tbeta\tse\tp\tdosage_attenuation_vs_unadj\tconverged\n"
with open(f"{OUT}/bmi_adjusted_test_eur.tsv", 'w') as f:
    f.write(HDR)
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
print(f"wrote {OUT}/bmi_adjusted_test_eur.tsv ({len(rows)} rows)")
for r in rows:
    print("\t".join(str(x) for x in r))
