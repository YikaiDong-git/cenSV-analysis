#!/usr/bin/env python3
"""Enlarged associations: the 107 AoU-testable SV x trait pairs (of 222), EUR + all-ancestry."""
import os, sys, gc, gzip, csv, numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
H = f"{BASE}/phewide_222_handoff"; R6 = f"{BASE}/aou_replication/results_vds/step6_222"; S3 = f"{BASE}/aou_replication/results/step3"
OUT = f"{BASE}/aou_replication/results_vds/step8_222"; os.makedirs(OUT, exist_ok=True)
COVDIR = os.path.expanduser("~/workspace/working/scripts/covariates/output")
DENOM = os.path.expanduser("~/workspace/working/PRS/results/phewas/aou_ehr_icd_denominator.tsv")


def ols_coef(X, y):
    XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ (X.T @ y); resid = y - X @ beta
    n, p = X.shape; s2 = float(resid @ resid) / (n - p); se = np.sqrt(np.diag(XtXi) * s2)
    return beta[1], se[1], 2 * tdist.sf(abs(beta[1] / se[1]), n - p)


def logit_coef(X, y, maxit=100, tol=1e-7):
    beta = np.zeros(X.shape[1]); converged = False
    for _ in range(maxit):
        mu = 1 / (1 + np.exp(-(X @ beta))); W = np.clip(mu * (1 - mu), 1e-9, None)
        d = np.linalg.solve((X * W[:, None]).T @ X, X.T @ (y - mu)); beta += d
        if np.max(np.abs(d)) < tol:
            converged = True; break
    mu = 1 / (1 + np.exp(-(X @ beta)))
    se = np.sqrt(np.diag(np.linalg.inv((X * (mu * (1 - mu))[:, None]).T @ X)))
    return beta[1], se[1], 2 * norm.sf(abs(beta[1] / se[1])), converged


def irnt(x):
    m = np.isfinite(x); out = np.full(x.shape, np.nan); r = rankdata(x[m]); out[m] = norm.ppf((r - 0.5) / m.sum()); return out


tmap = {r['ukb_trait']: (r['type'], r['loinc']) for r in csv.DictReader(open(f"{R6}/testable_trait_map.tsv"), delimiter='\t')}
pairs = [p for p in csv.DictReader(open(f"{H}/pairs_222.tsv"), delimiter='\t') if p['trait'] in tmap]
sv_order = list(csv.DictReader(open(f"{R6}/dosage_qc_222.tsv"), delimiter='\t'))
sv_order = [r['sv_locus'] for r in sv_order]
quant_loincs = sorted({tmap[p['trait']][1] for p in pairs if tmap[p['trait']][0] == 'quant'})
DISEASE = {'I9_VTE': 'vte', 'K11_CHOLELITH': 'cholelithiasis'}

dtm = {(r['vocabulary_id'], r['concept_code']): r['phenotype_id'] for r in csv.DictReader(open(f"{S3}/target_phecode_icd_map.3dz.tsv"), delimiter='\t')}
dcases = {'vte': set(), 'cholelithiasis': set()}
for r in csv.DictReader(open(f"{S3}/icd_evidence.3dz.tsv"), delimiter='\t'):
    ph = dtm.get((r['vocabulary_id'], r['concept_code']))
    if ph in dcases:
        dcases[ph].add(r['IID'])
denom = set(l.strip() for l in open(DENOM)); denom.discard('IID')


def run_mode(mode, covfile):
    base = ['age_at_recruitment_spline_1', 'age_at_recruitment_spline_2', 'age_at_recruitment_spline_3', 'sex', *[f'gPC{i}' for i in range(1, 17)]]
    use = ['IID', *base] + (['ancestry_pred'] if mode == 'all' else [])
    cov = pd.read_csv(covfile, sep='\t', dtype={'IID': str}, usecols=lambda c: c in use)
    cov['sex'] = cov['sex'].map({'Male': 1.0, 'Female': 0.0})
    for c in base:
        cov[c] = pd.to_numeric(cov[c], errors='coerce')
    cov = cov.dropna(subset=base).reset_index(drop=True)
    iids = cov['IID'].tolist(); idx = {i: k for k, i in enumerate(iids)}; N = len(iids)
    C = cov[base].to_numpy(np.float64)
    if mode == 'all':
        C = np.hstack([C, pd.get_dummies(cov['ancestry_pred'], drop_first=True).astype(np.float64).to_numpy()])
    del cov; gc.collect()
    D = np.full((N, len(sv_order)), np.nan, np.float32)
    with gzip.open(f"{R6}/sv_dosages_222.tsv.gz", 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t'); i = idx.get(p[0])
            if i is not None:
                D[i, :] = [float(x) for x in p[1:]]
    svcol = {s: k for k, s in enumerate(sv_order)}
    lv = {lc: np.full(N, np.nan, np.float32) for lc in quant_loincs}
    for r in csv.DictReader(open(f"{R6}/measurement_medians_222.tsv"), delimiter='\t'):
        i = idx.get(r['IID'])
        if i is not None and r['loinc'] in lv:
            lv[r['loinc']][i] = float(r['median_value'])
    lirnt = {lc: irnt(v.astype(float)) for lc, v in lv.items()}; del lv; gc.collect()
    dz = {}
    for d, cs in dcases.items():
        y = np.full(N, np.nan)
        for i, k in idx.items():
            if i in denom:
                y[k] = 1.0 if i in cs else 0.0
            elif i in cs:
                y[k] = 1.0
        dz[d] = y
    res = []
    for p in pairs:
        sv, trait = p['sv_locus'], p['trait']; typ, loinc = tmap[trait]
        x = D[:, svcol[sv]].astype(np.float64)
        if typ == 'disease':
            y = dz[DISEASE[trait]]; model = 'logistic'
        else:
            y = lirnt[loinc]; model = 'linear'
        mask = np.isfinite(y) & np.isfinite(x) & np.isfinite(C).all(1)
        n = int(mask.sum())
        if n < 100:
            res.append([sv, trait, p['tier'], p['cv_r2'], p['discovery_sign'], mode, model, n, '', '', '', '', 'n<100']); continue
        X = np.column_stack([np.ones(n), x[mask], C[mask]]).astype(np.float64)
        keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9; X = X[:, keep]
        yy = y[mask].astype(np.float64); ncase = int(yy.sum()) if model == 'logistic' else ''
        try:
            if model == 'linear':
                b, se, pv = ols_coef(X, yy)
            else:
                b, se, pv, conv = logit_coef(X, yy)
                if not conv:
                    res.append([sv, trait, p['tier'], p['cv_r2'], p['discovery_sign'], mode, model, n, ncase, f"{b:.5g}", f"{se:.5g}", f"{pv:.3g}", 'logit_not_converged']); continue
            asign = '+' if b > 0 else '-'
            res.append([sv, trait, p['tier'], p['cv_r2'], p['discovery_sign'], mode, model, n, ncase, f"{b:.5g}", f"{se:.5g}", f"{pv:.3g}", 'concordant' if asign == p['discovery_sign'] else 'discordant'])
        except Exception as e:
            res.append([sv, trait, p['tier'], p['cv_r2'], p['discovery_sign'], mode, model, n, ncase, '', '', '', f"err:{type(e).__name__}"])
        del X, yy
    del D, C, lirnt, dz; gc.collect()
    return res


HDR = "sv_locus\ttrait\ttier\tcv_r2\tdiscovery_sign\tancestry_mode\tmodel\tn\tn_cases\tbeta\tse\tp\tconcordance\n"
MODEFILE = {'eur': f"{COVDIR}/covariates_v3_eur.tsv", 'all': f"{COVDIR}/covariates_v3_all_ancestries.tsv"}
for mode in ([sys.argv[1]] if len(sys.argv) > 1 else ['eur', 'all']):
    print(f"=== {mode} ({len(pairs)} pairs) ===", flush=True)
    res = run_mode(mode, MODEFILE[mode])
    with open(f"{OUT}/results_222_{mode}.tsv", 'w') as f:
        f.write(HDR)
        for row in res:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"  wrote results_222_{mode}.tsv ({len(res)} pairs)", flush=True)
    gc.collect()
