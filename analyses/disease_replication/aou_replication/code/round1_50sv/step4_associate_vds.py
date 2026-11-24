#!/usr/bin/env python3
"""Steps 3d + 4: assemble phenotypes and run imputed-SV-dosage ~ trait associations."""
import os, sys, gc, gzip, csv, numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist


def ols_coef(X, y):
    XtXi = np.linalg.inv(X.T @ X)
    beta = XtXi @ (X.T @ y)
    resid = y - X @ beta
    n, p = X.shape
    sigma2 = float(resid @ resid) / (n - p)
    se = np.sqrt(np.diag(XtXi) * sigma2)
    tstat = beta[1] / se[1]
    return beta[1], se[1], 2 * tdist.sf(abs(tstat), n - p)


def logit_coef(X, y, maxit=50, tol=1e-7):
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(maxit):
        mu = 1.0 / (1.0 + np.exp(-(X @ beta)))
        W = np.clip(mu * (1 - mu), 1e-9, None)
        XtWX = (X * W[:, None]).T @ X
        delta = np.linalg.solve(XtWX, X.T @ (y - mu))
        beta += delta
        if np.max(np.abs(delta)) < tol:
            break
    mu = 1.0 / (1.0 + np.exp(-(X @ beta)))
    cov = np.linalg.inv((X * (mu * (1 - mu))[:, None]).T @ X)
    se = np.sqrt(np.diag(cov))
    z = beta[1] / se[1]
    return beta[1], se[1], 2 * norm.sf(abs(z))

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
H = f"{BASE}/aou_handoff"; S2 = f"{BASE}/aou_replication/results_vds/step2"; S3 = f"{BASE}/aou_replication/results/step3"
OUT = f"{BASE}/aou_replication/results_vds/step4"; os.makedirs(OUT, exist_ok=True)
COVDIR = os.path.expanduser("~/workspace/working/scripts/covariates/output")
DENOM = os.path.expanduser("~/workspace/working/PRS/results/phewas/aou_ehr_icd_denominator.tsv")

CODE2TRAIT = {'50_irnt':'height','HEIGHT_IRN':'height','21001_irnt':'BMI','BMI_IRN':'BMI','30600_irnt':'albumin',
 '30020_irnt':'haemoglobin','30080_irnt':'platelet_count','30110_irnt':'platelet_dist_width','30130_irnt':'monocyte_count',
 '30180_irnt':'lymphocyte_pct','30190_irnt':'monocyte_pct','30240_irnt':'reticulocyte_pct','30260_irnt':'retic_mean_vol',
 '30630_irnt':'ApoA','30710_irnt':'CRP','30720_irnt':'cystatinC','30770_irnt':'IGF1','30870_irnt':'triglycerides','30880_irnt':'urate'}
DISEASE = {'I9_VTE':'vte','K11_CHOLELITH':'cholelithiasis','E4_DIABINSIPIDUS':'diabetes_insipidus'}

meta = {r['sv_locus']: r for r in csv.DictReader(open(f"{H}/sv_metadata.tsv"), delimiter='\t')}
qc = {r['sv_locus']: r for r in csv.DictReader(open(f"{S2}/dosage_qc.tsv"), delimiter='\t')}
sv_order = list(qc)

tm = {(r['vocabulary_id'], r['concept_code']): r['phenotype_id'] for r in csv.DictReader(open(f"{S3}/target_phecode_icd_map.3dz.tsv"), delimiter='\t')}
disease_cases = {d: set() for d in DISEASE.values()}
for r in csv.DictReader(open(f"{S3}/icd_evidence.3dz.tsv"), delimiter='\t'):
    ph = tm.get((r['vocabulary_id'], r['concept_code']))
    if ph: disease_cases[ph].add(r['IID'])
denom = set(l.strip() for l in open(DENOM))
denom.discard('IID')

def irnt(x):
    m = np.isfinite(x); out = np.full(x.shape, np.nan)
    r = rankdata(x[m]); out[m] = norm.ppf((r - 0.5) / m.sum())
    return out

def run_mode(mode, covfile):
    base_cols = ['age_at_recruitment_spline_1', 'age_at_recruitment_spline_2', 'age_at_recruitment_spline_3', 'sex',
                 *[f'gPC{i}' for i in range(1, 17)]]
    usecols = ['IID', *base_cols] + (['ancestry_pred'] if mode == 'all' else [])
    cov = pd.read_csv(covfile, sep='\t', dtype={'IID': str}, usecols=lambda c: c in usecols)
    cov['sex'] = cov['sex'].map({'Male': 1.0, 'Female': 0.0})
    for c in base_cols:
        cov[c] = pd.to_numeric(cov[c], errors='coerce')
    cov = cov.dropna(subset=base_cols).reset_index(drop=True)
    iids = cov['IID'].tolist(); idx = {iid: i for i, iid in enumerate(iids)}; N = len(iids)
    C = cov[base_cols].to_numpy(np.float64)
    if mode == 'all' and 'ancestry_pred' in cov:
        anc = pd.get_dummies(cov['ancestry_pred'], drop_first=True).astype(np.float64).to_numpy()
        C = np.hstack([C, anc])
    del cov; gc.collect()
    D = np.full((N, len(sv_order)), np.nan, np.float32)
    with gzip.open(f"{S2}/sv_dosages.tsv.gz", 'rt') as f:
        hdr = f.readline().rstrip('\n').split('\t'); svc = {s: k for k, s in enumerate(hdr[1:])}
        for line in f:
            p = line.rstrip('\n').split('\t'); i = idx.get(p[0])
            if i is None: continue
            D[i, :] = [float(x) for x in p[1:]]
    tvals = {t: np.full(N, np.nan, np.float32) for t in set(CODE2TRAIT.values())}
    for r in csv.DictReader(open(f"{S3}/measurement_medians.tsv"), delimiter='\t'):
        i = idx.get(r['IID'])
        if i is not None and r['trait'] in tvals: tvals[r['trait']][i] = float(r['median_value'])
    tirnt = {t: irnt(v) for t, v in tvals.items()}
    dz = {}
    for d, cs in disease_cases.items():
        y = np.full(N, np.nan)
        for iid, i in idx.items():
            if iid in denom: y[i] = 1.0 if iid in cs else 0.0
            elif iid in cs: y[i] = 1.0
        dz[d] = y

    res = []
    for sv in sv_order:
        code = meta[sv]['rep_trait']
        if code in DISEASE:
            y = dz[DISEASE[code]]; model = 'logistic'; trait = DISEASE[code]
        elif code in CODE2TRAIT:
            y = tirnt[CODE2TRAIT[code]]; model = 'linear'; trait = CODE2TRAIT[code]
        else:
            continue
        x = D[:, sv_order.index(sv)].astype(np.float64)
        mask = np.isfinite(y) & np.isfinite(x) & np.isfinite(C).all(1)
        n = int(mask.sum())
        if n < 100:
            res.append((sv, code, trait, mode, model, n, '', '', '', '', 'n<100')); continue
        X = np.column_stack([np.ones(n), x[mask], C[mask]]).astype(np.float64)
        keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9; X = X[:, keep]
        yy = y[mask].astype(np.float64); ncase = int(yy.sum()) if model == 'logistic' else ''
        try:
            beta, se, p = ols_coef(X, yy) if model == 'linear' else logit_coef(X, yy)
            res.append((sv, code, trait, mode, model, n, ncase, f"{beta:.5g}", f"{se:.5g}", f"{p:.3g}", ''))
        except Exception as e:
            res.append((sv, code, trait, mode, model, n, ncase, '', '', '', f"fit_error:{type(e).__name__}"))
        del X, yy; gc.collect()
    del D, C, tvals, tirnt, dz
    return res

HDR = "sv_locus\trep_trait\ttrait\tancestry_mode\tmodel\tn\tn_cases\tbeta\tse\tp\tn_weight_snps_total\tn_weight_snps_in_pgen\tdosage_sd\tin_robust24\tnote\n"
def write_mode(mode, res):
    path = f"{OUT}/aou_results_{mode}.tsv"
    with open(path, 'w') as f:
        f.write(HDR)
        for (sv, code, trait, m, model, n, ncase, beta, se, p, note) in res:
            f.write(f"{sv}\t{code}\t{trait}\t{m}\t{model}\t{n}\t{ncase}\t{beta}\t{se}\t{p}\t{qc[sv]['n_total']}\t{qc[sv]['n_used']}\t{qc[sv]['dosage_sd']}\t{meta[sv]['in_robust24']}\t{note}\n")
    print(f"wrote {path} ({len(res)} rows)", flush=True)

MODEFILE = {'eur': f"{COVDIR}/covariates_v3_eur.tsv", 'all': f"{COVDIR}/covariates_v3_all_ancestries.tsv"}
modes = [sys.argv[1]] if len(sys.argv) > 1 else ['eur', 'all']
for mode in modes:
    print(f"=== mode {mode} ===", flush=True)
    res = run_mode(mode, MODEFILE[mode])
    write_mode(mode, res)
    gc.collect()
