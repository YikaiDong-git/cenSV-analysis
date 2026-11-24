#!/usr/bin/env python3
"""Disease/logistic-path positive control: Factor V Leiden (rs6025) -> VTE."""
import os, gzip, csv, numpy as np, pandas as pd
import statsmodels.api as sm
import pgenlib

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
S3 = f"{BASE}/aou_replication/results/step3"
PGEN = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen")
COVDIR = os.path.expanduser("~/workspace/working/scripts/covariates/output")
DENOM = os.path.expanduser("~/workspace/working/PRS/results/phewas/aou_ehr_icd_denominator.tsv")

POS, REF, ALT = '169549811', 'C', 'T'
offsets = [0]; target = None; vidx = 0
with open(f"{PGEN}/acaf_threshold.chr1.pvar") as f:
    for line in f:
        if line[0] == '#':
            continue
        c = line.rstrip('\n').split('\t', 5)
        alts = c[4].split(',')
        offsets.append(offsets[-1] + 1 + len(alts))
        if c[1] == POS and c[3] == REF and ALT in alts:
            target = (vidx, 1 + alts.index(ALT))
        vidx += 1
assert target is not None, "Leiden variant not found"
offsets = np.array(offsets, dtype=np.uintp)
iids = [l.rstrip('\n').split('\t')[0] for l in open(f"{PGEN}/acaf_threshold.chr1.psam")][1:]
N = len(iids)
r = pgenlib.PgenReader(f"{PGEN}/acaf_threshold.chr1.pgen".encode(), raw_sample_ct=N, variant_ct=vidx, allele_idx_offsets=offsets)
buf = np.empty(N, dtype=np.int32)
r.read(target[0], buf, allele_idx=target[1]); r.close()
leiden = buf.astype(np.float64); leiden[leiden < 0] = np.nan
freq = np.nanmean(leiden) / 2
print(f"Leiden T-allele freq (all samples): {freq:.4f}  (expect ~0.02-0.05 in EUR)")
leiden_by_iid = dict(zip(iids, leiden))

tm = {(x['vocabulary_id'], x['concept_code']): x['phenotype_id'] for x in csv.DictReader(open(f"{S3}/target_phecode_icd_map.3dz.tsv"), delimiter='\t')}
vte = set()
for x in csv.DictReader(open(f"{S3}/icd_evidence.3dz.tsv"), delimiter='\t'):
    if tm.get((x['vocabulary_id'], x['concept_code'])) == 'vte':
        vte.add(x['IID'])
denom = set(l.strip() for l in open(DENOM)); denom.discard('IID')

base = ['age_at_recruitment_spline_1', 'age_at_recruitment_spline_2', 'age_at_recruitment_spline_3', 'sex', *[f'gPC{i}' for i in range(1, 17)]]
for mode, cf in [('eur', f"{COVDIR}/covariates_v3_eur.tsv"), ('all', f"{COVDIR}/covariates_v3_all_ancestries.tsv")]:
    use = ['IID', *base] + (['ancestry_pred'] if mode == 'all' else [])
    cov = pd.read_csv(cf, sep='\t', dtype={'IID': str}, usecols=lambda c: c in use)
    cov['sex'] = cov['sex'].map({'Male': 1.0, 'Female': 0.0})
    for c in base:
        cov[c] = pd.to_numeric(cov[c], errors='coerce')
    cov['leiden'] = cov['IID'].map(leiden_by_iid)
    cov['y'] = cov['IID'].map(lambda i: 1.0 if i in vte else (0.0 if i in denom else np.nan))
    dcols = []
    if mode == 'all':
        dum = pd.get_dummies(cov['ancestry_pred'], drop_first=True).astype(float)
        dcols = list(dum.columns); cov = pd.concat([cov, dum], axis=1)
    cov = cov.dropna(subset=['leiden', 'y'] + base)
    y = cov['y'].to_numpy(float)
    raw = cov['leiden'].to_numpy(float)
    z = (raw - raw.mean()) / raw.std()
    X = np.column_stack([np.ones(len(y)), z, cov[base + dcols].to_numpy(float)])
    keep = np.ones(X.shape[1], bool)
    keep[2:] = X[:, 2:].std(0) > 1e-9
    X = X[:, keep]
    fit = sm.Logit(y, X).fit(disp=0, maxiter=200)
    b_sd, se_sd, p = fit.params[1], fit.bse[1], fit.pvalues[1]
    b_allele = b_sd / raw.std()
    print(f"{mode}: n={len(y)} n_cases={int(y.sum())} | per-SD logOR={b_sd:.3f} (OR/SD={np.exp(b_sd):.2f}) "
          f"| per-ALLELE OR={np.exp(b_allele):.2f} | p={p:.3g}")
