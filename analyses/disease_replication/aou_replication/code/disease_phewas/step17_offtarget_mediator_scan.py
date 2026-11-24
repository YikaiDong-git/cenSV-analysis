#!/usr/bin/env python3
"""Step 17: OFF-TARGET quantitative-mediator scan (EUR, VDS dosages, $0 — no new genotype reads)."""
import os, gc, gzip, csv, numpy as np, pandas as pd
from scipy.stats import rankdata, norm, t as tdist

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
RV = f"{BASE}/aou_replication/results_vds"
R1 = f"{RV}/step2/sv_dosages.tsv.gz"
R2 = f"{RV}/step6_222/sv_dosages_222.tsv.gz"
MEAS = f"{RV}/step6_222/measurement_medians_222.tsv"
H = f"{BASE}/phewide_222_handoff"
COV = os.path.expanduser("~/workspace/working/scripts/covariates/output/covariates_v3_eur.tsv")
OUT = f"{RV}/step17_mediator"; os.makedirs(OUT, exist_ok=True)

MED = {
    '39156-5': 'BMI', '8302-2': 'height', '62409-8': 'waist_circ', '29463-7': 'weight',
    '8480-6': 'SBP', '8462-4': 'DBP', '2345-7': 'glucose', '4548-4': 'HbA1c', '2085-9': 'HDL',
    '2093-3': 'total_cholesterol', '2571-8': 'triglycerides', '3084-1': 'urate',
    '33863-2': 'cystatinC', '3094-0': 'urea', '17861-6': 'calcium',
    '6690-2': 'WBC', '789-8': 'RBC', '731-0': 'lymphocyte_count', '742-7': 'monocyte_count',
    '751-8': 'neutrophil_count', '32623-1': 'mean_platelet_volume',
}
HITS = {'chr7_sv114_DEL_non_best8_edge1000', 'chr11_sv121_DEL_non_best8_fksize100k',
        'chr4_sv778_DEL_non_best8_fksize200k'}


def ols_coef(X, y):
    XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ (X.T @ y); resid = y - X @ beta
    n, p = X.shape; s2 = float(resid @ resid) / (n - p); se = np.sqrt(np.diag(XtXi) * s2)
    return beta[1], se[1], 2 * tdist.sf(abs(beta[1] / se[1]), n - p)


def irnt(x):
    m = np.isfinite(x); out = np.full(x.shape, np.nan)
    r = rankdata(x[m]); out[m] = norm.ppf((r - 0.5) / m.sum()); return out


tmap_loinc = {}
for src in (f"{RV}/step6_222/testable_trait_map.tsv",):
    for r in csv.DictReader(open(src), delimiter='\t'):
        tmap_loinc[r['ukb_trait']] = r.get('loinc', '')
on_target = set()
for pf in (f"{H}/pairs_222.tsv",):
    for r in csv.DictReader(open(pf), delimiter='\t'):
        lc = tmap_loinc.get(r['trait'], '')
        if lc:
            on_target.add((r['sv_locus'], lc))

base = ['age_at_recruitment_spline_1', 'age_at_recruitment_spline_2', 'age_at_recruitment_spline_3',
        'sex', *[f'gPC{i}' for i in range(1, 17)]]
cov = pd.read_csv(COV, sep='\t', dtype={'IID': str}, usecols=lambda c: c in ['IID', *base])
cov['sex'] = cov['sex'].map({'Male': 1.0, 'Female': 0.0})
for c in base:
    cov[c] = pd.to_numeric(cov[c], errors='coerce')
cov = cov.dropna(subset=base).reset_index(drop=True)
iids = cov['IID'].tolist(); idx = {i: k for k, i in enumerate(iids)}; N = len(iids)
C = cov[base].to_numpy(np.float64); del cov; gc.collect()
print(f"EUR covariate cohort N={N}", flush=True)

def header_svs(path):
    with gzip.open(path, 'rt') as f:
        return f.readline().rstrip('\n').split('\t')[1:]
svs1, svs2 = header_svs(R1), header_svs(R2)
sv_union = list(dict.fromkeys(svs1 + svs2))
svcol = {s: k for k, s in enumerate(sv_union)}
D = np.full((N, len(sv_union)), np.nan, np.float32)
for path, svs in ((R1, svs1), (R2, svs2)):
    take = [(j, svcol[s]) for j, s in enumerate(svs) if svcol[s] >= 0 and np.isnan(D[:, svcol[s]]).all()]
    with gzip.open(path, 'rt') as f:
        f.readline()
        for line in f:
            p = line.rstrip('\n').split('\t'); i = idx.get(p[0])
            if i is None:
                continue
            for j, col in take:
                D[i, col] = float(p[1 + j])
print(f"SV union = {len(sv_union)} SVs ({len(svs1)} round1 + {len(set(svs2)-set(svs1))} round2-only)", flush=True)

lv = {lc: np.full(N, np.nan, np.float32) for lc in MED}
for r in csv.DictReader(open(MEAS), delimiter='\t'):
    i = idx.get(r['IID'])
    if i is not None and r['loinc'] in lv:
        lv[r['loinc']][i] = float(r['median_value'])
lirnt = {lc: irnt(v.astype(float)) for lc, v in lv.items()}; del lv; gc.collect()

Cmask = np.isfinite(C).all(1)
res = []
for lc, tname in MED.items():
    y = lirnt[lc]
    for sv in sv_union:
        x = D[:, svcol[sv]].astype(np.float64)
        mask = np.isfinite(y) & np.isfinite(x) & Cmask
        n = int(mask.sum())
        ot = 1 if (sv, lc) in on_target else 0
        if n < 100:
            res.append([sv, tname, lc, n, '', '', '', ot]); continue
        X = np.column_stack([np.ones(n), x[mask], C[mask]]).astype(np.float64)
        keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9; X = X[:, keep]
        try:
            b, se, pv = ols_coef(X, y[mask].astype(np.float64))
            res.append([sv, tname, lc, n, f"{b:.5g}", f"{se:.5g}", f"{pv:.3g}", ot])
        except Exception as e:
            res.append([sv, tname, lc, n, '', '', f"err:{type(e).__name__}", ot])
        del X
print(f"ran {len(res)} SV x mediator tests", flush=True)

off = [(k, float(r[6])) for k, r in enumerate(res) if r[7] == 0 and r[6] not in ('', None) and not str(r[6]).startswith('err')]
q = {}
if off:
    ps = np.array([p for _, p in off]); order = np.argsort(ps); m = len(ps); prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]; prev = min(prev, ps[i] * m / (rank + 1)); q[off[i][0]] = prev

HDR = "sv_locus\tmediator\tloinc\tn\tbeta\tse\tp\ton_target\tq_bh_offtarget\n"
with open(f"{OUT}/mediator_scan_eur.tsv", 'w') as f:
    f.write(HDR)
    for k in sorted(range(len(res)), key=lambda k: (res[k][6] in ('', None), float(res[k][6]) if (res[k][6] not in ('', None) and not str(res[k][6]).startswith('err')) else 9)):
        r = res[k]
        f.write("\t".join(str(x) for x in r) + f"\t{q.get(k, float('nan')):.3g}\n")
print(f"wrote {OUT}/mediator_scan_eur.tsv", flush=True)

with open(f"{OUT}/mediator_scan_disease_hits_eur.tsv", 'w') as f:
    f.write(HDR)
    rows = [(k, res[k]) for k in range(len(res)) if res[k][0] in HITS]
    for k, r in sorted(rows, key=lambda kr: (kr[1][0], float(kr[1][6]) if (kr[1][6] not in ('', None) and not str(kr[1][6]).startswith('err')) else 9)):
        f.write("\t".join(str(x) for x in r) + f"\t{q.get(k, float('nan')):.3g}\n")
print(f"wrote {OUT}/mediator_scan_disease_hits_eur.tsv", flush=True)
