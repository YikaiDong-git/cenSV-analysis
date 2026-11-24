#!/usr/bin/env python3
"""Enlarged-run stats: dedupe synonym traits, BH-FDR, tier-wise Bonferroni, Meff, sign concordance."""
import csv, os, numpy as np, pandas as pd
BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
R6 = f"{BASE}/aou_replication/results_vds/step6_222"; R8 = f"{BASE}/aou_replication/results_vds/step8_222"
OUT = R8
tmap = {r['ukb_trait']: (r['type'], r['loinc']) for r in csv.DictReader(open(f"{R6}/testable_trait_map.tsv"), delimiter='\t')}

def loincof(trait):
    return tmap[trait][1] if trait in tmap else trait

eur = list(csv.DictReader(open(f"{R8}/results_222_eur.tsv"), delimiter='\t'))
allm = {(r['sv_locus'], r['trait']): r for r in csv.DictReader(open(f"{R8}/results_222_all.tsv"), delimiter='\t')}
def pf(x):
    try: return float(x)
    except: return None

for r in eur:
    r['loinc'] = loincof(r['trait'])
seen = {}
for r in eur:
    if pf(r['p']) is None:
        continue
    k = (r['sv_locus'], r['loinc'])
    if k not in seen or pf(r['p']) < pf(seen[k]['p']):
        seen[k] = r
uniq = list(seen.values())
print(f"pairs run: {len(eur)}  ->  unique (SV,LOINC) tests: {len(uniq)}  (dropped {len(eur)-len(uniq)} synonym duplicates)")

ps = np.array([pf(r['p']) for r in uniq]); m = len(ps); order = np.argsort(ps); q = np.empty(m); prev = 1.0
for k in range(m - 1, -1, -1):
    i = order[k]; prev = min(prev, ps[i] * m / (k + 1)); q[i] = prev
for r, qq in zip(uniq, q):
    r['q_bh'] = qq

loincs = sorted({r['loinc'] for r in uniq if tmap.get(r['trait'], ('quant',))[0] == 'quant'})
idxrow = {}; rows = []
for rr in csv.DictReader(open(f"{R6}/measurement_medians_222.tsv"), delimiter='\t'):
    if rr['loinc'] in loincs:
        if rr['IID'] not in idxrow:
            if len(idxrow) >= 120000:
                continue
            idxrow[rr['IID']] = len(rows); rows.append({})
        rows[idxrow[rr['IID']]][rr['loinc']] = float(rr['median_value'])
mat = pd.DataFrame(rows, columns=loincs).astype(float)
corr = mat.corr(min_periods=1000).to_numpy()
corr = np.nan_to_num(corr, nan=0.0); np.fill_diagonal(corr, 1.0)
corr = (corr + corr.T) / 2.0
M = corr.shape[0]
ev = np.clip(np.linalg.eigvalsh(corr), 0.0, None)
meff_liji = float(np.sum((ev >= 1).astype(float) + (ev - np.floor(ev))))
meff_nyholt = 1 + (M - 1) * (1 - np.var(ev, ddof=1) / M)
print(f"Meff over {M} AoU quant LOINCs: Li-Ji={meff_liji:.1f}, Nyholt-Cheverud={meff_nyholt:.1f} "
      f"(both ~M -> traits largely independent; BH-FDR is the operative control, Meff is only supporting)")

with open(f"{OUT}/results_222_unique_annotated.tsv", 'w') as f:
    f.write("sv_locus\ttrait\tloinc\ttier\tcv_r2\tdiscovery_sign\tn\tn_cases\tbeta\tse\tp\tq_bh\taou_sign\tconcordance\tp_all\n")
    for r in sorted(uniq, key=lambda r: pf(r['p'])):
        a = allm[(r['sv_locus'], r['trait'])]; aou = '+' if pf(r['beta']) > 0 else '-'
        f.write(f"{r['sv_locus']}\t{r['trait']}\t{r['loinc']}\t{r['tier']}\t{r['cv_r2']}\t{r['discovery_sign']}\t{r['n']}\t{r['n_cases']}\t{r['beta']}\t{r['se']}\t{r['p']}\t{r['q_bh']:.3g}\t{aou}\t{r['concordance']}\t{a['p']}\n")

from collections import Counter
tn = Counter(r['tier'] for r in uniq)
print(f"\nunique tests by tier: {dict(tn)}")
print("\n=== BH q<0.05 hits (unique) ===")
for r in sorted([r for r in uniq if r['q_bh'] < 0.05], key=lambda r: pf(r['p'])):
    a = allm[(r['sv_locus'], r['trait'])]
    print(f"  [{r['tier'][:4]}] {r['sv_locus']} x {r['trait']} (LOINC {r['loinc']}): p={pf(r['p']):.2g} q={r['q_bh']:.2g} beta={r['beta']} cv_r2={r['cv_r2']} disc={r['discovery_sign']} concord={r['concordance']} p_all={pf(a['p']):.2g}")

print("\n=== chr11 MEGA_PLEIOTROPIC (agnostic): breadth + directional concordance ===")
for sv in ['chr11_sv253_INS_ori_best8_midpoint500', 'chr11_sv252_INS_mon_best8_midpoint1000']:
    rr = [r for r in uniq if r['sv_locus'] == sv]
    nom = [r for r in rr if pf(r['p']) < 0.05]; qsig = [r for r in rr if r['q_bh'] < 0.05]
    conc = sum(1 for r in nom if r['concordance'] == 'concordant')
    print(f"  {sv}: {len(rr)} unique traits | nominal<0.05: {len(nom)} | BH q<0.05: {len(qsig)} | discovery-concordant among nominal: {conc}/{len(nom)}")
