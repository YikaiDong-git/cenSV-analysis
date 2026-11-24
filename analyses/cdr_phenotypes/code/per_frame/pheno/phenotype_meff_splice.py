#!/usr/bin/env python
import pandas as pd, numpy as np, os, sys
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

_D = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_D, os.path.dirname(_D)]
import sys as _sys
_sys.path[:0] = [CENSV_ROOT + "/analyses/cdr_qtl_full_panel/code"]
from phenotype_family import PRIMARY, SUPPLEMENTARY, ALL_SCANNED, NPH, denominators, alpha

ROOT = CENSV_ROOT + ""
CAN = CENSV_ROOT + "/analyses/cdr_phenotypes/output/compat/compat_root_splice"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"

corr_t = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
deep_t = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
key = ["sample", "hap_local", "chm13_chrom"]
d = corr_t.merge(deep_t[key + ["cdr_pos_offset_deepest"]], on=key, how="left")
sub = d[ALL_SCANNED].apply(pd.to_numeric, errors="coerce").dropna()
print(f"complete-case hap-chromosomes across all {len(ALL_SCANNED)} phenotypes: {len(sub)}")

Z = {c: stats.norm.ppf((stats.rankdata(sub[c]) - 0.5) / len(sub)) for c in ALL_SCANNED}


def meff(cols):
    R = np.corrcoef(np.column_stack([Z[c] for c in cols]).T)
    ev = np.maximum(np.sort(np.linalg.eigvalsh(R))[::-1], 0)
    m = float(sum((abs(l) >= 1) + (abs(l) - np.floor(abs(l))) for l in ev))
    return min(len(cols), max(1.0, m)), ev


R = np.corrcoef(np.column_stack([Z[c] for c in ALL_SCANNED]).T)
Rdf = pd.DataFrame(R, index=ALL_SCANNED, columns=ALL_SCANNED)
print("\n=== phenotype correlation matrix (RINT, Pearson, complete cases) ===")
print(Rdf.round(3).to_string())

print("\n=== the two supplementary twins, against their primaries ===")
for tw, (pr, why) in SUPPLEMENTARY.items():
    print(f"  {tw:16s} vs {pr:18s}  r = {Rdf.loc[tw, pr]:+.4f}   ({why})")

m_all, ev_all = meff(ALL_SCANNED)
m_pri, ev_pri = meff(PRIMARY)
print(f"\neigenvalues, all {len(ALL_SCANNED)}: {np.round(ev_all, 3).tolist()}")
print(f"eigenvalues, primary {len(PRIMARY)}: {np.round(ev_pri, 3).tolist()}")
print(f"\nLi & Ji Meff, all {len(ALL_SCANNED)} scanned : {m_all:.2f}")
print(f"Li & Ji Meff, the {len(PRIMARY)} primary   : {m_pri:.2f}")
verdict = "CONSERVATIVE" if NPH >= m_pri else "TOO GENEROUS -- the split would need defending"
print(f"\nmultiplier actually used = {NPH}  ->  {verdict} "
      f"(x{NPH} vs Meff {m_pri:.2f} over the same {len(PRIMARY)} phenotypes)")

vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
rows = []
print(f"\n=== alpha under each multiplier ===")
print(f"{'framing':34s} {'n':>7s} {'x8 PRIMARY':>13s} {'x10 all':>12s} "
      f"{'xMeff(pri)':>12s} {'xMeff(all)':>12s}")
for nm, n in denominators(vt):
    a8, a10 = alpha(n), 0.05 / (n * len(ALL_SCANNED))
    amp, ama = 0.05 / (n * m_pri), 0.05 / (n * m_all)
    print(f"{nm:34s} {n:7d} {a8:13.3e} {a10:12.3e} {amp:12.3e} {ama:12.3e}")
    rows.append(dict(framing=nm, n=n, alpha_primary8=a8, alpha_x10_all_scanned=a10,
                     alpha_Meff_primary=amp, alpha_Meff_all=ama,
                     Meff_primary=m_pri, Meff_all=m_all, n_primary=NPH))
os.makedirs(f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/recompute/meff_splice", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/recompute/meff_splice/multiplier_summary.tsv", sep="\t", index=False,
                          float_format="%.6g")
Rdf.round(4).to_csv(f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/recompute/meff_splice/phenotype_correlation.tsv", sep="\t")
print(f"\nwrote ${CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/recompute/meff_splice/")
print(f"\nNote: 0.05/(n x 8) is numerically the legacy x8 convention, so every threshold the "
      f"predecessor published stays directly comparable.")
