#!/usr/bin/env python
"""Is the cis-mQTL null CALIBRATED, or is the tail excess an artefact?"""
import sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
WU = f"{P1}/analyses/cdr_qtl_full_panel"
OUT = f"{WU}/output/mqtl"
SUF = sys.argv[1] if len(sys.argv) > 1 else "_bam"
RDIR = "bam" if SUF else "chain"

CHI2_MED = stats.chi2.ppf(0.5, 1)


def lam(t):
    t = np.asarray(t, float)
    t = t[np.isfinite(t)]
    return float(np.median(t ** 2) / CHI2_MED) if len(t) else np.nan


rows = []
for W in (10, 50, 100):
    p = f"{OUT}/scans/{RDIR}/scan_{W}kb{SUF}.tsv"
    try:
        S = pd.read_csv(p, sep="\t")
    except FileNotFoundError:
        print(f"[missing] {p}")
        continue
    S["phys"] = (S.chrom.astype(str) + ":" + S.pos.astype(str) + ":" +
                 S.svtype.astype(str) + ":" + S.svlen.astype(str))
    U = S.drop_duplicates("phys")
    for lab, D in (("all rows", S), ("distinct SVs", U)):
        r = dict(window_kb=W, set=lab, n=len(D),
                 lambda_meth=lam(D.t_meth), lambda_cov=lam(D.t_cov))
        for thr in (0.05, 0.01, 0.005):
            obs = int((D.p_meth <= thr).sum())
            r[f"meth_p<={thr}"] = f"{obs} vs {thr*len(D):.0f} ({obs/(thr*len(D)):.2f}x)"
            oc = int((D.p_cov <= thr).sum())
            r[f"cov_p<={thr}"] = f"{oc} vs {thr*len(D):.0f} ({oc/(thr*len(D)):.2f}x)"
        rows.append(r)

R = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(R[["window_kb", "set", "n", "lambda_meth", "lambda_cov"]].to_string(
    index=False, float_format=lambda x: f"{x:.3f}"))
print()
print(R[["window_kb", "set", "meth_p<=0.05", "meth_p<=0.01", "meth_p<=0.005"]].to_string(index=False))
print()
print(R[["window_kb", "set", "cov_p<=0.05", "cov_p<=0.01", "cov_p<=0.005"]].to_string(index=False))
R.to_csv(f"{OUT}/scans/{RDIR}/calibration{SUF}.tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {OUT}/scans/calibration{SUF}.tsv")
print("\nREADING GUIDE: lambda_meth ~ 1 with a tail excess = candidate polygenic signal worth the")
print("escalation. lambda_meth >> 1, or a comparable excess in the COVERAGE control, = the excess is")
print("technical and the tail count must not be read as biology.")
