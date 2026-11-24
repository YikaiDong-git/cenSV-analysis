#!/usr/bin/env python3
import pandas as pd, numpy as np, os
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BD = CENSV_ROOT + ""
SV = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
CHRS = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
STRAT = "best8"
OUT = Path(BD) / "results/postqtl_v2b"
OUT.mkdir(parents=True, exist_ok=True)

def bh_correct(qtype):
    base = f"{BD}/results/{qtype}_overlap/decompose/cosine_nj/{STRAT}_{SV}"
    peer = 60 if qtype == "eqtl" else 15
    rows = []
    for c in CHRS:
        f = f"{base}/chr{c}/sexTrue_gpc5_peer{peer}/trans/chr{c}_trans_results.tsv"
        if os.path.exists(f):
            df = pd.read_csv(f, sep="\t")
            df["sv_chr"] = c
            if "is_cis" in df.columns:
                df = df[~df["is_cis"]]
            rows.append(df)
    if not rows:
        print(f"[{qtype}] no trans files found")
        return
    d = pd.concat(rows, ignore_index=True)
    pcol = next((k for k in ("pval_nominal","pvalue","p","pval") if k in d.columns), None)
    if pcol is None:
        print(f"[{qtype}] no pval column found in {list(d.columns)}")
        return
    d = d.dropna(subset=[pcol]).reset_index(drop=True)
    _, q, _, _ = multipletests(d[pcol].values, method="fdr_bh")
    d["qval_bh"] = q
    sig01 = (q < 0.1).sum()
    sig05 = (q < 0.05).sum()
    print(f"[{qtype}] trans pairs: {len(d):,} | q<0.1: {sig01:,} | q<0.05: {sig05:,}")
    d.to_csv(OUT / f"trans_{qtype}_bh.tsv.gz", sep="\t", index=False, compression="gzip")
    d[q < 0.1].to_csv(OUT / f"trans_{qtype}_bh_sig_q0.1.tsv", sep="\t", index=False)

bh_correct("eqtl")
bh_correct("sqtl")
print(f"\nwrote to {OUT}/")
