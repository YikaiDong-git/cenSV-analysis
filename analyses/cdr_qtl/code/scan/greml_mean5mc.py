#!/usr/bin/env python3
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT = CENSV_ROOT + ""
DISC = f"{ROOT}/analyses/cdr_qtl"
ph = pd.read_csv(f"{DISC}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
fam = pd.read_csv(f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all.grm.id",
                  sep=r"\s+", header=None, names=["FID", "IID"])
ind = ph.groupby("sample").agg(mean_5mC=("mean_5mC", "mean"),
                               bc=("basecaller", "first")).reset_index()
ind["bc01"] = (ind["bc"] == "R1041").astype(int)
m = fam.merge(ind, left_on="IID", right_on="sample", how="inner")
m[["FID", "IID", "mean_5mC"]].to_csv(f"{DISC}/output/mean5mc.pheno", sep="\t",
                                     header=False, index=False)
m[["FID", "IID", "bc01"]].to_csv(f"{DISC}/output/bc.qcovar", sep="\t",
                                 header=False, index=False)
print(f"pheno n={len(m)}  mean_5mC mean={m.mean_5mC.mean():.2f} sd={m.mean_5mC.std():.2f}  "
      f"R1041={m.bc01.sum()}")
