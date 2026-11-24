#!/usr/bin/env python
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
OUT=CENSV_ROOT + "/analyses/cdr_qtl/canonical_geno"
PC=f"{OUT}/perchrom"
G=np.load(f"{OUT}/cdr6519_genotypes.npz"); G=G[list(G.keys())[0]].astype(float)
vt=pd.read_csv(f"{OUT}/cdr6519_variants.tsv",sep="\t")
samples=open(f"{OUT}/cdr6519_samples.txt").read().split()
assert G.shape==(len(vt),len(samples)), (G.shape,len(vt),len(samples))
nmiss_before=int((G<0).sum())
for i in range(G.shape[0]):
    row=G[i,:]; m=row<0
    if m.any():
        obs=row[~m]
        row[m]=obs.mean() if len(obs) else 0.0
        G[i,:]=row
print(f"imputed {nmiss_before} missing cells (locus-mean)")

for ch,sub in vt.groupby("chrom"):
    idx=sub.index.values
    d=f"{PC}/{ch}"; os.makedirs(d,exist_ok=True)
    gsub=G[idx,:]
    np.savez_compressed(f"{d}/{ch}_sv_genotypes.npz", genotypes=gsub)
    open(f"{d}/{ch}_sv_samples.txt","w").write("\n".join(samples)+"\n")
    sub.reset_index(drop=True).to_csv(f"{d}/{ch}_sv_variants.tsv",sep="\t",index=False)
print(f"wrote per-chrom files for {vt.chrom.nunique()} chroms under {PC}")
print(vt.chrom.value_counts().sort_index().to_string())
