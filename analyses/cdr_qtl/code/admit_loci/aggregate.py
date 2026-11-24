#!/usr/bin/env python
import glob, os, numpy as np, pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
DISC=CENSV_ROOT + "/analyses/cdr_qtl"
PERCHROM=f"{DISC}/results/canonical_scan_1e5/perchrom"
OUT=f"{DISC}/results/canonical_admit"
os.makedirs(OUT,exist_ok=True)
N_PERM=100000
PHENOS_SCANNED=["cdr_pos_offset","cdr_pos_offset_deepest"]
N_ADMIT=6519

def bh(p):
    p=np.asarray(p); n=len(p); o=np.argsort(p); r=np.empty(n)
    r[o]=(p[o]*n)/(np.arange(1,n+1))
    q=np.minimum.accumulate(r[o][::-1])[::-1]; out=np.empty(n); out[o]=np.minimum(q,1)
    return out

files=sorted(glob.glob(f"{PERCHROM}/*_perm.tsv"))
print(f"[load] {len(files)} perchrom files")
frames=[]
for f in files:
    d=pd.read_csv(f,sep="\t"); d["stratum"]=Path(f).stem.replace("_perm",""); frames.append(d)
alld=pd.concat(frames,ignore_index=True)
print(f"[load] strata={alld.stratum.nunique()} phenotypes={sorted(alld.phenotype.unique())} tests={len(alld)}")

alld["q_within"]=np.nan
for s,g in alld.groupby("stratum"):
    alld.loc[g.index,"q_within"]=bh(g["p_emp_perm"].values)
alld["is_floor"]=alld["count_perm_ge_obs"]==0
n_floor=int(alld.is_floor.sum())
print(f"\n=== floor hits @1e5 (p={1/(N_PERM+1):.2e}): {n_floor} ===")
print(alld[alld.is_floor].groupby("phenotype").size().to_string())

floor_cells=alld[alld.is_floor][["chrom","phenotype"]].drop_duplicates().sort_values(["phenotype","chrom"])
print(f"\n=== floor-saturated (chrom,phenotype) cells to escalate to 1e7: {len(floor_cells)} ===")
print(floor_cells.to_string(index=False))
floor_cells.to_csv(f"{OUT}/canonical_floor_cells.txt",sep="\t",index=False,header=False)

n_actual=len(alld)
frame_b=N_ADMIT*8
bonf_actual=0.05/n_actual
bonf_claim=0.05/frame_b
floor_p=1/(N_PERM+1); floor_p7=1/(10_000_000+1)
print(f"\n=== multiple-testing frames ===")
print(f"  (a) actual scan-row burden      n={n_actual:>7d}  Bonferroni={bonf_actual:.2e}")
print(f"  (b) canonical claim 6519x8       n={frame_b:>7d}  Bonferroni={bonf_claim:.2e}")
print(f"  1e5 floor p={floor_p:.2e} (ABOVE both -> unresolved); 1e7 floor p={floor_p7:.2e}")
print(f"  => 1e7 floor {'CLEARS' if floor_p7<bonf_claim else 'does NOT clear'} the 6519x8 claim frame")

cols=["stratum","chrom","pos","svtype","svlen","af","phenotype","t_obs","count_perm_ge_obs","p_emp_perm","q_within","n_indiv"]
cols=[c for c in cols if c in alld.columns]
sig=alld[alld.q_within<0.10].sort_values("p_emp_perm")
print(f"\n=== top canonical hits (q_within<0.10: {len(sig)}) ===")
print(sig.head(30)[cols].to_string(index=False))

print("\n=== chr8 pericentromeric cluster (both phenotypes) ===")
c8=alld[(alld.chrom=="chr8")&(alld.pos.between(44240000,44470000))&alld.is_floor].sort_values(["phenotype","p_emp_perm"])
print(c8[cols].to_string(index=False))

alld[alld.is_floor].sort_values("p_emp_perm").to_csv(f"{OUT}/canonical_floor_hits.tsv",sep="\t",index=False)
alld.to_csv(f"{OUT}/canonical_all_tests.tsv.gz",sep="\t",index=False,compression="gzip")
print(f"\n[written] {OUT}/canonical_floor_hits.tsv ({n_floor}), canonical_all_tests.tsv.gz ({len(alld)}), canonical_floor_cells.txt ({len(floor_cells)})")
