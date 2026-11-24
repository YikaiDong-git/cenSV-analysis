"""Is the AFR deficit a property of AFRICAN CENTROMERES, or of the CATALOGUE?"""
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT=CENSV_ROOT + ""
WU=f"{ROOT}/analyses/cdr_qtl_full_panel"
G=np.load(f"{WU}/geno/full_genotypes.npz"); G=G[list(G.keys())[0]].astype(float); G[G<0]=np.nan
vt=pd.read_csv(f"{WU}/geno/full_variants.tsv",sep="\t")
samples=np.array(open(f"{WU}/geno/full_samples.txt").read().split())
rep=vt.reset_index().sort_values("index").drop_duplicates("geno_group",keep="first")
ri=rep["index"].values; Gu=G[ri,:]; vu=vt.iloc[ri].reset_index(drop=True)
ped=pd.read_csv(f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt",sep=r"\s+")
sc=[c for c in ped.columns if c.lower() in ("sampleid","sample")][0]
spc=[c for c in ped.columns if "super" in c.lower()][0]
POP=np.array([ped.set_index(sc)[spc].to_dict().get(s,"NA") for s in samples])
POPS=["AFR","AMR","EAS","EUR","SAS"]

print("=== 1. The frequency floor ===")
print(f"  min cohort af in the substrate : {vu.af.min():.4f}")
print(f"  = {vu.af.min()*376:.0f} alt alleles out of 376 -- a hard ascertainment floor.")
print(f"  variants with af < 0.05: {(vu.af<0.05).sum()} of {len(vu)} ({100*(vu.af<0.05).mean():.0f}%)")
print("  -> the rare tail, where African-specific variation concentrates, is ABSENT by construction.")

print("\n=== 2. Does AFR-specific frequency track cv_r2 (taggability) INSIDE the retained range? ===")
frq={p:np.nanmean(Gu[:,POP==p],axis=1)/2 for p in POPS}
nonafr=np.nanmean(np.column_stack([frq[p] for p in POPS if p!="AFR"]),axis=1)
delta=frq["AFR"]-nonafr
cv=vu.cv_r2.values
ok=np.isfinite(cv)&np.isfinite(delta)
r,p_=stats.spearmanr(delta[ok],cv[ok])
print(f"  spearman(AFR-enrichment, cv_r2) = {r:+.3f}   p = {p_:.2g}   n = {ok.sum()}")
q=pd.qcut(pd.Series(cv[ok]),5,labels=["cv_r2 Q1 (worst tagged)","Q2","Q3","Q4","Q5 (best tagged)"])
t=pd.DataFrame({"q":q.values,"d":delta[ok]}).groupby("q",observed=True).d.agg(["mean","median","size"])
print(t.round(4).to_string())
print("  -> a NEGATIVE trend means AFR-enriched SVs are systematically harder to tag, so the")
print("     cv_r2 >= 0.1 substrate filter removes them preferentially.")

print("\n=== 3. AFR-enriched vs non-AFR-enriched variants: how many of each survive into the substrate? ===")
print(f"  variants commoner in AFR than in the non-AFR mean : {(delta>0).sum()} "
      f"({100*(delta>0).mean():.1f}%)")
print(f"  variants commoner outside AFR                     : {(delta<0).sum()} "
      f"({100*(delta<0).mean():.1f}%)")
print("  (a 50/50 split is what an ancestry-neutral catalogue would give)")
for lab,m in [("AFR-enriched",delta>0),("non-AFR-enriched",delta<0)]:
    print(f"    {lab:18s} median cv_r2 = {np.nanmedian(cv[m]):.3f}   median af = {np.nanmedian(vu.af.values[m]):.3f}")
