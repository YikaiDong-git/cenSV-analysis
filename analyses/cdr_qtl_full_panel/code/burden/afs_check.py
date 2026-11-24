"""Why does AFR carry FEWER alt alleles but MORE private variants?"""
import numpy as np, pandas as pd
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

af=np.nanmean(Gu,axis=1)/2.0
bins=[(0,0.01,"<1%"),(0.01,0.05,"1-5%"),(0.05,0.20,"5-20%"),(0.20,0.50,"20-50%"),(0.50,1.01,">50%")]
print("Mean alt-allele burden per individual, by COHORT allele-frequency bin:\n")
print(f"{'AF bin':10s} {'n_groups':>9s} " + "".join(f"{p:>9s}" for p in POPS) + f"{'AFR/EUR':>9s}")
tot={p:0.0 for p in POPS}
for lo,hi,lab in bins:
    m=(af>=lo)&(af<hi)
    row=f"{lab:10s} {int(m.sum()):9d} "
    vals={}
    for p in POPS:
        v=float(np.nanmean(np.nansum(Gu[m][:,POP==p],axis=0))); vals[p]=v; tot[p]+=v
        row+=f"{v:9.1f}"
    row+=f"{vals['AFR']/max(vals['EUR'],1e-9):9.2f}"
    print(row)
print(f"{'TOTAL':10s} {len(af):9d} "+"".join(f"{tot[p]:9.1f}" for p in POPS)
      +f"{tot['AFR']/tot['EUR']:9.2f}")

print("\nDistinct genotype-groups CARRIED (>=1 alt allele), equal-n bootstrap already in F7 panel B.")
print("Here: number of groups where each superpop has the HIGHEST allele frequency:")
best={p:0 for p in POPS}
for j in range(Gu.shape[0]):
    fr={p:np.nanmean(Gu[j,POP==p])/2 for p in POPS}
    if max(fr.values())>0: best[max(fr,key=fr.get)]+=1
print("  "+"  ".join(f"{p}={best[p]}" for p in POPS))

print("\nMean cohort AF of the variants each superpop carries (dosage-weighted):")
for p in POPS:
    b=np.nansum(Gu[:,POP==p],axis=1)
    print(f"  {p}: {np.average(af,weights=np.nan_to_num(b)):.4f}")
print("\n-> if AFR's deficit sits in the COMMON bins and its excess in the RARE bins, the pattern is")
print("   the out-of-Africa frequency spectrum plus a taggability-ascertained substrate, NOT")
print("   'AFR has less centromeric structural variation'.")
