#!/usr/bin/env python3
"""Phase E: aggregate genome-wide cis result, simulate power at N=49, negative-control permutation."""
import glob, json, re, os
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P=CENSV_ROOT + ""
D=f"{P}/analyses/qtl_power"
W=f"{D}/output/qtl_c_ribo49/win10000000"
SUF="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
OUT=f"{D}/output/phase_e_n49"; os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(42)

bed=pd.read_csv(f"{D}/data/c_layers/ribo_n49.bed.gz",sep="\t")
pheno=bed.set_index("ID")[[c for c in bed.columns if c.startswith("NA")]]
cov=pd.read_csv(f"{D}/data/c_layers/ribo_n49.cov.tab.gz",sep="\t",index_col=0)
cov.loc["sex"]=cov.loc["sex"].map({"male":1,"female":0,"XY":1,"XX":0}).fillna(0)
cov=cov.apply(pd.to_numeric,errors="coerce").fillna(0)
samples=[s for s in pheno.columns if s in cov.columns]
pheno=pheno[samples]; cov=cov[samples]
C=np.column_stack([np.ones(len(samples)), cov.T.values.astype(float)])
n=len(samples); df=n-np.linalg.matrix_rank(C)-1
print(f"n={n}, covariate rank={np.linalg.matrix_rank(C)}, residual df={df}")

def resid(v):
    beta,_,_,_=np.linalg.lstsq(C,v,rcond=None); return v - C@beta
def partial_p(d,y):
    dr,yr=resid(d),resid(y)
    sd,sy=dr.std(),yr.std()
    if sd<1e-9 or sy<1e-9: return np.nan,0.0
    r=np.corrcoef(dr,yr)[0,1]; r=max(min(r,0.999999),-0.999999)
    t=r*np.sqrt(df/(1-r*r)); return 2*stats.t.sf(abs(t),df), r*r

perm=pd.concat([pd.read_csv(f,sep="\t") for f in sorted(glob.glob(f"{W}/*/chr*_cis_permutation.tsv"))],ignore_index=True)
perm=perm.dropna(subset=["pval_beta"]).copy()
perm["chrom"]=perm["variant_id"].map(lambda v: None)
perm["q_bh"]=multipletests(perm["pval_beta"].values,method="fdr_bh")[1]
perm.sort_values("pval_beta").to_csv(f"{OUT}/genomewide_cis_perm.tsv",sep="\t",index=False)
agg={"n_genes_tested":int(len(perm)),"min_pval_beta":float(perm.pval_beta.min()),
     "n_q05":int((perm.q_bh<=0.05).sum()),"n_q10":int((perm.q_bh<=0.10).sum()),
     "top_gene":str(perm.sort_values('pval_beta').iloc[0]['phenotype_id']),
     "top_variant":str(perm.sort_values('pval_beta').iloc[0]['variant_id']),
     "top_q":float(perm.sort_values('pval_beta').iloc[0]['q_bh'])}
alpha_bonf=0.05/len(perm)

def lead_dosage(chrom_file, sv_index):
    dos=pd.read_csv(chrom_file,sep="\t",index_col=0)
    if sv_index in dos.index:
        row=dos.loc[sv_index]
        return row.reindex(samples).astype(float).values
    return None
lead=[]
for f in sorted(glob.glob(f"{W}/*/chr*_cis_permutation.tsv")):
    chrom=re.search(r'(chr[0-9XY]+)_cis_permutation',f).group(1)
    t=pd.read_csv(f,sep="\t").dropna(subset=["pval_beta"])
    dfile=f"{D}/imputation/best8_{SUF}/{chrom}/{chrom}_imputed_dosages.tsv"
    for _,r in t.iterrows():
        m=re.match(r'sv_(\d+)_',str(r["variant_id"]))
        if m: lead.append((chrom,int(m.group(1)),r["phenotype_id"],dfile))
print(f"lead pairs: {len(lead)}")

idx=rng.choice(len(lead),size=min(200,len(lead)),replace=False)
R2grid=[0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]; M=400
dos_cache={}
power={r2:[] for r2 in R2grid}
mafs=[]
for j in idx:
    chrom,svi,gene,dfile=lead[j]
    if dfile not in dos_cache: dos_cache[dfile]=pd.read_csv(dfile,sep="\t",index_col=0)
    dos=dos_cache[dfile]
    if svi not in dos.index: continue
    d=dos.loc[svi].reindex(samples).astype(float).values
    dr=resid(d)
    if dr.std()<1e-9: continue
    dz=dr/dr.std()
    maf=min(np.mean(d)/2,1-np.mean(d)/2); mafs.append(maf)
    for r2 in R2grid:
        hits=0
        for _ in range(M):
            er=resid(rng.standard_normal(n))
            sd=er.std()
            if sd<1e-9: continue
            er=er/sd
            y=np.sqrt(r2)*dz+np.sqrt(1-r2)*er
            p,_=partial_p(d,y)
            if p<alpha_bonf: hits+=1
        power[r2].append(hits/M)
powsum={f"R2_{r2}":{"median_power":float(np.median(power[r2])),"mean_power":float(np.mean(power[r2]))} for r2 in R2grid}
md=[(r2,np.median(power[r2])) for r2 in R2grid]
detect=[r2 for r2,p in md if p>=0.8]
min_detect=min(detect) if detect else None

perm_samples=list(samples); rng.shuffle(perm_samples)
permmap=dict(zip(samples,perm_samples))
negp=[]
for chrom,svi,gene,dfile in lead:
    if gene not in pheno.index: continue
    if dfile not in dos_cache: dos_cache[dfile]=pd.read_csv(dfile,sep="\t",index_col=0)
    dos=dos_cache[dfile]
    if svi not in dos.index: continue
    d=dos.loc[svi].reindex(samples).astype(float).values
    y=pheno.loc[gene, [permmap[s] for s in samples]].astype(float).values
    p,_=partial_p(d,y)
    if not np.isnan(p): negp.append(p)
negp=np.array(negp)
ks=stats.kstest(negp,"uniform")
neg_q=multipletests(negp,method="fdr_bh")[1] if len(negp) else np.array([])
negctrl={"n":int(len(negp)),"ks_uniform_stat":float(ks.statistic),"ks_uniform_p":float(ks.pvalue),
         "n_q05_under_permutation":int((neg_q<=0.05).sum()) if len(neg_q) else 0,
         "frac_p_lt_0.05":float(np.mean(negp<0.05)) if len(negp) else None}

res={"aggregate":agg,"alpha_bonferroni":alpha_bonf,
     "power_sim":{"n_pairs_simulated":int(len(mafs)),"reps":M,"median_maf":float(np.median(mafs)),
                  "by_R2":powsum,"min_R2_detectable_at_80pct_median_power":min_detect},
     "negative_control":negctrl}
json.dump(res,open(f"{OUT}/phase_e_summary.json","w"),indent=2)
print(json.dumps(res,indent=2))

med=[float(np.median(power[r2])) for r2 in R2grid]
pd.DataFrame({"R2":R2grid,"median_power":med}).to_csv(
    f"{OUT}/power_curve.tsv",sep="\t",index=False)
pd.DataFrame({"perm_p":negp}).to_csv(f"{OUT}/negcontrol_perm_pvals.tsv",sep="\t",index=False)
print(f"wrote power_curve.tsv ({len(R2grid)} points) and "
      f"negcontrol_perm_pvals.tsv ({len(negp)} p-values)")

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
ax[0].plot([r*100 for r in R2grid],med,"o-",color="#2c7fb8")
ax[0].axhline(0.8,ls="--",c="grey",lw=1.0)
ax[0].text(5.5,0.82,"80% power",fontsize=8.5,color="grey")
ax[0].set_xlabel("True cis variance explained by the SV (%)")
ax[0].set_ylabel(r"Median power at $\alpha_{\mathrm{Bonferroni}}$ = %.1e"%alpha_bonf)
ax[0].set_title("A   Power is unreachable at this sample size",loc="left",fontweight="bold")
ax[0].set_ylim(0,1); ax[0].grid(alpha=0.25); ax[0].set_axisbelow(True)
ax[1].hist(negp,bins=20,color="#bdbdbd",edgecolor="k",linewidth=0.6)
ax[1].axhline(len(negp)/20.0,ls="--",c="#c0392b",lw=1.0)
ax[1].text(0.02,len(negp)/20.0*1.06,"uniform expectation",fontsize=8.5,color="#c0392b")
ax[1].set_xlabel("Permuted-label nominal $P$ (lead SV per gene)")
ax[1].set_ylabel("Genes")
ax[1].set_title(r"B   The permuted null is uniform (KS $P$ = %.2f)"%ks.pvalue,
                loc="left",fontweight="bold")
ax[1].grid(axis="y",alpha=0.25); ax[1].set_axisbelow(True)
fig.tight_layout(); fig.savefig(f"{OUT}/power_and_negcontrol.jpg",dpi=200,
                                bbox_inches="tight",pil_kwargs={"quality":93})
print("PHASE_E_DONE")
