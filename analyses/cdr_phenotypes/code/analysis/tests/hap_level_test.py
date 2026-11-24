import numpy as np, pandas as pd, os
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
P=CENSV_ROOT + ""
CAN=P+"/analyses/cdr_qtl/canonical_geno/decompose"
WU=P+"/analyses/cdr_qtl_full_panel"
RD=P+"/analyses/cdr_phenotypes"

hg=pd.read_csv(CAN+"/chr8_all_haplogroups.tsv",sep="\t")
hg=hg.reset_index(drop=True)
idx2hap={i:(r.sample_id,r.haplotype) for i,r in hg.iterrows()}
print("haplotype index table: %d entries (%d individuals)"%(len(hg),hg.sample_id.nunique()))

dec=pd.read_csv(CAN+"/chr8_all_decomposition.tsv",sep="\t")
sv=pd.read_csv(CAN+"/chr8_all_sv_variants.tsv",sep="\t")
lead=sv[(sv.pos==44300221)&(sv.svtype=="DEL")]
print("\nlead rows in chr8_all_sv_variants.tsv: %d"%len(lead))
print(lead[["pos","svtype","svlen","af","mapped_haplogroup","sv_category","r2_score"]].to_string(index=False))

G=np.load(WU+"/geno/full_genotypes.npz"); G=G[list(G.keys())[0]].astype(float)
vt=pd.read_csv(WU+"/geno/full_variants.tsv",sep="\t")
samp=np.array(open(WU+"/geno/full_samples.txt").read().split())
LEAD="chr8_sv30_DEL_mon_best8_midpoint1000"
row=vt.index[vt.sv_locus==LEAD][0]
dip_dos=pd.Series(G[row,:],index=samp).replace(-1,np.nan)
print("\nscan lead %s: dosage counts %s"%(LEAD,dip_dos.value_counts().to_dict()))

best=None
for r in dec.itertuples():
    ci=str(r.carrier_indices)
    if ci in ("nan",""): continue
    idxs=[int(x) for x in ci.split(",") if x.strip().isdigit()]
    cnt={}
    for i in idxs:
        s,h=idx2hap[i]; cnt[s]=cnt.get(s,0)+1
    vec=pd.Series({s:cnt.get(s,0) for s in samp})
    ok=(vec.reindex(samp).values==dip_dos.reindex(samp).values)
    agree=np.nanmean(ok[~np.isnan(dip_dos.reindex(samp).values)])
    if best is None or agree>best[0]:
        best=(agree,r,idxs)
agree,r,idxs=best
print("\nBEST-MATCHING decomposition row: branch=%s n_carriers=%s  diploid agreement with the scan genotype = %.4f"
      %(getattr(r,'branch_id',''),r.n_carriers,agree))
cnt={}
for i in idxs:
    s_,h_=idx2hap[i]; cnt[s_]=cnt.get(s_,0)+1
recon=pd.Series({s_:cnt.get(s_,0) for s_ in samp})
concord=recon.reindex(samp).values==dip_dos.reindex(samp).values
SAFE=set(np.array(samp)[concord & ~np.isnan(dip_dos.reindex(samp).values)])
print("  concordant individuals kept: %d of %d (dropped %d where B1 != scan lead)"
      %(len(SAFE),len(samp),len(samp)-len(SAFE)))

carrier=set(idxs)
rows=[]
for i,(s,h) in idx2hap.items():
    rows.append(dict(sample=s,hap_local=h,carr=1 if i in carrier else 0))
H=pd.DataFrame(rows)

ph=pd.read_csv(RD+"/output/phenotypes/per_hap_chrom_phenotypes_frames.tsv",sep="\t",engine="python")
ph=ph[(ph.chm13_chrom=="chr8")&ph.cdr_pos_offset_splice.notna()]
m=ph.merge(H,on=["sample","hap_local"],how="inner")
m=m[m["sample"].isin(SAFE)]
print("\npaired rows: %d haplotypes, %d individuals"%(len(m),m["sample"].nunique()))

w=m.pivot_table(index="sample",columns="hap_local",values="cdr_pos_offset_splice")
c=m.pivot_table(index="sample",columns="hap_local",values="carr")
both=w.dropna().index.intersection(c.dropna().index)
w,c=w.loc[both],c.loc[both]
het=c[(c.hap1+c.hap2)==1].index
print("\n=== WITHIN-INDIVIDUAL PAIRED TEST (heterozygotes only) ===")
print("informative individuals (exactly one carrier haplotype): %d"%len(het))
carr_pos=[]; non_pos=[]
for s in het:
    if c.loc[s,"hap1"]==1: carr_pos.append(w.loc[s,"hap1"]); non_pos.append(w.loc[s,"hap2"])
    else:                  carr_pos.append(w.loc[s,"hap2"]); non_pos.append(w.loc[s,"hap1"])
carr_pos=np.array(carr_pos); non_pos=np.array(non_pos); d=carr_pos-non_pos
t,p=stats.ttest_rel(carr_pos,non_pos)
ws,wp=stats.wilcoxon(carr_pos,non_pos)
print("  carrier hap   mean dip position: %+.4f"%carr_pos.mean())
print("  non-carrier   mean dip position: %+.4f"%non_pos.mean())
print("  paired difference (carrier - non): mean %+.4f, median %+.4f, SD %.4f"%(d.mean(),np.median(d),d.std(ddof=1)))
print("  paired t = %+.3f   P = %.3g"%(t,p))
print("  Wilcoxon  P = %.3g"%wp)
print("  sign: %d of %d individuals have the dip closer to the array START on the carrier haplotype"
      %(int((d<0).sum()),len(d)))
np.save(RD+"/output/phenotypes/hap_paired_diff.npy",d)
print("\n[saved] output/hap_paired_diff.npy")
