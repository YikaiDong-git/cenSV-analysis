# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT=Path(CENSV_ROOT + "")
SUF="mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
BEST8=ROOT/f"results/imputation/decompose/cosine_nj/best8_{SUF}"
ENED=ROOT/f"results/elastic_net/decompose/cosine_nj/edge1000_{SUF}"
OUT=ROOT/"results/abstract_numbers"; OUT.mkdir(parents=True,exist_ok=True)
CHROMS=[f"chr{i}" for i in list(range(1,23))+["X"]]

b8=[]; me=[]
for c in CHROMS:
    f1=BEST8/c/f"{c}_best_strategy_selection.tsv"
    if f1.exists():
        d=pd.read_csv(f1,sep="\t"); d["chrom"]=c; b8.append(d[["chrom","sv_index","best_cv_r2"]])
    f2=ENED/c/"sv_metrics.tsv"
    if f2.exists():
        d=pd.read_csv(f2,sep="\t",low_memory=False); d["chrom"]=c
        me.append(d[["chrom","sv_index","sv_category","svfreq","svtype","svlen","n_nonzero_coef","original_sv_id"]])
b8=pd.concat(b8,ignore_index=True); me=pd.concat(me,ignore_index=True)
for x in (b8,me): x["sv_index"]=pd.to_numeric(x["sv_index"],errors="coerce").astype("Int64")
df=b8.merge(me,on=["chrom","sv_index"],how="inner")
for col in ["best_cv_r2","svfreq","svlen","n_nonzero_coef"]: df[col]=pd.to_numeric(df[col],errors="coerce")
df.to_csv(OUT/"persv_cat_cvr2_freq.tsv",sep="\t",index=False)

BINS=[0.0,0.05,0.10,0.20,0.30,1.01]
LAB=["<0.05","0.05-0.10","0.10-0.20","0.20-0.30","≥0.30"]
df["fbin"]=pd.cut(df.svfreq,BINS,right=False,labels=LAB)
rows=[]
for cls in ["monophyletic","non_decomposed"]:
    for i,lab in enumerate(LAB):
        s=df[(df.sv_category==cls)&(df.fbin==lab)]["best_cv_r2"].dropna()
        rows.append(dict(freq_bin=lab,lo=BINS[i],hi=BINS[i+1],
                         cls=("monophyletic" if cls=="monophyletic" else "unresolved"),
                         mean_cv_r2=round(float(s.mean()),4) if len(s) else np.nan,
                         median_cv_r2=round(float(s.median()),4) if len(s) else np.nan,n=int(len(s))))
cbin=pd.DataFrame(rows); cbin.to_csv(OUT/"cvr2_by_freqbin_class.tsv",sep="\t",index=False)

sub=df[df.sv_category.isin(["monophyletic","non_decomposed"])].dropna(subset=["best_cv_r2","svfreq","svlen"]).copy()
sub=sub[sub.svlen.abs()>0]
X=np.column_stack([np.ones(len(sub)),sub.svfreq,np.log10(sub.svlen.abs()),(sub.sv_category=="monophyletic").astype(float)])
y=sub.best_cv_r2.values; beta,_,_,_=np.linalg.lstsq(X,y,rcond=None)
resid=y-X@beta; n,k=X.shape; se=np.sqrt(np.diag(((resid@resid)/(n-k))*np.linalg.inv(X.T@X))); tt=beta/se
orig=df[df.sv_category=="original_decomposed"][["chrom","original_sv_id","best_cv_r2","svfreq"]].rename(columns={"best_cv_r2":"o_cv","svfreq":"o_f"})
mo=df[df.sv_category=="monophyletic"].groupby(["chrom","original_sv_id"]).agg(m_cv=("best_cv_r2","mean"),m_f=("svfreq","mean")).reset_index()
pr=orig.merge(mo,on=["chrom","original_sv_id"]).dropna()
with open(OUT/"resolution_freq_control.txt","w") as fh:
    fh.write("Resolution-vs-frequency control (best8, n=%d)\n"%len(df))
    fh.write("Freq-stratified mono vs unresolved (mean cv_r2):\n%s\n\n"%cbin.to_string(index=False))
    fh.write("OLS cv_r2 ~ svfreq + log10|svlen| + is_mono (mono+unresolved):\n")
    for nm,b_,t_ in zip(["intercept","svfreq","log10_svlen","is_monophyletic"],beta,tt):
        fh.write("  %-16s beta=%+.4f t=%+.1f\n"%(nm,b_,t_))
    w,pw=stats.wilcoxon(pr.m_cv,pr.o_cv)
    fh.write("\nPaired original->monophyletic (same parent, n=%d): orig cv=%.3f freq=%.4f ; mono cv=%.3f freq=%.4f ; component rarer in %.0f%% ; Wilcoxon p=%.1e\n"
             %(len(pr),pr.o_cv.mean(),pr.o_f.mean(),pr.m_cv.mean(),pr.m_f.mean(),100*(pr.m_f<pr.o_f).mean(),pw))
    for cov in ["svfreq","n_nonzero_coef"]:
        a=df.dropna(subset=[cov,"best_cv_r2"]); fh.write("corr(%s, cv_r2)=%+.3f\n"%(cov,np.corrcoef(a[cov],a.best_cv_r2)[0,1]))
print("WROTE:"); print(cbin.to_string(index=False))
print("\nOK best8 n=%d"%len(df))
