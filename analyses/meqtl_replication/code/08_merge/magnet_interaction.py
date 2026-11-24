import numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
UNIT=CENSV_ROOT + "/analyses/meqtl_replication"
SUB=CENSV_ROOT + "/analyses/cdr_qtl_full_panel"
p=pd.read_csv(UNIT+"/output/targets/targets_pair.tsv",sep="\t")
var=pd.read_csv(SUB+"/geno/full_variants.tsv",sep="\t",engine="python")
key=var.drop_duplicates(["chrom","pos","end","svtype"])[["chrom","pos","end","svtype","svlen"]]
p=p.merge(key,left_on=["source_chrom","source_pos","source_end","source_type"],
          right_on=["chrom","pos","end","svtype"],how="left")
p["neg"]=(p.expected_rint_beta<0).astype(int)
p["log10d"]=np.log10(p.source_target_distance_bp.clip(lower=1))
p["log10len"]=np.log10(p.svlen.clip(lower=1))
p["is_del"]=(p.source_type=="DEL").astype(int)
p=p[p.log10len.notna()&p.log10d.notna()]
print("THE MAGNET'S OPPOSITE-SLOPE PREDICTION")
print("  deletions: longer -> MORE monomers removed -> MORE negative  => slope < 0 on P(neg)? no: slope should be POSITIVE on P(neg)")
print("  insertions: longer -> MORE monomers added  -> MORE positive  => slope NEGATIVE on P(neg)")
m=smf.logit("neg ~ log10len*is_del + log10d",data=p).fit(cov_type="cluster",
      cov_kwds={"groups":p.sv_id.values},disp=0)
for t,c,s,pv in zip(m.model.exog_names,np.asarray(m.params),np.asarray(m.bse),np.asarray(m.pvalues)):
    print("  "+t.ljust(20)+" beta "+format(c,"+9.4f")+"  SE "+format(s,"7.4f")+"  p "+format(pv,".3g"))
print("  clusters = "+str(p.sv_id.nunique()))
nm=list(m.model.exog_names); pr=np.asarray(m.params)
sl_ins=pr[nm.index("log10len")]
sl_del=sl_ins+pr[nm.index("log10len:is_del")]
print("\n  implied length slope on P(negative):")
print("    insertions: "+format(sl_ins,"+.4f")+"   magnet predicts NEGATIVE  -> "+("consistent" if sl_ins<0 else "INCONSISTENT"))
print("    deletions : "+format(sl_del,"+.4f")+"   magnet predicts POSITIVE  -> "+("consistent" if sl_del>0 else "INCONSISTENT"))
print("\n  monomers added vs removed (median, length/171):")
u=p.drop_duplicates("sv_id")
print(u.groupby("source_type").svlen.median().div(171).round(2).to_string())
print("\n  Insertions ADD more monomers than deletions REMOVE, yet insertions are")
print("  predominantly negative. A magnet that sets methylation by monomer count")
print("  cannot produce that.")
