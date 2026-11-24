import os, numpy as np, pandas as pd, pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
P=CENSV_ROOT + ""
U=P+"/analyses/meqtl_replication"
OUT=U+"/output/headline_svs"
T=pd.read_csv(OUT+"/null/headline_targets_raw.tsv",sep="\t")
assert T.target_class.notna().all() and set(T.target_class)=={"discovery","matched_null"}
fa=pysam.FastaFile(P+"/resources/chm13v2.0.fa.gz")
b=[];c=[];dn=[]
for t in T.itertuples(index=False):
    p=int(t.pos); base=fa.fetch(t.chrom,p,p+1).upper(); b.append(base)
    cc=p if base=="C" else (p-1 if base=="G" else -1); c.append(cc)
    dn.append(fa.fetch(t.chrom,cc,cc+2).upper() if cc>=0 else "")
fa.close()
T["chm13_base"]=b; T["plus_c_pos"]=c; T["chm13_dinuc"]=dn
T["anchor_ok"]=T.chm13_base.isin(["C","G"])&(T.chm13_dinuc=="CG")
print("anchor resolved "+str(int(T.anchor_ok.sum()))+"/"+str(len(T)))
print(pd.crosstab(T.chm13_base,T.plus_c_pos-T.pos).to_string())
T=T[T.anchor_ok].copy()
T["cpg_key"]=T.chrom.astype(str)+":"+T.plus_c_pos.astype(str)
disc=set(T[T.target_class=="discovery"].cpg_key)
leak=T[(T.target_class=="matched_null")&T.cpg_key.isin(disc)]
print("strand leaks removed: "+str(len(leak)))
T=T[~((T.target_class=="matched_null")&T.cpg_key.isin(disc))]
T["cls_rank"]=(T.target_class!="discovery").astype(int)
T=T.sort_values(["cpg_key","cls_rank"])
rep=T.groupby("cpg_key",as_index=False).first()
rep["member_target_ids"]=rep.cpg_key; rep["n_members"]=1
rep["target_id"]="h"+pd.Series(np.arange(len(rep))).astype(str).str.zfill(6)
print("unique CpGs to lift: "+str(len(rep))+"  "+str(dict(rep.target_class.value_counts())))
o=rep[["cpg_key","target_id","target_class","feature_id","chrom","plus_c_pos",
       "chm13_base","member_target_ids","n_members","sv_name"]]
tmp=OUT+"/null/headline_targets_lift.tsv.tmp"
o.to_csv(tmp,sep="\t",index=False); os.replace(tmp,OUT+"/null/headline_targets_lift.tsv")
print("wrote headline_targets_lift.tsv")
