import numpy as np, os, json, csv
BASE="/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds/step30_phewide_ext2"
DM=os.path.join(BASE,"disease_matrices"); DF=os.path.join(BASE,"disease_matrices_firth")

svs=[]; tiers=[]
for ln in open(os.path.join(DM,"disease_svorder.txt")):
    p=ln.rstrip("\n").split("\t")
    svs.append(p[0]); tiers.append(p[1] if len(p)>1 else "NA")
svs=np.array(svs); tiers=np.array(tiers)

rows=list(csv.DictReader(open(os.path.join(DM,"disease_meta.tsv")),delimiter="\t"))
hdr=list(rows[0].keys())
print("META COLS:",hdr)
phe=np.array([r["phecode"] for r in rows])
lab=np.array([r.get("label","") for r in rows])
tested=np.array([int(r["tested"]) for r in rows])
degcol=[c for c in hdr if "degen" in c.lower()]
print("degenerate-like cols:",degcol)
keep_col=(tested==1)
if degcol:
    dg=np.array([int(r[degcol[0]]) for r in rows]); keep_col&= (dg==0)
print("n_sv=%d n_phe=%d tested_cols=%d"%(len(svs),len(phe),keep_col.sum()))

P=np.load(os.path.join(DF,"disease_Pv_final.npy"),mmap_mode='r')
print("Pv_final shape",P.shape,P.dtype)
assert P.shape==(len(svs),len(phe))

colidx=np.where(keep_col)[0]

def bh(pvals):
    p=np.asarray(pvals,dtype=np.float64)
    n=p.size
    o=np.argsort(p,kind='mergesort')
    ps=p[o]
    q=ps*n/np.arange(1,n+1)
    q=np.minimum.accumulate(q[::-1])[::-1]
    q=np.minimum(q,1.0)
    out=np.empty(n,dtype=np.float64); out[o]=q
    return out

res={}
for tier in ["null_ext","signal"]:
    ridx=np.where(tiers==tier)[0]
    sub=np.array(P[ridx,:][:,colidx],dtype=np.float64)
    fin=np.isfinite(sub)&(sub>=0)&(sub<=1)
    nvalid=int(fin.sum())
    pv=sub[fin]
    q=bh(pv)
    Q=np.full(sub.shape,np.nan); Q[fin]=q
    res[tier]=dict(n_sv=len(ridx),n_cells=int(sub.size),n_valid=nvalid,
                   min_p=float(np.nanmin(pv)) if nvalid else None,
                   min_q=float(np.nanmin(q)) if nvalid else None,
                   n_q_lt_010=int((q<0.10).sum()),
                   n_q_lt_005=int((q<0.05).sum()))
    w=np.where(svs[ridx]=="chr8_sv30_DEL_mon_best8_midpoint1000")[0]
    if len(w):
        r=w[0]
        rowq=Q[r,:]; rowp=sub[r,:]
        j=int(np.nanargmin(rowp))
        res[tier]["chr8_sv30_best_p"]=float(rowp[j])
        res[tier]["chr8_sv30_best_q"]=float(rowq[j])
        res[tier]["chr8_sv30_best_phecode"]=phe[colidx][j]
        res[tier]["chr8_sv30_best_label"]=lab[colidx][j]
        res[tier]["chr8_sv30_n_valid_cells"]=int(np.isfinite(rowp).sum())
        res[tier]["chr8_sv30_min_q_row"]=float(np.nanmin(rowq))
        ordr=np.argsort(np.where(np.isfinite(rowp),rowp,2.0))[:5]
        res[tier]["chr8_sv30_top5"]=[(phe[colidx][k],lab[colidx][k],float(rowp[k]),float(rowq[k])) for k in ordr]
    else:
        gi=np.where(svs=="chr8_sv30_DEL_mon_best8_midpoint1000")[0][0]
        rowp=np.array(P[gi,colidx],dtype=np.float64)
        j=int(np.nanargmin(rowp)); pstar=float(rowp[j])
        allp=np.sort(pv)
        n=allp.size
        k=int(np.searchsorted(allp,pstar,side='right'))
        qraw=pstar*n/max(k,1)
        tail=allp[k:] if k<n else np.array([])
        cand=np.min(tail*n/np.arange(k+1,n+1)) if tail.size else 1.0
        res[tier]["chr8_sv30_best_p"]=pstar
        res[tier]["chr8_sv30_best_phecode"]=phe[colidx][j]
        res[tier]["chr8_sv30_best_label"]=lab[colidx][j]
        res[tier]["chr8_sv30_q_vs_this_denom"]=float(min(1.0,min(qraw,cand)))
        res[tier]["chr8_sv30_rank_in_tier"]=k
    del sub,fin,Q,pv,q
    import gc; gc.collect()

gi=np.where(svs=="chr8_sv30_DEL_mon_best8_midpoint1000")[0][0]
rowp=np.array(P[gi,colidx],dtype=np.float64)
j=int(np.nanargmin(rowp))
res["chr8_sv30_global"]=dict(row_index=int(gi),tier=str(tiers[gi]),
    best_p=float(rowp[j]),best_phecode=phe[colidx][j],best_label=lab[colidx][j],
    n_tested_cells_finite=int(np.isfinite(rowp).sum()))
print(json.dumps(res,indent=2,default=str))
json.dump(res,open("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds/step97_discovery_requests/logs/step97_bh_verify.json","w"),indent=2,default=str)
