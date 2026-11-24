import csv, gzip, numpy as np, sys
sys.path.insert(0,"/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/code/phewide")
import step30g_numpy_disease_scoretest as S
AR="/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
D=f"{AR}/results_vds/step30_phewide"; MAT=f"{D}/disease_matrices"
OUT=f"{AR}/send_back_to_discovery/aou_persv_summary.tsv"

Pv=np.load(f"{MAT}/disease_Pv.npy"); Be=np.load(f"{MAT}/disease_Be.npy")
meta=list(csv.DictReader(open(f"{MAT}/disease_meta.tsv"),delimiter="\t"))
phe=[m["phecode"] for m in meta]; lab=[m["label"] for m in meta]
cases=np.array([int(m["cases"]) for m in meta]); 
tested=np.array([m["tested"]=="1" for m in meta]); degen=np.array([m.get("degenerate","0")=="1" for m in meta])
svrows=[l.rstrip("\n").split("\t") for l in open(f"{MAT}/disease_svorder.txt") if l.strip()]
sv=[r[0] for r in svrows]; tiers=np.array([r[1] if len(r)>1 else "" for r in svrows]); M,J=Pv.shape
sv_ix={s:i for i,s in enumerate(sv)}
refit={}
for r in csv.DictReader(open(f"{D}/phewide_disease_full_firthrefit.tsv"),delimiter="\t"):
    pf=r["p_firth"]; bf=r["beta_firth"]
    refit[(r["sv_locus"],r["phecode"])]=(float(pf) if pf not in("","nan") else np.nan,
                                          float(bf) if bf not in("","nan") else np.nan)
q=np.full((M,J),np.nan,np.float32)
validcol=tested&~degen
for tn in ("signal","null_ext"):
    rows=np.where(tiers==tn)[0]
    sub=Pv[np.ix_(rows,np.where(validcol)[0])]
    flat=sub.ravel(); fin=np.isfinite(flat)
    qf=S.bh(np.where(fin,flat,np.nan))
    qsub=np.full(sub.shape,np.nan,np.float32); qsub.ravel()[:]=qf
    q[np.ix_(rows,np.where(validcol)[0])]=qsub
dz_nbh=np.nansum(q<0.05,axis=1).astype(int)
Pvv=Pv.copy(); Pvv[:,~validcol]=np.nan
dz=[]
for i in range(M):
    row=Pvv[i]; 
    if not np.isfinite(row).any(): dz.append((np.nan,"",np.nan,0,0)); continue
    j=int(np.nanargmin(row)); sp=float(row[j]); bscore=float(Be[i,j])
    fp,fb=refit.get((sv[i],phe[j]),(np.nan,np.nan))
    mp=fp if np.isfinite(fp) else sp; lo=fb if np.isfinite(fb) else bscore
    ntest=int(np.isfinite(row).sum())
    dz.append((mp,phe[j],lo,int(cases[j]),ntest))

qt_minp={}; qt_arg={}; qt_beta={}; qt_ntest={}
allp=[]; allsv=[]; alltier=[]
tmap={"signal":0,"null_ext":1}
with gzip.open(f"{D}/phewide_quant_full_all.tsv.gz","rt") as f:
    rd=csv.reader(f,delimiter="\t"); next(rd)
    for s,ti,loinc,n,beta,p in rd:
        try: pv=float(p)
        except: continue
        if not np.isfinite(pv) or pv<=0: continue
        i=sv_ix.get(s)
        if i is None: continue
        if (s not in qt_minp) or pv<qt_minp[s]:
            qt_minp[s]=pv; qt_arg[s]=loinc; qt_beta[s]=float(beta)
        qt_ntest[s]=qt_ntest.get(s,0)+1
        allp.append(pv); allsv.append(i); alltier.append(tmap.get(ti,0))
allp=np.array(allp); allsv=np.array(allsv); alltier=np.array(alltier)
qt_nbh=np.zeros(M,int)
for t in (0,1):
    m=alltier==t
    qf=S.bh(allp[m]); sig=allsv[m][qf<0.05]
    for i in sig: qt_nbh[i]+=1
print("quant rows used:",len(allp),flush=True)

with open(OUT,"w") as f:
    f.write("sv_locus\ttier\tdz_min_p\tdz_argmin_phecode\tdz_argmin_logOR\tdz_argmin_cases\tdz_n_bh_sig\t"
            "qt_min_p\tqt_argmin_loinc\tqt_argmin_beta\tqt_n_bh_sig\tn_dz_tested\tn_qt_tested\n")
    for i,s in enumerate(sv):
        mp,ph,lo,ca,ndz=dz[i]
        f.write(f"{s}\t{tiers[i]}\t{mp:.4g}\t{ph}\t{lo:.4g}\t{ca}\t{dz_nbh[i]}\t"
                f"{qt_minp.get(s,float('nan')):.4g}\t{qt_arg.get(s,'')}\t{qt_beta.get(s,float('nan')):.4g}\t"
                f"{qt_nbh[i]}\t{ndz}\t{qt_ntest.get(s,0)}\n")
print("wrote",OUT,flush=True)

print("\n=== CHECK 1: chr11 hit-SVs across COAGULATION/thrombosis phecodes vs the genetic-testing phecode ===",flush=True)
coag=[j for j in range(J) if any(k in lab[j].lower() for k in ("thromb","embol","coagul","venous","prothrombin","hypercoag"))]
for s in ("chr11_sv200_DEL_non_best6_alt_edge2000","chr11_sv400_INS_non_best8_fksize100k","chr11_sv148_INS_non_best8_fksize500k"):
    i=sv_ix[s]; print(f" {s}",flush=True)
    for j in coag:
        p=Pv[i,j]
        if np.isfinite(p): print(f"    {phe[j]:11s} {lab[j][:34]:34s} cases={cases[j]:>6} p={p:.2g}",flush=True)
print("\n=== CHECK 2: chr16 4 SVs x MB_280 sign coherence (logOR sign) ===",flush=True)
mb=[j for j in range(J) if phe[j] in ("MB_280","MB_280.1","MB_280.11","MB_280.8","MB_282")]
for s in ("chr16_sv1029_DEL_ori_best6_alt_edge2000","chr16_sv559_DEL_non_best6_alt_edge1000","chr16_sv1017_INS_non_best8_fksize200k","chr16_sv37_INS_ori_best8_midpoint1000"):
    i=sv_ix.get(s)
    if i is None: 
        cand=[x for x in sv if x.startswith(s[:12])][:1]; i=sv_ix.get(cand[0]) if cand else None
    if i is None: continue
    signs=" ".join(f"{phe[j].split('_')[1]}:{'+' if Be[i,j]>0 else '-'}{abs(Be[i,j]):.2f}" for j in mb if np.isfinite(Be[i,j]))
    print(f" {sv[i][:34]:34s} {signs}",flush=True)
print("CHECKS_DONE",flush=True)
