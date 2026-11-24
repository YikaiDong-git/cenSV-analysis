import csv, numpy as np, sys
sys.path.insert(0,"/home/jupyter/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/code/phewide")
import step30g_numpy_disease_scoretest as S
AR="/home/jupyter/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication"
MAT=f"{AR}/results_vds/step30_phewide/disease_matrices"; D=f"{AR}/results_vds/step30_phewide"
Pv=np.load(f"{MAT}/prefix_degenerate/disease_Pv.npy")
meta=list(csv.DictReader(open(f"{MAT}/prefix_degenerate/disease_meta.tsv"),delimiter="\t"))
phe=[m["phecode"] for m in meta]; cases=np.array([int(m["cases"]) for m in meta])
tested=np.array([m["tested"]=="1" for m in meta]); J=len(phe)
with np.errstate(invalid="ignore"):
    frac4=np.array([(np.nan_to_num(Pv[:,j],nan=1.0)<1e-4).mean() for j in range(J)])
susp=[j for j in range(J) if tested[j] and frac4[j]>0.004]
susp=sorted(susp,key=lambda j:-frac4[j])
print(f"suspect phenos (frac@1e-4>0.004): {len(susp)}",flush=True)
iids,iidx=S.load_order(); X,_=S.load_X(iidx); N=len(iids)
hdr=open(f"{S.INP}/disease.pheno").readline().rstrip("\n").split("\t")
cols={phe[j]:hdr.index(phe[j]) for j in susp if phe[j] in hdr}
Y={ph:np.full(N,-1,np.int8) for ph in cols}
for line in open(f"{S.INP}/disease.pheno"):
    p=line.rstrip("\n").split("\t"); kk=iidx.get(p[1])
    if kk is None or p[0]=="FID": continue
    for ph,c in cols.items():
        if p[c]=="0": Y[ph][kk]=0
        elif p[c]=="1": Y[ph][kk]=1
print("phecode      cases  frac@1e-4  |Sigma_r|  conv  DEGEN(|Sr|>1)?",flush=True)
degen=[]
for j in susp:
    ph=phe[j]
    if ph not in Y: continue
    y=Y[ph]; present=y>=0; yy=np.where(present,y,0).astype(float)
    _,_,r,conv=S.fit_null(X,yy,present); sr=abs(float(r.sum()))
    dg = sr>S.IDENT_TOL
    if dg: degen.append(ph)
    print(f"{ph:12s} {int(cases[j]):>5}  {frac4[j]:.4f}    {sr:9.2f}  {str(conv):5s} {'YES' if dg else 'no'}",flush=True)
print(f"\nTRUE degenerate (|Sr|>1) among suspects: {len(degen)} -> {degen}",flush=True)

print("\n=== Firth confirmation gated on CONVERGENCE (from firthrefit) ===",flush=True)
rf=list(csv.DictReader(open(f"{D}/phewide_disease_full_firthrefit.tsv"),delimiter="\t"))
rf=[r for r in rf if r["phecode"]!="PP_907.3"]
def b(r):
    try: return int(r["cases"])
    except: return 0
for lo,hi in [(50,200),(200,500),(500,5000),(5000,10**9)]:
    s=[r for r in rf if lo<=b(r)<hi]
    if not s: continue
    conf=sum(1 for r in s if r["firth_confirmed"]=="1")
    convd=sum(1 for r in s if r["firth_converged"]=="1")
    confconv=sum(1 for r in s if r["firth_confirmed"]=="1" and r["firth_converged"]=="1")
    print(f" cases[{lo},{hi}): n={len(s):4d}  confirmed={100*conf/len(s):5.1f}%  firth_converged={100*convd/len(s):5.1f}%  confirmed&converged={100*confconv/len(s):5.1f}%",flush=True)
alln=len(rf); ac=sum(1 for r in rf if r["firth_confirmed"]=="1"); acc=sum(1 for r in rf if r["firth_confirmed"]=="1" and r["firth_converged"]=="1")
print(f" ALL non-degenerate: confirmed={100*ac/alln:.1f}%  confirmed&converged={100*acc/alln:.1f}%  ({acc}/{alln})",flush=True)
print("AUDIT_DONE",flush=True)
