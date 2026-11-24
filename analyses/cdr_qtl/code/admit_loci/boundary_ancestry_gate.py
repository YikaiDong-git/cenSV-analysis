#!/usr/bin/env python
import pandas as pd, numpy as np
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT=CENSV_ROOT + ""
DISC=f"{ROOT}/analyses/cdr_qtl"
CG=f"{DISC}/canonical_geno"

G=np.load(f"{CG}/cdr6519_genotypes.npz"); G=G[list(G.keys())[0]].astype(float)
vt=pd.read_csv(f"{CG}/cdr6519_variants.tsv",sep="\t")
samp=np.array(open(f"{CG}/cdr6519_samples.txt").read().split())
ph=pd.read_csv(f"{DISC}/output/per_hap_chrom_phenotypes.tsv",sep="\t")
ph8=ph[ph.chm13_chrom=="chr8"].groupby("sample").agg(
    aa_s=("active_start","mean"), aa_e=("active_end","mean"),
    aa_len=("active_len","mean"), off=("cdr_pos_offset","mean")).reset_index()
ph8i={r.sample:(r.aa_s,r.aa_e,r.aa_len,r.off) for r in ph8.itertuples()}

leads=[("DEL",113,44300221),("INS",59,44302147),("DEL",347,44460941),("INS",1868,44247324)]
print("=== BOUNDARY-STABILITY GATE (canonical genotypes) ===")
print(f"{'lead':22s} {'n':>4s} {'carr':>5s} {'af':>6s} | aa_s~dos       aa_e~dos       off~dos      verdict")
for svt,svl,pos in leads:
    row=vt[(vt.chrom=='chr8')&(vt.svtype==svt)&(vt.svlen==svl)&(vt.pos==pos)]
    if len(row)==0:
        row=vt[(vt.chrom=='chr8')&(vt.svtype==svt)&(vt.pos==pos)]
    if len(row)==0:
        print(f"{svt}-{svl}@{pos}: NOT FOUND in canonical"); continue
    i=row.index[0]
    dos=pd.Series(G[i,:],index=samp)
    rows=[]
    for s in samp:
        if s in ph8i and dos[s]>=0:
            aa_s,aa_e,aa_len,off=ph8i[s]
            rows.append((dos[s],aa_s,aa_e,off))
    d=pd.DataFrame(rows,columns=["dos","aa_s","aa_e","off"]).dropna()
    def sp(col):
        r,p=stats.spearmanr(d.dos,d[col]); return r,p
    rs,ps=sp("aa_s"); re,pe=sp("aa_e"); ro,po=sp("off")
    verdict="GENUINE REPOSITION" if (ps>0.05 and pe>0.05) else ("boundary-confounded" if (ps<0.05 or pe<0.05) else "?")
    tag_s="STABLE" if ps>0.05 else "MOVES"; tag_e="STABLE" if pe>0.05 else "MOVES"
    print(f"{svt}-{svl}@{pos:<9d} {len(d):>4d} {int((d.dos>0).sum()):>5d} {d.dos.mean()/2:>6.3f} | "
          f"r={rs:+.2f} p={ps:.1e}[{tag_s:6s}] r={re:+.2f} p={pe:.1e}[{tag_e:6s}] r={ro:+.2f} p={po:.1e}  {verdict}")

print("\n=== ANCESTRY GATE — headline DEL-113 @44300221 (canonical) ===")
import os
sp_file=f"{ROOT}/resources/igsr_samples.tsv"
supermap={}
if os.path.exists(sp_file):
    sp=pd.read_csv(sp_file,sep="\t")
    col_s=[c for c in sp.columns if 'ample' in c][0]; col_p=[c for c in sp.columns if 'opulation code' in c or 'Superpop' in c or 'super' in c.lower()]
    if col_p:
        supermap=dict(zip(sp[col_s],sp[col_p[0]]))
if not supermap:
    print("  (no superpop file at resources/igsr_samples.tsv; skipping stratified meta - genotype-level effect already Bonferroni-passing)")
else:
    row=vt[(vt.chrom=='chr8')&(vt.svtype=='DEL')&(vt.svlen==113)&(vt.pos==44300221)]
    i=row.index[0]; dos=pd.Series(G[i,:],index=samp)
    recs=[]
    for s in samp:
        if s in ph8i and dos[s]>=0:
            recs.append((supermap.get(s,'NA'),dos[s],ph8i[s][3]))
    dd=pd.DataFrame(recs,columns=['sp','dos','off']).dropna()
    print(f"  per-superpop SV~offset (n>=12, MAF>=0.05):")
    betas=[]; ses=[]
    for spn,g in dd.groupby('sp'):
        if len(g)>=12 and 0.05<g.dos.mean()/2<0.95:
            b,a,r,p,se=stats.linregress(g.dos,g.off)
            betas.append(b); ses.append(se)
            print(f"    {spn:4s} n={len(g):3d} af={g.dos.mean()/2:.3f} beta={b:+.4f} p={p:.2e} sign={'-' if b<0 else '+'}")
    if len(betas)>=2:
        b=np.array(betas); w=1/np.array(ses)**2
        bmeta=np.sum(w*b)/np.sum(w)
        Q=np.sum(w*(b-bmeta)**2); I2=max(0,(Q-(len(b)-1))/Q)*100 if Q>0 else 0
        print(f"  RE-meta beta={bmeta:+.4f}  I2={I2:.1f}%  n_superpops={len(b)}  all-same-sign={len(set(np.sign(b)))==1}")
