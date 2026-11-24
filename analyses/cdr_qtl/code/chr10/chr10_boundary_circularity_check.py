"""Decisive check for the boundary-circularity concern on the chr10 DEL-340 cdr_pos_offset lead."""
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT=CENSV_ROOT + ""
DISC=f"{ROOT}/analyses/cdr_qtl"
SVD=f"{ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200/chr10"
POS, SVT, SVL = 41048073, "DEL", 340

ph=pd.read_csv(f"{DISC}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
ph=ph[ph.chm13_chrom=="chr10"]
ind=ph.groupby("sample").agg(aa_s=("active_start","mean"), aa_e=("active_end","mean"),
                             aa_len=("active_len","mean"), off=("cdr_pos_offset","mean")).reset_index()
ind["aa_mid"]=(ind.aa_s+ind.aa_e)/2; ind["aa_half"]=ind.aa_len/2
ind["cdr_from_start_kb"]=(ind.aa_half*(1+ind.off))/1000.0
ind["cdr_mid"]=ind.aa_s+ind.aa_half*(1+ind.off)

npz=np.load(f"{SVD}/chr10_sv_genotypes.npz"); G=npz[list(npz.keys())[0]]
samp=open(f"{SVD}/chr10_sv_samples.txt").read().split()
if G.shape[0]!=len(samp) and G.shape[1]==len(samp): G=G.T
svv=pd.read_csv(f"{SVD}/chr10_sv_variants.tsv", sep="\t")
cand=svv[(svv.pos==POS)&(svv.svtype.astype(str).str.startswith(SVT))&(svv.svlen.between(SVL-3,SVL+3))]
dos=pd.Series(np.rint(G[:,cand.index[0]].astype(float)),index=samp).clip(0,2)
d=ind.merge(dos.rename("dos"),left_on="sample",right_index=True)
print(f"n={len(d)}  DEL pos={POS:,}  carriers(dos>0)={int((d.dos>0).sum())}")

print("\n-- DEL position relative to active array [aa_s, aa_e] --")
print(f"  median aa_s={d.aa_s.median():,.0f}  aa_e={d.aa_e.median():,.0f}  DEL@{POS:,}")
frac=((POS-d.aa_s)/(d.aa_e-d.aa_s)).median()
print(f"  DEL fractional position in array (0=start,1=end): {frac:.3f}  -> {'INTERIOR' if 0.05<frac<0.95 else 'BOUNDARY-PROXIMAL'}")
print(f"  distance DEL->aa_s: {(POS-d.aa_s).median()/1e6:+.2f} Mb ; DEL->aa_e: {(POS-d.aa_e).median()/1e6:+.2f} Mb")

print("\n-- what moves with dosage? (Spearman rho, p) --")
for col in ["aa_s","aa_e","aa_len","cdr_mid","cdr_from_start_kb","off"]:
    rho,p=stats.spearmanr(d.dos, d[col])
    g=d.groupby("dos")[col].median()
    gm=" ".join(f"d{int(k)}={v/1e6:.2f}M" if col!="off" and col!="cdr_from_start_kb" else f"d{int(k)}={v:.3f}" for k,v in g.items())
    print(f"  {col:18s} rho={rho:+.3f} p={p:.1e}   [{gm}]")

print("\nVERDICT heuristic:")
rs,ps=stats.spearmanr(d.dos,d.aa_s); re_,pe=stats.spearmanr(d.dos,d.aa_e)
print(f"  aa_s~dos p={ps:.2f} ({'STABLE' if ps>0.05 else 'MOVES'}) ; aa_e~dos p={pe:.2f} ({'STABLE' if pe>0.05 else 'MOVES'})")
print("  -> if aa_s STABLE and the CDR/offset still move, cdr_from_start is NOT circular (start-anchored).")
