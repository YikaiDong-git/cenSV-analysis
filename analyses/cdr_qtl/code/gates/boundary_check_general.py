import sys, numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT=CENSV_ROOT + ""
DISC=f"{ROOT}/analyses/cdr_qtl"
SVD=f"{ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
CHROM, POS, SVT, SVL = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])

ph=pd.read_csv(f"{DISC}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
ph=ph[ph.chm13_chrom==CHROM]
ind=ph.groupby("sample").agg(aa_s=("active_start","mean"), aa_e=("active_end","mean"),
                             aa_len=("active_len","mean"), off=("cdr_pos_offset","mean")).reset_index()
ind["aa_half"]=ind.aa_len/2
ind["cdr_mid"]=ind.aa_s+ind.aa_half*(1+ind.off)
d0=f"{SVD}/{CHROM}"
npz=np.load(f"{d0}/{CHROM}_sv_genotypes.npz"); G=npz[list(npz.keys())[0]]
samp=open(f"{d0}/{CHROM}_sv_samples.txt").read().split()
if G.shape[0]!=len(samp) and G.shape[1]==len(samp): G=G.T
svv=pd.read_csv(f"{d0}/{CHROM}_sv_variants.tsv", sep="\t")
cand=svv[(svv.pos==POS)&(svv.svtype.astype(str).str.startswith(SVT))&(svv.svlen.between(SVL-3,SVL+3))]
if len(cand)==0:
    print(f"NO MATCH for {CHROM} {SVT}-{SVL}@{POS}; nearby {SVT} within 5kb:")
    print(svv[(svv.svtype.astype(str).str.startswith(SVT))&(svv.pos.between(POS-5000,POS+5000))][["pos","svtype","svlen","af"]].to_string(index=False)); sys.exit()
dos=pd.Series(np.rint(G[:,cand.index[0]].astype(float)),index=samp).clip(0,2)
d=ind.merge(dos.rename("dos"),left_on="sample",right_index=True)
frac=((POS-d.aa_s)/(d.aa_e-d.aa_s)).median()
print(f"{CHROM} {SVT}-{SVL}@{POS:,}  n={len(d)} carriers={int((d.dos>0).sum())}  af={d.dos.mean()/2:.3f}")
print(f"  DEL/SV array-fraction={frac:.3f} ({'INTERIOR' if 0.05<frac<0.95 else 'BOUNDARY-PROXIMAL'})")
print("  --- what moves with dosage (Spearman rho, p) ---")
for col in ["aa_s","aa_e","aa_len","cdr_mid","off"]:
    rho,p=stats.spearmanr(d.dos, d[col])
    tag=""
    if col in ("aa_s","aa_e"): tag=" <-- BOUNDARY " + ("STABLE" if p>0.05 else "MOVES")
    print(f"    {col:8s} rho={rho:+.3f} p={p:.1e}{tag}")
