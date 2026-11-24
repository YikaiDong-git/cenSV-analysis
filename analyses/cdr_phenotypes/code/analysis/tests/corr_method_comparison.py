import pandas as pd, numpy as np
from scipy.stats import spearmanr, pearsonr, kendalltau
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
U=CENSV_ROOT + "/analyses/cdr_phenotypes"
ORDER=["chr%d"%i for i in range(1,23)]+["chrX"]
dp=pd.read_csv(U+"/output/phenotypes/per_dip_positions.tsv",sep="\t")
dp=dp[dp.chm13_chrom.isin(ORDER)].copy()
dp["rel"]=(dp.pos_splice+1)/2; dp["L"]=dp.active_len_sum/1e6
hap=(dp.groupby(["sample","hap_local","chm13_chrom"])
       .agg(rel=("rel","mean"),L=("L","first")).reset_index())
A=0.05/23
print("WHAT IS CURRENTLY PLOTTED: rho = SPEARMAN (scipy.stats.spearmanr), rank-based, already.")
print("Below: Pearson r and Kendall tau for comparison, all per haplotype (n = 4,795).")
print()
print("%-7s %10s %10s %10s %9s %s" % ("chrom","Pearson r","Spearman","Kendall","|S|-|P|","sig"))
d=[]
for c in ORDER:
    h=hap[hap.chm13_chrom==c]
    r,pr=pearsonr(h.L,h.rel); s,ps=spearmanr(h.L,h.rel); t,pt=kendalltau(h.L,h.rel)
    d.append((c,r,pr,s,ps,t,pt))
    sig="".join(x for x,pv in (("P",pr),("S",ps),("K",pt)) if pv<A)
    print("%-7s %+10.3f %+10.3f %+10.3f %+9.3f %s" % (c,r,s,t,abs(s)-abs(r),sig))
D=pd.DataFrame(d,columns=["c","r","pr","s","ps","t","pt"])
print()
print("significant at Bonferroni 0.05/23:  Pearson %d   Spearman %d   Kendall %d"
      % ((D.pr<A).sum(),(D.ps<A).sum(),(D.pt<A).sum()))
print("Spearman |rho| exceeds Pearson |r| on %d of 23 chromosomes; mean gain %+.3f"
      % ((D.s.abs()>D.r.abs()).sum(), (D.s.abs()-D.r.abs()).mean()))
print()
fixed=["chr7","chr13","chr14","chr15","chr16"]
sub=D[D.c.isin(fixed)]
print("on the FIXED-ABSOLUTE-OFFSET chromosomes (%s), where rel = const / L is a hyperbola:"
      % ", ".join(fixed))
print("   mean Pearson |r| %.3f   mean Spearman |rho| %.3f   gain %+.3f"
      % (sub.r.abs().mean(), sub.s.abs().mean(), (sub.s.abs()-sub.r.abs()).mean()))
print()
print("skew of array length (Pearson assumes neither variable is heavily skewed):")
print("   pooled skew %.2f ; chromosomes with |skew| > 1: %d of 23"
      % (hap.L.skew(), sum(abs(hap[hap.chm13_chrom==c].L.skew())>1 for c in ORDER)))
