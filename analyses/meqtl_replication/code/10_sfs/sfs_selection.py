#!/usr/bin/env python3
"""Q6: the site-frequency-spectrum test of the negative-selection claim."""
import os
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
OUT = UNIT + "/output/arms/sfs_selection"
os.makedirs(OUT, exist_ok=True)


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


A = pd.read_csv(M + "/results/meqtl_cis/chm13_observed_sv_robust_mac8_pc10/summary/"
                "association_results.tsv.gz", sep="\t")
A["dist"] = (A.methyl_start - A.sv_pos).abs()
print("MAC8 association rows " + str(len(A)) + "   tested SVs " + str(A.sv_id.nunique()))

loc = A[A.dist < 5e6].copy()
loc["neg"] = (loc.rint_beta < 0).astype(int)
d = (loc.groupby(["sv_id", "sv_chrom", "sv_pos", "sv_end", "svtype", "svlen"], as_index=False)
     .agg(n_cpg=("neg", "size"), frac_neg=("neg", "mean"),
          median_beta=("rint_beta", "median")))
d = d[d.n_cpg >= 3].copy()
d["direction"] = np.where(d.frac_neg > 0.5, "decreaser",
                          np.where(d.frac_neg < 0.5, "increaser", "tied"))
print("\nSVs with >=3 local CpGs: " + str(len(d)))
print(d.direction.value_counts().to_string())
print("\nby type:")
print(pd.crosstab(d.svtype, d.direction).to_string())

var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
key = (var.groupby(["chrom", "pos", "end", "svtype"], as_index=False)
       .agg(hprc_af=("af", "median"), n_rows=("af", "size"),
            af_spread=("af", lambda s: float(s.max() - s.min()))))
d = d.merge(key, left_on=["sv_chrom", "sv_pos", "sv_end", "svtype"],
            right_on=["chrom", "pos", "end", "svtype"], how="left")
print("\nSVs matched to the HPRC catalog: " + str(int(d.hprc_af.notna().sum()))
      + " / " + str(len(d)))
print("  multi-branch loci (af spread > 0.05): "
      + str(int((d.af_spread > 0.05).sum())) + "  (af is the branch median for these)")

d = d[d.hprc_af.notna()].copy()
d["maf"] = np.minimum(d.hprc_af, 1 - d.hprc_af)

print("\n" + "=" * 74)
print("1. THE DIRECT TEST: are methylation-INCREASING SVs rarer?")
print("=" * 74)
inc = d[d.direction == "increaser"]
dec = d[d.direction == "decreaser"]
print("  increasers n=" + str(len(inc)) + "   decreasers n=" + str(len(dec)))
for lab, col in [("HPRC AF", "hprc_af"), ("HPRC MAF", "maf")]:
    if len(inc) > 3 and len(dec) > 3:
        u = stats.mannwhitneyu(inc[col], dec[col], alternative="less")
        print("  " + lab.ljust(9) + "  increaser median "
              + format(float(inc[col].median()), ".4f")
              + "   decreaser median " + format(float(dec[col].median()), ".4f")
              + "   MWU (increaser<decreaser) p = " + format(u.pvalue, ".4g"))

print("\n  WITHIN TYPE (not confounded by DEL/INS spectra):")
for t in ["DEL", "INS"]:
    i2 = inc[inc.svtype == t]
    d2 = dec[dec.svtype == t]
    if len(i2) >= 3 and len(d2) >= 3:
        u = stats.mannwhitneyu(i2.hprc_af, d2.hprc_af, alternative="less")
        print("    " + t + ": increaser n=" + str(len(i2)) + " median AF "
              + format(float(i2.hprc_af.median()), ".4f")
              + "   decreaser n=" + str(len(d2)) + " median AF "
              + format(float(d2.hprc_af.median()), ".4f")
              + "   p = " + format(u.pvalue, ".4g"))
    else:
        print("    " + t + ": too few (" + str(len(i2)) + " increasers, "
              + str(len(d2)) + " decreasers)")

print("\n" + "=" * 74)
print("2. CONTINUOUS: does degree of methylation-lowering track frequency?")
print("=" * 74)
for lab, sub in [("all", d), ("DEL only", d[d.svtype == "DEL"]),
                 ("INS only", d[d.svtype == "INS"])]:
    if len(sub) > 6:
        r1, p1 = stats.spearmanr(sub.frac_neg, sub.hprc_af)
        r2, p2 = stats.spearmanr(sub.median_beta, sub.hprc_af)
        print("  " + lab.ljust(9) + " n=" + str(len(sub))
              + "  Spearman(frac_neg, AF) = " + format(r1, "+.4f") + " p=" + format(p1, ".3g")
              + "   Spearman(median_beta, AF) = " + format(r2, "+.4f")
              + " p=" + format(p2, ".3g"))
print("\n  Under the claim, methylation-INCREASING SVs are purged, so AF should RISE with")
print("  frac_neg. A null or negative correlation is evidence against.")

print("\n" + "=" * 74)
print("3. RARE-VARIANT FRACTION")
print("=" * 74)
for thr in (0.05, 0.10):
    tab = np.array([[int((inc.maf < thr).sum()), int((inc.maf >= thr).sum())],
                    [int((dec.maf < thr).sum()), int((dec.maf >= thr).sum())]])
    orr, pv = stats.fisher_exact(tab, alternative="greater")
    a = float((inc.maf < thr).mean()) if len(inc) else np.nan
    b = float((dec.maf < thr).mean()) if len(dec) else np.nan
    print("  MAF < " + str(thr) + ":  increasers " + format(a, ".3f")
          + "   decreasers " + format(b, ".3f")
          + "   Fisher OR " + format(orr, ".3f") + "  p = " + format(pv, ".4g"))

print("\n" + "=" * 74)
print("4. LENGTH- AND TYPE-MATCHED COMPARISON")
print("=" * 74)
rows = []
for r in inc.itertuples(index=False):
    pool = dec[(dec.svtype == r.svtype)
               & (np.abs(np.log10(dec.svlen.clip(lower=1))
                         - np.log10(max(r.svlen, 1))) <= 0.3)]
    if len(pool) >= 2:
        rows.append(dict(sv_id=r.sv_id, svtype=r.svtype, svlen=int(r.svlen),
                         inc_af=float(r.hprc_af),
                         matched_dec_median_af=float(pool.hprc_af.median()),
                         n_matched=len(pool)))
Mt = pd.DataFrame(rows)
if len(Mt) >= 5:
    diff = Mt.inc_af - Mt.matched_dec_median_af
    w = stats.wilcoxon(Mt.inc_af, Mt.matched_dec_median_af)
    print("  matched increasers: " + str(len(Mt)))
    print("  median (increaser AF - matched decreaser AF) = "
          + format(float(diff.median()), "+.4f"))
    print("  increasers rarer than their matches in " + str(int((diff < 0).sum()))
          + " / " + str(len(diff)) + " cases")
    print("  Wilcoxon signed-rank p = " + format(w.pvalue, ".4g"))
    print("\n  The claim predicts a NEGATIVE median difference and most increasers rarer.")
else:
    print("  too few matched increasers (" + str(len(Mt)) + ") -- not estimable")

print("\n" + "=" * 74)
print("5. BASELINE vs the whole catalog")
print("=" * 74)
allsv = var.drop_duplicates(["chrom", "pos", "end", "svtype"])
print("  full catalog AF by type:")
print(allsv.groupby("svtype").af.describe()[["count", "25%", "50%", "75%"]].round(4).to_string())
print("  tested-and-classified SVs, by type and direction:")
print(d.groupby(["svtype", "direction"]).hprc_af
      .describe()[["count", "25%", "50%", "75%"]].round(4).to_string())

write_atomic(d, OUT + "/s1_sv_direction_af.tsv")
if len(Mt):
    write_atomic(Mt, OUT + "/s2_matched_comparison.tsv")

print("\n" + "=" * 74)
print("CAVEAT, to be stated wherever this is reported")
print("=" * 74)
print("  Every SV here passed MAC >= 8 to enter the scan, so the rarest variants -- the ones")
print("  purifying selection would push hardest on -- are absent by construction. This tests")
print("  the claim among variants common enough to be tested, which is where the claim was")
print("  made, but it cannot see the low-frequency tail.")
