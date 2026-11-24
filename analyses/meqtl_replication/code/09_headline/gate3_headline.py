#!/usr/bin/env python3
"""GATE 3 for the manuscript's two headline meQTL SVs."""
import glob
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
H = UNIT + "/output/headline_svs"

files = sorted(glob.glob(H + "/per_hap/*.tsv.gz"))
want = [t.strip() for t in open(UNIT + "/output/a1_lift/all_tags.txt") if t.strip()]
got = [os.path.basename(f).replace(".tsv.gz", "") for f in files]
print("haplotypes " + str(len(got)) + " / " + str(len(want)))
missing = sorted(set(want) - set(got))
assert not missing, "incomplete lift, missing: " + str(missing[:8])

D = pd.concat([pd.read_csv(f, sep="\t",
                           usecols=["meth_tag", "cpg_key", "target_class", "hit_anchor",
                                    "n_valid", "pct_mod", "lift_status"])
               for f in files], ignore_index=True)
assert D.target_class.notna().all()
D["sample"] = D.meth_tag.str.rsplit("_", n=1).str[0]
D["measured"] = (D.hit_anchor == 1) & (D.n_valid > 0)
print("rows " + str(len(D)) + "   CpGs " + str(D.cpg_key.nunique())
      + "   " + str(dict(D.drop_duplicates('cpg_key').target_class.value_counts())))

T = pd.read_csv(H + "/null/headline_targets_lift.tsv", sep="\t")
N = pd.read_csv(H + "/null/headline_null_pairs.tsv", sep="\t")
B = pd.read_csv(H + "/inputs/branches.tsv", sep="\t")
cpg_sv = dict(zip(T.cpg_key, T.sv_name))
D["sv_name"] = D.cpg_key.map(cpg_sv)
print("unmapped CpGs: " + str(int(D.sv_name.isna().sum())))

G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
gs = np.array(open(SUB + "/geno/full_samples.txt").read().split())
common = sorted(set(D["sample"]) & set(gs))
sidx = [int(np.where(gs == s)[0][0]) for s in common]
print("samples with genotype and methylation: " + str(len(common)))

meas = D[D["sample"].isin(common)].groupby(["sample", "cpg_key"], as_index=False).measured.max()
MP = meas.pivot(index="cpg_key", columns="sample", values="measured").fillna(False)[common]

rows = []
for b in B.itertuples(index=False):
    g = G[int(b.row_index), sidx].astype(float)
    ok = g >= 0
    if ok.sum() < 20 or np.std(g[ok]) == 0:
        print("skip branch " + str(b.geno_group) + ": no variation")
        continue
    gg = g[ok]
    sub = T[T.sv_name == b.sv_name]
    for r in sub.itertuples(index=False):
        if r.cpg_key not in MP.index:
            continue
        y = MP.loc[r.cpg_key].values.astype(float)[ok]
        if y.std() == 0:
            continue
        sl, ic, rr, p, se = stats.linregress(gg, y)
        rows.append(dict(sv_name=b.sv_name, geno_group=int(b.geno_group), af=float(b.af),
                         cpg_key=r.cpg_key, arm=r.target_class, n=int(ok.sum()),
                         frac_measured=float(y.mean()), slope=float(sl), p=float(p)))
R = pd.DataFrame(rows)
if len(R):
    rej, q, _, _ = multipletests(R.p.values, alpha=0.05, method="fdr_bh")
    R["q"] = q
    R["sig"] = rej

print("\n" + "=" * 78)
print("GATE 3 -- missingness ~ dosage, per SV per branch, discovery vs matched null")
print("=" * 78)
z = (R.groupby(["sv_name", "geno_group", "af", "arm"])
     .agg(n_cpg=("slope", "size"), median_slope=("slope", "median"),
          frac_sig=("sig", "mean"), frac_negative=("slope", lambda s: float((s < 0).mean())),
          median_frac_measured=("frac_measured", "median")).round(4))
print(z.to_string())

print("\n=== discovery-minus-null contrast, per branch ===")
for (sv, gg), sub in R.groupby(["sv_name", "geno_group"]):
    d = sub[sub.arm == "discovery"].slope
    n = sub[sub.arm == "matched_null"].slope
    if len(d) < 5 or len(n) < 5:
        print("  " + sv + " grp " + str(gg) + ": too few (" + str(len(d)) + "/" + str(len(n)) + ")")
        continue
    u = stats.mannwhitneyu(d, n, alternative="two-sided")
    print("  " + sv.ljust(13) + " grp " + str(gg).ljust(5)
          + " disc median " + format(float(d.median()), "+.4f")
          + "  null " + format(float(n.median()), "+.4f")
          + "  diff " + format(float(d.median() - n.median()), "+.4f")
          + "  MWU p = " + format(u.pvalue, ".3g")
          + "   (n " + str(len(d)) + "/" + str(len(n)) + ")")

print("\n=== does the conclusion hold across genotype branches? ===")
for sv, sub in R.groupby("sv_name"):
    per = sub[sub.arm == "discovery"].groupby("geno_group").slope.median()
    print("  " + sv + ": discovery median slope by branch = "
          + ", ".join(str(k) + ":" + format(float(v), "+.4f") for k, v in per.items()))
    print("    range across branches: " + format(float(per.max() - per.min()), ".4f")
          + ("   -> branch-robust" if per.max() * per.min() > 0 else
             "   -> SIGN FLIPS ACROSS BRANCHES; not a conclusion"))

print("\n=== reference: the MAC8 cis layer measured earlier ===")
print("  discovery   median slope -0.1584   46.3% significant")
print("  matched null median slope +0.0136    2.6% significant")

tmp = H + "/results/gate3_headline.tsv.tmp"
R.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
os.replace(tmp, H + "/results/gate3_headline.tsv")
print("\nwrote " + H + "/results/gate3_headline.tsv  (" + str(len(R)) + " rows)")
