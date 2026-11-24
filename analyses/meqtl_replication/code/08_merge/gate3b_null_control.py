#!/usr/bin/env python3
"""GATE 3b: is the genotype-dependent missingness SPECIFIC to discovery CpGs, or REGIONAL?"""
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
OUT = UNIT + "/output/merge"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


files = sorted(glob.glob(UNIT + "/output/v2_lift/per_hap/*.tsv.gz"))
D = pd.concat([pd.read_csv(f, sep="\t",
                           usecols=["meth_tag", "cpg_key", "target_class",
                                    "hit_anchor", "n_valid"])
               for f in files], ignore_index=True)
assert D.target_class.notna().all()
D["sample"] = D.meth_tag.str.rsplit("_", n=1).str[0]
D["measured"] = (D.hit_anchor == 1) & (D.n_valid > 0)
print("rows " + str(len(D)) + "  CpGs " + str(D.cpg_key.nunique()))

canon = pd.read_csv(UNIT + "/output/null_design/b0_sv_units/b0_sv_unit_map.tsv", sep="\t")
G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
gs = np.array(open(SUB + "/geno/full_samples.txt").read().split())

common = sorted(set(D["sample"]) & set(gs))
sidx = [int(np.where(gs == s)[0][0]) for s in common]
print("samples with genotype and methylation: " + str(len(common)))

meas = D[D["sample"].isin(common)].groupby(["sample", "cpg_key"], as_index=False).measured.max()
MP = meas.pivot(index="cpg_key", columns="sample", values="measured").fillna(False)
MP = MP[common]

pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
tgt = pd.read_csv(UNIT + "/output/targets/targets_v2_lift.tsv", sep="\t")
f2k = dict(zip(tgt.feature_id, tgt.cpg_key))
pairs["cpg_key"] = pairs.feature_id.map(f2k)
disc_map = pairs[pairs.cpg_key.notna()][["sv_id", "cpg_key"]].drop_duplicates()
disc_map["arm"] = "discovery"

nul = pd.read_csv(UNIT + "/output/null_design/b2_matched_null/b2_null_pairs_dedup.tsv", sep="\t")
null_map = nul[["sv_id", "cpg_key"]].drop_duplicates()
null_map["arm"] = "matched_null"

M = pd.concat([disc_map, null_map], ignore_index=True)
print("test units: " + str(len(M)) + "  (" + str(dict(M.arm.value_counts())) + ")")

rows = []
for sv_id, blk in M.groupby("sv_id"):
    cr = canon[canon.sv_id == sv_id]
    if not len(cr):
        continue
    g = G[int(cr.row_index.iloc[0]), sidx].astype(float)
    ok = g >= 0
    if ok.sum() < 20 or np.std(g[ok]) == 0:
        continue
    gg = g[ok]
    for r in blk.itertuples(index=False):
        if r.cpg_key not in MP.index:
            continue
        y = MP.loc[r.cpg_key].values.astype(float)[ok]
        if y.std() == 0:
            continue
        sl, ic, rr, p, se = stats.linregress(gg, y)
        rows.append(dict(sv_id=sv_id, cpg_key=r.cpg_key, arm=r.arm, n=int(ok.sum()),
                         frac_measured=float(y.mean()), slope=float(sl), p=float(p)))
R = pd.DataFrame(rows)
rej, q, _, _ = multipletests(R.p.values, alpha=0.05, method="fdr_bh")
R["q"] = q
R["sig"] = rej

print("\n=== missingness ~ dosage, BY ARM ===")
z = R.groupby("arm").agg(n=("slope", "size"),
                         median_slope=("slope", "median"),
                         mean_slope=("slope", "mean"),
                         frac_sig=("sig", "mean"),
                         frac_sig_negative=("slope", lambda s: float((s < 0).mean())),
                         median_frac_measured=("frac_measured", "median")).round(4)
print(z.to_string())

d = R[R.arm == "discovery"].slope
n = R[R.arm == "matched_null"].slope
u = stats.mannwhitneyu(d, n, alternative="two-sided")
print("\nMann-Whitney discovery vs null slopes: U=" + format(u.statistic, ".6g")
      + "  p=" + format(u.pvalue, ".3g"))
print("  median slope  discovery " + format(float(d.median()), ".5f")
      + "   null " + format(float(n.median()), ".5f")
      + "   difference " + format(float(d.median() - n.median()), ".5f"))

print("\n=== WITHIN-SV paired contrast (the decisive comparison) ===")
w = R.groupby(["sv_id", "arm"]).slope.median().unstack()
w = w.dropna()
print("  source SVs with both arms: " + str(len(w)))
if len(w) >= 6:
    diff = w["discovery"] - w["matched_null"]
    ws = stats.wilcoxon(w["discovery"], w["matched_null"])
    print("  median within-SV difference (discovery - null): "
          + format(float(diff.median()), ".5f"))
    print("  SVs where discovery is MORE negative: "
          + str(int((diff < 0).sum())) + " / " + str(len(diff)))
    print("  Wilcoxon signed-rank p = " + format(ws.pvalue, ".4g"))

print("\nINTERPRETATION")
print("  If the two arms match, the missingness is REGIONAL: carrier haplotypes lift worse")
print("  across the pericentromeric region, the matched null measures it, and it can be")
print("  subtracted. If discovery is markedly more negative, the discovery CpGs are")
print("  selected where genotype predicts measurability and the association is a collider.")
write_atomic(R, OUT + "/c4_gate3b_missingness_by_arm.tsv")
