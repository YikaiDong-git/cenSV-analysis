#!/usr/bin/env python3
import os
import glob
import sys
import numpy as np
import pandas as pd
from scipy import stats
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


os.makedirs(OUT, exist_ok=True)

want_tags = [t.strip() for t in open(UNIT + "/output/a1_lift/all_tags.txt") if t.strip()]
files = sorted(glob.glob(UNIT + "/output/v2_lift/per_hap/*.tsv.gz"))
got_tags = [os.path.basename(f).replace(".tsv.gz", "") for f in files]
print("haplotypes expected " + str(len(want_tags)) + ", present " + str(len(got_tags)))
missing = sorted(set(want_tags) - set(got_tags))
extra = sorted(set(got_tags) - set(want_tags))
if missing or extra:
    print("  missing: " + str(missing[:10]))
    print("  extra  : " + str(extra[:10]))
    sys.exit("ABORT: tag set mismatch -- a partial lift silently biases every downstream "
             "statistic, and v1's loss was biased toward the scarce R10.4.1 stratum")
print("  tag set-equality: OK")

D = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
assert D.target_class.notna().all(), "target_class NaN -- NA-string collision"
print("rows " + str(len(D)) + "  haplotypes " + str(D.meth_tag.nunique())
      + "  CpGs " + str(D.cpg_key.nunique()))
print("chemistry: " + str(dict(D.drop_duplicates('meth_tag').chemistry.value_counts())))

D["measured"] = (D.hit_anchor == 1) & (D.n_valid > 0)
rec = (D.groupby(["cpg_key", "target_class"], as_index=False)
       .agg(n_hap=("measured", "size"), n_measured=("measured", "sum")))
rec["frac_measured"] = rec.n_measured / rec.n_hap
print("\n=== recoverability by target class (fraction of 311 haplotypes measured) ===")
print(rec.groupby("target_class").frac_measured.describe()[
    ["count", "mean", "25%", "50%", "75%"]].round(4).to_string())
write_atomic(rec, OUT + "/c1_recoverability_by_cpg.tsv")

d_med = float(rec[rec.target_class == "discovery"].frac_measured.median())
n_med = float(rec[rec.target_class == "matched_null"].frac_measured.median())
print("\nmedian recoverability: discovery " + str(round(d_med, 4))
      + "   matched_null " + str(round(n_med, 4))
      + "   ratio " + str(round(n_med / max(d_med, 1e-9), 3)))
print("-> b_hat is estimated on the null and applied to discovery; if these differ the")
print("   null must be reweighted onto the discovery recoverability distribution.")

D["sample"] = D.meth_tag.str.rsplit("_", n=1).str[0]
m = D[D.measured]
dip = (m.groupby(["sample", "cpg_key", "target_class"], as_index=False)
       .agg(n_hap_measured=("n_mod", "size"), sum_mod=("n_mod", "sum"),
            sum_valid=("n_valid", "sum"), mean_of_pct=("pct_mod", "mean")))
dip["pct_merged"] = 100.0 * dip.sum_mod / dip.sum_valid
print("\n=== diploid merge ===")
print("  sample x CpG cells: " + str(len(dip))
      + "   both haplotypes measured: " + str(int((dip.n_hap_measured == 2).sum())))
both = dip[dip.n_hap_measured == 2]
if len(both):
    delta = (both.pct_merged - both.mean_of_pct).abs()
    print("  |count-weighted - mean-of-percentages| : median "
          + format(float(delta.median()), ".4g") + ", max " + format(float(delta.max()), ".4g"))
    print("  cells where the two differ by >5 points: " + str(int((delta > 5).sum()))
          + " / " + str(len(both)))
    print("  -> the mean-of-percentages shortcut is wrong wherever coverage is unequal;")
    print("     v1's schema could not have computed the left-hand quantity at all.")
write_atomic(dip, OUT + "/c2_diploid_methylation.tsv")

print("\n=== GATE 3: is missingness genotype-dependent? ===")
pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
tgt = pd.read_csv(UNIT + "/output/targets/targets_v2_lift.tsv", sep="\t")
canon = pd.read_csv(UNIT + "/output/null_design/b0_sv_units/b0_sv_unit_map.tsv", sep="\t")
var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
gsamples = np.array(open(SUB + "/geno/full_samples.txt").read().split())

t_by_feature = dict(zip(tgt.feature_id, tgt.cpg_key))
pairs["cpg_key"] = pairs.feature_id.map(t_by_feature)
pairs = pairs[pairs.cpg_key.notna()]

samples_meth = sorted(D["sample"].unique())
common = [s for s in samples_meth if s in set(gsamples)]
sidx = {s: int(np.where(gsamples == s)[0][0]) for s in common}
print("samples with genotype AND methylation: " + str(len(common)))

meas = (D[D["sample"].isin(common)]
        .groupby(["sample", "cpg_key"], as_index=False).measured.max())
meas_p = meas.pivot(index="cpg_key", columns="sample", values="measured").fillna(False)

rows = []
for sv_id, blk in pairs.groupby("sv_id"):
    cr = canon[canon.sv_id == sv_id]
    if not len(cr):
        continue
    ri = int(cr.row_index.iloc[0])
    g = G[ri, [sidx[s] for s in common]].astype(float)
    ok = g >= 0
    if ok.sum() < 20 or np.nanstd(g[ok]) == 0:
        continue
    for ck in blk.cpg_key.unique():
        if ck not in meas_p.index:
            continue
        y = meas_p.loc[ck, common].values.astype(float)
        yy, gg = y[ok], g[ok]
        if yy.std() == 0:
            rows.append(dict(sv_id=sv_id, cpg_key=ck, n=int(ok.sum()),
                             frac_measured=float(yy.mean()), slope=0.0, p=1.0,
                             note="no_variation_in_missingness"))
            continue
        sl, ic, r, p, se = stats.linregress(gg, yy)
        rows.append(dict(sv_id=sv_id, cpg_key=ck, n=int(ok.sum()),
                         frac_measured=float(yy.mean()), slope=float(sl), p=float(p),
                         note="ok"))
g3 = pd.DataFrame(rows)
if len(g3):
    tested = g3[g3.note == "ok"]
    print("  pairs tested: " + str(len(tested)) + " (of " + str(len(g3)) + " considered)")
    if len(tested):
        from statsmodels.stats.multitest import multipletests
        rej, q, _, _ = multipletests(tested.p.values, alpha=0.05, method="fdr_bh")
        tested = tested.assign(q=q, sig=rej)
        print("  genotype-dependent missingness at BH q<0.05: "
              + str(int(tested.sig.sum())) + " / " + str(len(tested))
              + "  (" + str(round(100.0 * float(tested.sig.mean()), 2)) + "%)")
        print("  slope sign among significant: "
              + str(dict(np.sign(tested[tested.sig].slope).value_counts())))
        print("  median slope (all tested): " + format(float(tested.slope.median()), ".5g"))
        print("\n  READING: a NEGATIVE slope means carriers are LESS likely to be measured --")
        print("  the collider this unit exists to detect. A near-zero rate clears GATE 3.")
        g3 = g3.merge(tested[["sv_id", "cpg_key", "q", "sig"]], on=["sv_id", "cpg_key"],
                      how="left")
    write_atomic(g3, OUT + "/c3_gate3_missingness_vs_genotype.tsv")
else:
    print("  no testable pairs -- investigate before proceeding")
