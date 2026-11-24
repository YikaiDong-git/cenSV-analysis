#!/usr/bin/env python3
"""A1 step 1: freeze the replication target set."""
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
REPORT = (CENSV_COLLAB_ROOT + "/exploration/meQTL/"
          "reports/chm13_observed_sv_meqtl_mac8_files")
OUT = UNIT + "/output/targets"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


os.makedirs(OUT, exist_ok=True)

pairs = pd.read_csv(REPORT + "/all_discovery_conditioning_consensus.tsv", sep="\t")
print("report pairs: " + str(len(pairs)))
assert len(pairs) == 465, len(pairs)

tg = pairs[["feature_id", "target_chrom", "target_start", "target_end"]].drop_duplicates()
assert (tg.target_end - tg.target_start == 1).all(), "target intervals are not 1 bp"

fid_chrom = tg.feature_id.str.split(":").str[0]
fid_start = tg.feature_id.str.split(":").str[1].str.split("-").str[0].astype(np.int64)
assert (fid_chrom == tg.target_chrom).all(), "feature_id chrom disagrees with target_chrom"
assert (fid_start == tg.target_start).all(), "feature_id start disagrees with target_start"

tg = tg.rename(columns={"target_chrom": "chm13_chrom", "target_start": "chm13_pos"})
tg = tg.drop(columns=["target_end"])
tg["target_id"] = ("t" + pd.Series(np.arange(len(tg)), index=tg.index).astype(str).str.zfill(5))
tg = tg.sort_values(["chm13_chrom", "chm13_pos"]).reset_index(drop=True)
tg["target_id"] = "t" + pd.Series(np.arange(len(tg))).astype(str).str.zfill(5)
print("distinct target CpGs: " + str(len(tg)))
write_atomic(tg[["target_id", "feature_id", "chm13_chrom", "chm13_pos"]],
             OUT + "/targets_cpg.tsv")

var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t")
pairs["source_type"] = pairs.sv_id.str.split("_").str[3]

rows = []
for sv_id, blk in pairs.groupby("sv_id"):
    s = blk.iloc[0]
    hit = var[(var.chrom == s.source_chrom) & (var.svtype == s.source_type)
              & (var.pos == s.source_pos) & (var.end == s.source_end)]
    assert len(hit), "no 0 bp HPRC match for " + sv_id
    for gg in sorted(hit.geno_group.unique()):
        rows.append(dict(sv_id=sv_id, geno_group=int(gg)))
svmap = pd.DataFrame(rows)
print("SV -> genotype-group rows: " + str(len(svmap))
      + "  (" + str(svmap.sv_id.nunique()) + " SVs, "
      + str(svmap.geno_group.nunique()) + " groups)")
write_atomic(svmap.sort_values(["sv_id", "geno_group"]), OUT + "/targets_sv_genogroup.tsv")

pr = pairs.merge(tg[["target_id", "feature_id"]], on="feature_id", how="left")
assert pr.target_id.notna().all()
pr["cohortA_sign"] = np.sign(pr.expected_rint_beta).astype(int)
assert (pr.cohortA_sign != 0).all(), "a cohort-A beta is exactly zero"
pr["dist_class"] = pd.cut(pr.source_target_distance_bp.fillna(-1),
                          [-2, -0.5, 1e5, 1e6, 5e6, 2.5e7, np.inf],
                          labels=["trans", "<=100kb", "0.1-1Mb", "1-5Mb",
                                  "5-25Mb", ">25Mb"])
keep = ["pair_id", "sv_id", "target_id", "feature_id", "scope", "discovery_scope",
        "source_chrom", "source_pos", "source_end", "source_type",
        "chm13_chrom_target", "source_target_distance_bp", "dist_class",
        "expected_rint_beta", "expected_rint_p", "expected_raw_beta", "discovery_q",
        "cohortA_sign", "conditioning_class_scope", "label", "gene_symbol"]
pr["chm13_chrom_target"] = pr.target_chrom
pr = pr[[c for c in keep if c in pr.columns]]
pr = pr.sort_values(["pair_id"]).reset_index(drop=True)
write_atomic(pr, OUT + "/targets_pair.tsv")

exp = (pr.groupby("dist_class", observed=True)
       .agg(n_pairs=("pair_id", "size"),
            cohortA_neg_frac=("cohortA_sign", lambda x: float((x < 0).mean())),
            cohortA_median_beta=("expected_rint_beta", "median")).reset_index())
print("\nPRE-REGISTERED cohort-A expectation (the numbers cohort B must reproduce):")
print(exp.to_string(index=False))
write_atomic(exp, OUT + "/prereg_expectation.tsv")
