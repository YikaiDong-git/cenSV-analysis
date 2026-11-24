#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
OUT = UNIT + "/output/targets"
ND = UNIT + "/output/null_design"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


T = pd.read_csv(OUT + "/targets_v2_normalised.tsv", sep="\t")
assert T.target_class.notna().all(),     "target_class contains NaN -- a class token collided with a pandas NA string"
assert set(T.target_class.unique()) == {"discovery", "matched_null"},     "unexpected target_class values: " + str(sorted(T.target_class.unique()))
print("class counts as read: " + str(dict(T.target_class.value_counts())))

T["cpg_key"] = T.chrom.astype(str) + ":" + T.plus_c_pos.astype(str)
print("targets in: " + str(len(T)) + "   unique CpGs: " + str(T.cpg_key.nunique()))

disc_keys = set(T[T.target_class == "discovery"].cpg_key)
leak = T[(T.target_class == "matched_null") & (T.cpg_key.isin(disc_keys))]
print("\n=== LEAK AUDIT ===")
print("null targets that ARE a discovery CpG on the opposite strand: " + str(len(leak)))
if len(leak):
    print("  distinct CpGs affected: " + str(leak.cpg_key.nunique()))
    print(leak.head(6)[["target_id", "feature_id", "chrom", "pos", "plus_c_pos",
                        "chm13_base"]].to_string(index=False))
    print("  -> these are removed: a null cannot be the CpG it is a null for.")

T = T[~((T.target_class == "matched_null") & (T.cpg_key.isin(disc_keys)))].copy()

T["cls_rank"] = (T.target_class != "discovery").astype(int)
T = T.sort_values(["cpg_key", "cls_rank", "target_id"])
rep = T.groupby("cpg_key", as_index=False).first()
members = (T.groupby("cpg_key").target_id.apply(lambda s: ";".join(sorted(s)))
           .rename("member_target_ids").reset_index())
rep = rep.merge(members, on="cpg_key")
rep["n_members"] = rep.member_target_ids.str.count(";") + 1

print("\nafter leak removal and collapse:")
print("  unique CpGs: " + str(len(rep)))
print(rep.groupby("target_class").size().to_string())
print("  CpGs represented by 2 original targets: " + str(int((rep.n_members > 1).sum())))

write_atomic(rep[["cpg_key", "target_id", "target_class", "feature_id", "chrom",
                  "plus_c_pos", "chm13_base", "member_target_ids", "n_members"]],
             OUT + "/targets_v2_lift.tsv")

np_ = pd.read_csv(ND + "/b2_matched_null/b2_null_pairs.tsv", sep="\t")
nt = pd.read_csv(ND + "/b2_matched_null/b2_null_targets_cpg.tsv", sep="\t")
nt["cpg_key_src"] = nt.null_chrom.astype(str) + ":" + nt.null_pos.astype(str)

key_by_feature = dict(zip(T.feature_id, T.cpg_key))
np_["cpg_key"] = np_.null_feature_id.map(key_by_feature)
dropped = int(np_.cpg_key.isna().sum())
print("\nnull pair rows whose CpG was removed as a leak/unresolved: " + str(dropped))
np_ = np_[np_.cpg_key.notna()].copy()

before = np_.groupby("matched_to_pair_id").size()
np_ = np_.drop_duplicates(["matched_to_pair_id", "cpg_key"])
after = np_.groupby("matched_to_pair_id").size()
print("nulls per discovery pair, before collapse: min " + str(int(before.min()))
      + " median " + str(int(before.median())))
print("nulls per discovery pair, after  collapse: min " + str(int(after.min()))
      + " median " + str(int(after.median())))
print("  discovery pairs now below 10 nulls: " + str(int((after < 10).sum()))
      + " / " + str(len(after)))
print("  discovery pairs now below 5 nulls : " + str(int((after < 5).sum())))

write_atomic(np_, ND + "/b2_matched_null/b2_null_pairs_dedup.tsv")
print("\nfinal: " + str(len(np_)) + " null pair rows over "
      + str(np_.cpg_key.nunique()) + " distinct null CpGs for "
      + str(np_.matched_to_pair_id.nunique()) + " discovery pairs")
