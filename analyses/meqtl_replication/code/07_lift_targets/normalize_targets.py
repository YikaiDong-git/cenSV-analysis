#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import pysam
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
REF = PROJ + "/resources/chm13v2.0.fa.gz"
OUT = UNIT + "/output/targets"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


disc = pd.read_csv(UNIT + "/output/targets/targets_cpg.tsv", sep="\t")
disc = disc.rename(columns={"chm13_chrom": "chrom", "chm13_pos": "pos"})
disc["target_class"] = "discovery"
print("discovery targets: " + str(len(disc)))

nul = pd.read_csv(UNIT + "/output/null_design/b2_matched_null/b2_null_targets_cpg.tsv", sep="\t")
nul = nul.rename(columns={"null_target_id": "target_id", "null_feature_id": "feature_id",
                          "null_chrom": "chrom", "null_pos": "pos"})
nul["target_class"] = "matched_null"
print("null targets     : " + str(len(nul)))

T = pd.concat([disc[["target_id", "feature_id", "chrom", "pos", "target_class"]],
               nul[["target_id", "feature_id", "chrom", "pos", "target_class"]]],
              ignore_index=True)
assert T.target_id.is_unique, "target_id collision between discovery and null"
overlap = set(disc.feature_id) & set(nul.feature_id)
assert not overlap, "a null CpG is also a discovery CpG: " + str(list(overlap)[:5])
print("combined         : " + str(len(T)) + "  (ids unique, no discovery/null overlap)")

fa = pysam.FastaFile(REF)
bases, plus_c, dinuc = [], [], []
for t in T.itertuples(index=False):
    p = int(t.pos)
    b = fa.fetch(t.chrom, p, p + 1).upper()
    bases.append(b)
    if b == "C":
        c = p
    elif b == "G":
        c = p - 1
    else:
        c = -1
    plus_c.append(c)
    dinuc.append(fa.fetch(t.chrom, c, c + 2).upper() if c >= 0 else "")
fa.close()

T["chm13_base"] = bases
T["plus_c_pos"] = plus_c
T["chm13_dinuc"] = dinuc

print("\nCHM13 base at the recorded target position:")
print(T.chm13_base.value_counts().to_string())
print("\ndinucleotide at the normalised plus-strand C:")
print(T.chm13_dinuc.value_counts().head(6).to_string())

bad_base = T[~T.chm13_base.isin(["C", "G"])]
bad_cg = T[(T.chm13_base.isin(["C", "G"])) & (T.chm13_dinuc != "CG")]
print("\ntargets whose CHM13 base is neither C nor G: " + str(len(bad_base)))
print("targets that normalise to a non-CG dinucleotide: " + str(len(bad_cg)))
if len(bad_cg):
    print(bad_cg.head(8)[["target_id", "target_class", "chrom", "pos",
                          "chm13_base", "plus_c_pos", "chm13_dinuc"]].to_string(index=False))

T["anchor_ok"] = (T.chm13_base.isin(["C", "G"])) & (T.chm13_dinuc == "CG")
print("\nanchor resolved for " + str(int(T.anchor_ok.sum())) + " / " + str(len(T))
      + " (" + str(round(100.0 * float(T.anchor_ok.mean()), 2)) + "%)")
print("by class:")
print(T.groupby("target_class").anchor_ok.agg(["size", "sum", "mean"]).to_string())

T = T[T.anchor_ok].copy()
T["offset_from_recorded"] = T.plus_c_pos - T.pos
print("\noffset from the recorded position (should be 0 for C, -1 for G):")
print(pd.crosstab(T.chm13_base, T.offset_from_recorded).to_string())

T = T.sort_values(["chrom", "plus_c_pos", "target_id"]).reset_index(drop=True)
write_atomic(T[["target_id", "target_class", "feature_id", "chrom", "pos",
                "chm13_base", "plus_c_pos", "chm13_dinuc", "offset_from_recorded"]],
             OUT + "/targets_v2_normalised.tsv")

dup = T[T.duplicated(["chrom", "plus_c_pos"], keep=False)]
print("\ntargets sharing a normalised CpG with another target: " + str(len(dup)))
if len(dup):
    print("  distinct CpGs involved: " + str(dup.groupby(['chrom', 'plus_c_pos']).ngroups))
    print("  these are ONE measurement, not two -- downstream must not double-count them")
    print(dup.groupby("target_class").size().to_string())
