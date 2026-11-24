#!/usr/bin/env python3
"""Step 1 of the genotype-complete refit: assemble the pair list and the shard manifest."""
import json
import os

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
CIS = f"{M}/results/meqtl_cis/chm13_observed_sv_robust_mac8_pc10/summary"
TRA = f"{M}/results/meqtl_trans/chm13_observed_sv_all_family_pc10_mac8"
MAT = f"{M}/results/meqtl_cis/chm13_native_same_chrom/matrices"
U = (CENSV_ROOT + "/analyses/"
     "meqtl_replication")
OUT = f"{U}/output/arms/refit_genotype_complete"
os.makedirs(OUT, exist_ok=True)

C = pd.read_csv(f"{CIS}/association_results.tsv.gz", sep="\t", low_memory=False)
C = C[C.rint_hc3_bh_q_primary_t < 0.01].copy()
C["distance_bp"] = (C.methyl_start - C.sv_pos).abs()
C["arm"] = np.where(C.distance_bp <= 5_000_000, "cis_local", "cis_distal")
cis = C[["sv_id", "feature_id", "methyl_chrom", "methyl_start", "arm", "distance_bp",
         "rint_beta", "rint_hc3_t", "rint_hc3_p_t", "n_samples_analyzed"]].rename(
    columns={"sv_id": "family_id", "rint_beta": "pub_rint_beta",
             "rint_hc3_t": "pub_t", "rint_hc3_p_t": "pub_p", "n_samples_analyzed": "pub_n"})
print(f"cis discoveries q<0.01 : {len(cis)}  "
      f"(local {int((cis.arm=='cis_local').sum())}, distal {int((cis.arm=='cis_distal').sum())})")

T = pd.read_csv(f"{TRA}/summary/global_pair_discoveries_physical_family.tsv.gz",
                sep="\t", low_memory=False)
tra = T[["representative_family_id", "feature_id", "methyl_chrom", "methyl_start",
         "pc10_rint_beta", "pc10_rint_hc3_t", "pc10_rint_hc3_p", "n_samples_analyzed"]].rename(
    columns={"representative_family_id": "family_id", "pc10_rint_beta": "pub_rint_beta",
             "pc10_rint_hc3_t": "pub_t", "pc10_rint_hc3_p": "pub_p",
             "n_samples_analyzed": "pub_n"})
tra["arm"] = "trans"
tra["distance_bp"] = np.nan
print(f"trans discoveries      : {len(tra)}")

P = pd.concat([cis, tra], ignore_index=True)
print(f"total pairs            : {len(P)}   distinct CpGs: {P.feature_id.nunique()}")
print(P.arm.value_counts().to_string())

fam = pd.read_csv(f"{TRA}/manifests/family_to_profile.tsv", sep="\t")
f2p = dict(zip(fam.physical_family_id, fam.profile_id))
P["profile_id"] = P.family_id.map(f2p)
miss = P.profile_id.isna().sum()
print(f"pairs with no profile mapping: {miss}")
assert miss == 0, "every discovery SV must map to a dosage profile"

shards = sorted(f for f in os.listdir(MAT) if f.endswith(".methylation.tsv.gz"))
sh = []
for s in shards:
    base = s.replace(".methylation.tsv.gz", "")
    c, a, b = base.rsplit("_", 2)
    sh.append((s, c, int(a), int(b)))
S = pd.DataFrame(sh, columns=["file", "chrom", "start", "end"])


def shard_of(chrom, pos):
    m = S[(S.chrom == chrom) & (S.start <= pos) & (pos < S.end)]
    return m.file.iloc[0] if len(m) else None


P["shard"] = [shard_of(c, p) for c, p in zip(P.methyl_chrom, P.methyl_start)]
nosh = P.shard.isna().sum()
print(f"pairs with no shard: {nosh}")
assert nosh == 0, "every target CpG must fall in a matrix shard"

need = sorted(P.shard.unique())
print(f"shards to read: {len(need)} of {len(S)}")

P.to_csv(f"{OUT}/r1_pairs.tsv", sep="\t", index=False, float_format="%.10g")
with open(f"{OUT}/r2_shards.txt", "w") as fh:
    fh.write("\n".join(need) + "\n")
for s, g in P.groupby("shard"):
    base = s.replace(".methylation.tsv.gz", "")
    c = base.split("_")[0]
    os.makedirs(f"{OUT}/wanted/{c}", exist_ok=True)
    with open(f"{OUT}/wanted/{c}/{base}.ids", "w") as fh:
        fh.write("\n".join(sorted(g.feature_id.unique())) + "\n")

d = np.load(f"{M}/results/meqtl_hierarchical_pc_sensitivity/shared/canonical_pc_design.npz",
            allow_pickle=True)
names = list(d["covariate_names"])
print(f"\ncovariate_names ({len(names)}): {names[:20]} ...")
print(f"pc_counts={list(d['pc_counts'])} covariate_counts={list(d['covariate_counts'])} "
      f"primary_pc={int(d['primary_pc'])}")
z = np.load(f"{TRA}/manifests/profile_dosages.npz", allow_pickle=True)
print(f"design samples == dosage samples: "
      f"{list(d['sample_ids']) == list(z['sample_ids'])}")
json.dump({"n_pairs": int(len(P)), "n_cpgs": int(P.feature_id.nunique()),
           "n_shards": len(need), "arms": P.arm.value_counts().to_dict()},
          open(f"{OUT}/r3_manifest.json", "w"), indent=2)
print(f"\nwrote {OUT}/r1_pairs.tsv, r2_shards.txt, wanted/*.ids")
