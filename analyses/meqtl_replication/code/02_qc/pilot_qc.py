#!/usr/bin/env python3
"""Pilot QC: explain the two anomalies in the 8-haplotype lift before scaling to 311."""
import os
import glob
import subprocess
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
REF = PROJ + "/resources/chm13v2.0.fa.gz"
ACTIVE_BED = PROJ + "/resources/chm13v2.0_censat_v2.1.active.bed"
OUT = UNIT + "/output/a1_lift/pilot_qc"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(UNIT + "/output/a1_lift/per_hap_pilot/*.tsv.gz"))
print("pilot haplotypes: " + str(len(files)))
L = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
print("rows: " + str(len(L)))
L["hit"] = L.cpg_hit_offset.notna()

tg = pd.read_csv(UNIT + "/output/targets/targets_cpg.tsv", sep="\t")
pr = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
dmin = (pr.groupby("target_id")
        .agg(min_dist=("source_target_distance_bp", "min"),
             dist_class=("dist_class", "first"),
             scope=("scope", "first")).reset_index())
T = tg.merge(dmin, on="target_id", how="left")

ab = pd.read_csv(ACTIVE_BED, sep="\t", header=None, usecols=[0, 1, 2],
                 names=["chrom", "start", "end"])
byc = {}
for c, b in ab.groupby("chrom"):
    byc[c] = (np.sort(b.start.values), np.sort(b.end.values))


def in_active(chrom, p):
    if chrom not in byc:
        return False
    st, en = byc[chrom]
    j = int(np.searchsorted(st, p, side="right")) - 1
    return bool(j >= 0 and p < en[j])


T["in_active_hor"] = [in_active(c, p) for c, p in zip(T.chm13_chrom, T.chm13_pos)]
print("\ntargets inside active HOR: " + str(int(T.in_active_hor.sum())) + "/" + str(len(T)))

M = L.merge(T[["target_id", "in_active_hor", "min_dist", "dist_class", "scope"]],
            on="target_id", how="left")
ok = M[M.lift_status == "LIFTED"]

h1 = (ok.groupby("in_active_hor")
      .agg(n=("hit", "size"), hit_rate=("hit", "mean")).reset_index())
print("\n=== H1: hit rate by active-HOR membership ===")
print(h1.to_string(index=False))
write_atomic(h1, OUT + "/h1_hit_by_active_hor.tsv")

h1b = (ok.groupby(["dist_class", "in_active_hor"], observed=True)
       .agg(n=("hit", "size"), hit_rate=("hit", "mean")).reset_index())
print("\n  by distance class:")
print(h1b.to_string(index=False))
write_atomic(h1b, OUT + "/h1_hit_by_distance.tsv")

uniq = T[["target_id", "chm13_chrom", "chm13_pos"]].drop_duplicates()
regions = "\n".join(r.chm13_chrom + ":" + str(int(r.chm13_pos) + 1) + "-"
                    + str(int(r.chm13_pos) + 2) for r in uniq.itertuples(index=False))
reg_file = OUT + "/h2_regions.txt"
tmp = reg_file + ".tmp." + str(os.getpid())
with open(tmp, "w") as fh:
    fh.write(regions + "\n")
os.replace(tmp, reg_file)
p = subprocess.run(["samtools", "faidx", REF, "-r", reg_file],
                   capture_output=True, text=True)
if p.returncode != 0:
    raise RuntimeError("samtools faidx failed: " + p.stderr[:500])
bases, cur = [], None
for line in p.stdout.splitlines():
    if line.startswith(">"):
        cur = line[1:]
    else:
        bases.append(line.strip().upper())
assert len(bases) == len(uniq), (len(bases), len(uniq))
uniq = uniq.copy()
uniq["ref_dinuc"] = bases
uniq["ref_base"] = uniq.ref_dinuc.str[0]
print("\n=== H2: CHM13 reference base at the target position ===")
print(uniq.ref_base.value_counts().to_string())
print("dinucleotide (base at p and p+1):")
print(uniq.ref_dinuc.value_counts().head(6).to_string())

M2 = ok.merge(uniq[["target_id", "ref_base", "ref_dinuc"]], on="target_id", how="left")
xt = pd.crosstab(M2[M2.hit].ref_base, M2[M2.hit].cpg_hit_offset.astype(int))
print("\n  offset that hit, by CHM13 reference base (forward+reverse pooled):")
print(xt.to_string())
write_atomic(xt.reset_index(), OUT + "/h2_offset_by_refbase.tsv")

fwd = M2[M2.hit & (~M2.is_reverse.astype(bool))]
xtf = pd.crosstab(fwd.ref_base, fwd.cpg_hit_offset.astype(int))
print("\n  FORWARD alignments only:")
print(xtf.to_string())

h2 = (M2.groupby("ref_base").agg(n=("hit", "size"), hit_rate=("hit", "mean")).reset_index())
print("\n  hit rate by reference base:")
print(h2.to_string(index=False))
write_atomic(h2, OUT + "/h2_hit_by_refbase.tsv")

inact = ok[ok.in_active_hor.astype(bool)]
print("\n=== WHAT IS ACTUALLY TESTABLE ===")
print("lifted target-haplotype rows      : " + str(len(ok)))
print("of those, inside active HOR       : " + str(len(inact))
      + "  (hit rate " + str(round(float(inact.hit.mean()) * 100, 1)) + "%)")
per_t = (ok.groupby("target_id")
         .agg(in_active=("in_active_hor", "first"), n_hap=("hit", "size"),
              n_hit=("hit", "sum")).reset_index())
print("distinct targets with >=1 hit     : " + str(int((per_t.n_hit > 0).sum()))
      + " / " + str(len(per_t)))
print("distinct targets hit in ALL 8 hap : " + str(int((per_t.n_hit == per_t.n_hap).sum())))
write_atomic(per_t, OUT + "/target_hit_summary.tsv")

pairs_ok = pr[pr.target_id.isin(set(per_t[per_t.n_hit > 0].target_id))]
print("\npairs whose target is recoverable, by distance class:")
print(pairs_ok.groupby("dist_class", observed=True).size().to_string())
print("\n(compare with the full pre-registered set)")
print(pr.groupby("dist_class", observed=True).size().to_string())
write_atomic(pairs_ok.groupby("dist_class", observed=True).size().reset_index(name="n_recoverable"),
             OUT + "/recoverable_pairs_by_distclass.tsv")
