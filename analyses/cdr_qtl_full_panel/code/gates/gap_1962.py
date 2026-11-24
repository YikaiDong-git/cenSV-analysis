#!/usr/bin/env python
"""Why does the paper say 1,962 monophyletic components when the decomposition log has 2,245?"""
import os

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC = ROOT + "/sv_extraction/decompose/" + SUF
EN = ROOT + "/results/elastic_net/decompose/cosine_nj/edge1000_" + SUF
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

log, var, met = [], [], []
for c in CHROMS:
    d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, c, c), sep="\t")
    d["chrom"] = c
    log.append(d[d.type == "monophyletic"])

    v = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, c, c), sep="\t", low_memory=False,
                    usecols=["chrom", "id", "svtype", "svlen", "ac", "af", "maf",
                             "original_sv_id", "mapped_haplogroup", "r2_score", "sv_category"])
    var.append(v)

    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=["sv_id", "sv_category", "cv_r2", "svfreq", "svtype", "svlen"])
    m["chrom"] = c
    met.append(m)

log = pd.concat(log, ignore_index=True)
var = pd.concat(var, ignore_index=True)
met = pd.concat(met, ignore_index=True)

print("1  decomposition log, type=monophyletic : %5d" % len(log))
print("2  sv_variants, sv_category=monophyletic: %5d" % int((var.sv_category == "monophyletic").sum()))
print("3  sv_metrics,  sv_category=monophyletic: %5d   <- the paper's 1,962" % int((met.sv_category == "monophyletic").sum()))
print()
print("all categories, stage 2 vs stage 3:")
t = pd.concat([var.sv_category.value_counts().rename("sv_variants"),
               met.sv_category.value_counts().rename("sv_metrics")], axis=1)
t["lost"] = t.sv_variants - t.sv_metrics
print(t.to_string())

vm = var[var.sv_category == "monophyletic"]
mm = met[met.sv_category == "monophyletic"]
in_met = set(zip(mm.chrom, mm.sv_id))
vm = vm.assign(reached=[(c, i) in in_met for c, i in zip(vm.chrom, vm.id)])
miss = vm[~vm.reached]
print("\n=== the %d monophyletic components that never reached imputation ===" % len(miss))
print("allele count (ac) of those lost vs those kept:")
kept = vm[vm.reached]
for nm, s in (("lost", miss.ac), ("kept", kept.ac)):
    print("  %-5s n=%4d  min %2d  median %4.0f  max %5d" % (nm, len(s), s.min(), s.median(), s.max()))
print("\nlost, by allele count:")
print(miss.ac.value_counts().sort_index().head(12).to_string())
print("\nkept, smallest allele counts:")
print(kept.ac.value_counts().sort_index().head(8).to_string())

thr = kept.ac.min()
print("\nlowest ac that reached imputation: %d" % thr)
print("lost with ac >= %d: %d of %d" % (thr, int((miss.ac >= thr).sum()), len(miss)))
print("lost by chromosome:")
print(miss.groupby("chrom").size().reindex([c for c in CHROMS]).fillna(0).astype(int).to_string())

key = set(zip(miss.chrom, miss.id))
log["cid"] = log.sv_id.astype(str) + "_" + log.branch_id.astype(str)
lost_iter = []
for r in log.itertuples():
    cid = "%s_%s" % (r.sv_id, r.branch_id)
    lost_iter.append(r.iteration)
log["it"] = lost_iter
vm2 = vm.merge(log[["chrom", "sv_id", "branch_id", "iteration", "n_carriers", "r2"]],
               left_on=["chrom", "original_sv_id", "mapped_haplogroup"],
               right_on=["chrom", "sv_id", "branch_id"], how="left")
print("\nof the lost components, which decomposition iteration produced them:")
lo = vm2[~vm2.reached]
print(lo.iteration.value_counts().sort_index().to_string())
print("\nof the kept components:")
print(vm2[vm2.reached].iteration.value_counts().sort_index().to_string())
