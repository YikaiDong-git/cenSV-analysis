#!/usr/bin/env python
"""Of the 138 SVs with more than one origin, how many actually got a model for EACH origin?"""
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC = ROOT + "/sv_extraction/decompose/" + SUF
EN = ROOT + "/results/elastic_net/decompose/cosine_nj/edge1000_" + SUF
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

var, met = [], []
for c in CHROMS:
    v = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, c, c), sep="\t", low_memory=False,
                    usecols=["chrom", "id", "ac", "original_sv_id", "sv_category"])
    var.append(v)
    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=["sv_id", "sv_category", "cv_r2"])
    m["chrom"] = c
    met.append(m)
var = pd.concat(var, ignore_index=True)
met = pd.concat(met, ignore_index=True)

vm = var[var.sv_category == "monophyletic"].copy()
mm = met[met.sv_category == "monophyletic"]
in_met = set(zip(mm.chrom, mm.sv_id))
vm["modelled"] = [(c, i) in in_met for c, i in zip(vm.chrom, vm.id)]

per = vm.groupby(["chrom", "original_sv_id"]).agg(
    origins=("id", "size"), modelled=("modelled", "sum")).reset_index()
print("SVs with >=1 monophyletic component : %d" % len(per))
print("SVs with >1 (multiple origins)      : %d" % int((per.origins > 1).sum()))

multi = per[per.origins > 1]
print("\n=== of the %d multi-origin SVs, how many of their origins got a model ===" % len(multi))
tab = multi.groupby(["origins", "modelled"]).size().rename("SVs").reset_index()
print(tab.to_string(index=False))

both = multi[multi.modelled == multi.origins]
some = multi[(multi.modelled > 0) & (multi.modelled < multi.origins)]
none = multi[multi.modelled == 0]
print("\n  every origin modelled : %3d SVs" % len(both))
print("  only some modelled    : %3d SVs" % len(some))
print("  none modelled         : %3d SVs" % len(none))
print("\n  -> the claim 'a model for each origin' holds for %d of %d multi-origin SVs (%.0f%%)"
      % (len(both), len(multi), 100.0 * len(both) / len(multi)))

print("\ncomponents, not SVs:")
print("  monophyletic components from multi-origin SVs : %d" % int(multi.origins.sum()))
print("  of those, modelled                            : %d" % int(multi.modelled.sum()))
first = vm.groupby(["chrom", "original_sv_id"]).ac.idxmax()
print("\nlater-origin components (total minus one per SV): %d"
      % (len(vm) - per.shape[0]))
print("  of those, modelled: %d" % (int(vm.modelled.sum()) - int(
    per.assign(f=per.modelled.clip(upper=1)).f.sum())))

print("\nSVs whose ONLY modelled component is a single origin, i.e. the SV enters the")
print("analysis looking monophyletic when the tree said it was not: %d" % int(
    (multi.modelled == 1).sum()))
