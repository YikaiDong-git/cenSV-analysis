#!/usr/bin/env python
"""The set difference said 296 lost; the counts say 2,245 - 1,962 = 283."""
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
                    usecols=["chrom", "id", "ac", "sv_category"])
    var.append(v)
    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=["sv_id", "sv_category", "svfreq"])
    m["chrom"] = c
    met.append(m)
var = pd.concat(var, ignore_index=True)
met = pd.concat(met, ignore_index=True)

print("duplicate (chrom, id) in sv_variants  : %d" % int(var.duplicated(["chrom", "id"]).sum()))
print("duplicate (chrom, sv_id) in sv_metrics: %d" % int(met.duplicated(["chrom", "sv_id"]).sum()))

vm = var[var.sv_category == "monophyletic"]
mm = met[met.sv_category == "monophyletic"]
print("\nmonophyletic rows      : sv_variants %d, sv_metrics %d" % (len(vm), len(mm)))
print("monophyletic UNIQUE ids: sv_variants %d, sv_metrics %d"
      % (len(set(zip(vm.chrom, vm.id))), len(set(zip(mm.chrom, mm.sv_id)))))

sv_ids, mt_ids = set(zip(vm.chrom, vm.id)), set(zip(mm.chrom, mm.sv_id))
only_v, only_m = sv_ids - mt_ids, mt_ids - sv_ids
print("\nin sv_variants only: %d" % len(only_v))
print("in sv_metrics only : %d" % len(only_m))
if only_m:
    print("  those %d -- what category do they carry in sv_variants?" % len(only_m))
    allv = var.set_index(["chrom", "id"])
    seen = [allv.loc[k].sv_category if k in allv.index else "ABSENT ENTIRELY" for k in only_m]
    print("  %s" % pd.Series([s if isinstance(s, str) else "DUPLICATED" for s in seen])
          .value_counts().to_string())
    for k in list(only_m)[:8]:
        print("    %s  %s" % (k[0], str(k[1])[:70]))

print("\nreconciliation:")
print("  sv_variants monophyletic rows            %5d" % len(vm))
print("  minus those whose id reaches sv_metrics  %5d" % len(sv_ids & mt_ids))
print("  = ids present only in sv_variants        %5d" % len(only_v))
print("  sv_metrics monophyletic rows             %5d" % len(mm))
print("  count difference (2245 - 1962)           %5d" % (len(vm) - len(mm)))
print("  set difference                           %5d" % len(only_v))
print("  the two differ by                        %5d" % (len(only_v) - (len(vm) - len(mm))))
print("  explained by ids in sv_metrics only:     %5d" % len(only_m))
assert len(only_v) - (len(vm) - len(mm)) == len(only_m), "still unreconciled"
print("\n  reconciled: set-difference = count-difference + (sv_metrics-only ids)")
