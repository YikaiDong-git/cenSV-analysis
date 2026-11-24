#!/usr/bin/env python
"""Two loose ends in the 2,245 -> 1,962 gap."""
import re

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC = ROOT + "/sv_extraction/decompose/" + SUF
EN = ROOT + "/results/elastic_net/decompose/cosine_nj/edge1000_" + SUF
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
STRIP = re.compile(r"_(B\d+|residual|original)$")

var, met = [], []
for c in CHROMS:
    v = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, c, c), sep="\t", low_memory=False,
                    usecols=["chrom", "pos", "end", "id", "svtype", "svlen", "ac",
                             "original_sv_id", "sv_category"])
    var.append(v)
    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=["sv_id", "sv_category", "svpos", "svfreq"])
    m["chrom"] = c
    met.append(m)
var = pd.concat(var, ignore_index=True)
met = pd.concat(met, ignore_index=True)

vm = var[var.sv_category == "monophyletic"]
mm = met[met.sv_category == "monophyletic"]
only_m = set(zip(mm.chrom, mm.sv_id)) - set(zip(vm.chrom, vm.id))

print("=== (a) the %d ids present only in sv_metrics ===" % len(only_m))
print("chromosomes: %s" % sorted({c for c, _ in only_m}))
parents = {(c, STRIP.sub("", str(i))) for c, i in only_m}
have_parent = {(c, p) for c, p in zip(var.chrom, var.original_sv_id.astype(str))}
hit = sum(1 for k in parents if k in have_parent)
print("their parent SV is present in sv_variants: %d of %d" % (hit, len(parents)))
b_only = {re.search(r"_(B\d+)$", str(i)).group(1) for c, i in only_m
          if re.search(r"_(B\d+)$", str(i))}
print("branches involved: %s" % sorted(b_only))
c10 = var[(var.chrom == "chr10") & (var.sv_category == "monophyletic")]
print("chr10 monophyletic branches in the reconstructed sv_variants: %d distinct"
      % c10.id.map(lambda s: STRIP.search(str(s)).group(1) if STRIP.search(str(s)) else "?")
        .nunique())
m10 = mm[mm.chrom == "chr10"]
print("chr10 monophyletic in sv_metrics %d, in sv_variants %d" % (len(m10), len(c10)))
shared = set(m10.sv_id) & set(c10.id)
print("  ids shared: %d | metrics-only: %d | variants-only: %d"
      % (len(shared), len(set(m10.sv_id) - set(c10.id)), len(set(c10.id) - set(m10.sv_id))))

print("\n=== (b) what drops the 186 with ac >= 10 ===")
in_met = set(zip(mm.chrom, mm.sv_id))
vm = vm.assign(reached=[(c, i) in in_met for c, i in zip(vm.chrom, vm.id)])
miss = vm[~vm.reached & (vm.ac >= 10)]
kept = vm[vm.reached]
print("dropped with ac >= 10: %d   kept: %d" % (len(miss), len(kept)))
print("\nby chromosome (dropped ac>=10 / kept):")
t = pd.concat([miss.groupby("chrom").size().rename("dropped"),
               kept.groupby("chrom").size().rename("kept")], axis=1).fillna(0).astype(int)
t["pct_dropped"] = (100.0 * t.dropped / (t.dropped + t.kept)).round(1)
print(t.reindex([c for c in CHROMS if c in t.index]).to_string())

print("\nposition: do the dropped sit outside the window the kept occupy?")
for c in ["chr1", "chr19", "chr10"]:
    k = kept[kept.chrom == c]
    d = miss[miss.chrom == c]
    if not len(d):
        continue
    print("  %-6s kept pos %d..%d | dropped pos %d..%d | dropped inside kept range: %d of %d"
          % (c, k.pos.min(), k.pos.max(), d.pos.min(), d.pos.max(),
             int(((d.pos >= k.pos.min()) & (d.pos <= k.pos.max())).sum()), len(d)))

print("\nis a dropped component's PARENT also dropped? (whole-SV exclusion vs per-component)")
kept_par = set(zip(kept.chrom, kept.original_sv_id.astype(str)))
miss_par = [(c, p) in kept_par for c, p in zip(miss.chrom, miss.original_sv_id.astype(str))]
print("  dropped components whose SV has another component that DID reach: %d of %d"
      % (sum(miss_par), len(miss)))
print("  -> %d dropped belong to SVs excluded entirely" % (len(miss) - sum(miss_par)))

print("\nsvtype of dropped ac>=10 vs kept:")
print(pd.concat([miss.svtype.value_counts().rename("dropped"),
                 kept.svtype.value_counts().rename("kept")], axis=1).fillna(0).astype(int)
      .to_string())
