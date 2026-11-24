#!/usr/bin/env python
"""Which haplotype order do the recorded carrier indices belong to?"""
import importlib.util
import os

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
VCF = os.path.join(ROOT, "variant_calling/svim_asm", "cohort_svim_asm_%s.vcf.gz" % SUF)
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")
C = "chr9"

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)

matrix, variants, samples = hc.parse_vcf_region(VCF, C, 2, 1.0, None)
res = hc.perform_clustering(matrix, samples, (2, 30), 5, "cosine_nj")
HAPS = res["samples"]
hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
TAB = hg["sample"].astype(str).tolist()

print("parsed samples %d, clustered %d, haplogroup table %d"
      % (len(samples), len(HAPS), len(TAB)))
print("same set: %s" % (set(HAPS) == set(TAB)))
print("same order: %s" % (HAPS == TAB))
bad = [(i, a, b) for i, (a, b) in enumerate(zip(HAPS, TAB)) if a != b]
print("positions differing: %d" % len(bad))
for i, a, b in bad[:8]:
    print("   %4d  clustered=%-16s table=%-16s" % (i, a, b))
print("\nfirst 6 clustered: %s" % HAPS[:6])
print("first 6 in table : %s" % TAB[:6])
print("dropped by the presence filter: %s" % sorted(set(samples) - set(HAPS)))
