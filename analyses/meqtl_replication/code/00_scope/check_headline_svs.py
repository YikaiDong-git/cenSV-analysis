#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")

pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
sv = pairs.drop_duplicates("sv_id")[["sv_id", "source_chrom", "source_pos", "source_end",
                                     "source_type"]]
print("source SVs tested in this unit: " + str(len(sv)))
print("by chromosome: " + str(dict(sv.source_chrom.value_counts())))

HEAD = [("chr12 DEL-333 (SV-specific meQTL)", "chr12", 37183837),
        ("chr8 DEL-113 (strongest cis meQTL; CDR lead)", "chr8", 44300221),
        ("chr3 DEL-336 (CENP-B box, trans to chr17)", "chr3", None)]

print("\n=== are the manuscript's headline SVs in this unit's tested set? ===")
for name, chrom, pos in HEAD:
    sub = sv[sv.source_chrom == chrom]
    if pos is None:
        print("\n" + name)
        print("  (position not stated in the results text; listing all " + chrom
              + " SVs tested here)")
        print("  " + chrom + " SVs in unit: " + str(len(sub)))
        if len(sub):
            print(sub[["sv_id", "source_pos", "source_end", "source_type"]].to_string(index=False))
        continue
    exact = sub[(sub.source_pos == pos)]
    near = sub[(sub.source_pos - pos).abs() < 100000]
    print("\n" + name + "   " + chrom + ":" + format(pos, ",d"))
    print("  exact position match : " + str(len(exact)))
    print("  within 100 kb        : " + str(len(near)))
    if len(near):
        n2 = near.assign(dist=(near.source_pos - pos).abs())
        print(n2[["sv_id", "source_pos", "source_end", "source_type", "dist"]]
              .sort_values("dist").to_string(index=False))
    else:
        print("  -> NOT tested by this unit")

print("\n=== what this unit's 47 SVs actually are ===")
print(sv.groupby(["source_chrom", "source_type"]).size().to_string())
