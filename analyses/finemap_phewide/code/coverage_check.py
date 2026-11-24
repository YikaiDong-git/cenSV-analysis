#!/usr/bin/env python3
import sys, os
from collections import Counter

NEW = sys.argv[1]
N_TRAITS = 464

blocks = []
with open(os.path.join(NEW, "manifests/full_tasks_v2.tsv")) as fh:
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 8 or p[3] != "hitbundle":
            continue
        blocks.append((p[7], p[4], int(p[5]), int(p[6]), p[1], p[2]))

def diag_path(region, arm, cs, ce):
    shard = "s%02d" % (cs // 500)
    return os.path.join(NEW, "output/full_v2", region, arm, shard, "hitbundle_blk%06d_%d_diag.tsv" % (cs, ce))

def rows(f):
    if not os.path.exists(f):
        return -1
    with open(f) as fh:
        return sum(1 for _ in fh) - 1

missing, examples, complete_ratio = [], [], []
for (region, arm, cs, ce, bun, zm) in blocks:
    exp = (ce - cs) * N_TRAITS
    got = rows(diag_path(region, arm, cs, ce))
    if got < 0 or got < exp - N_TRAITS:
        missing.append((region, arm, cs, ce, bun, zm, got, exp))
        if len(examples) < 8:
            examples.append((region.split('_')[0], arm, cs, ce, got, exp))
    else:
        complete_ratio.append(got / exp if exp else 1.0)

print(f"TOTAL_HITBUNDLE_BLOCKS={len(blocks)}  MISSING_OR_PARTIAL={len(missing)}")
if complete_ratio:
    cr = sorted(complete_ratio)
    print(f"complete-block got/exp: min={cr[0]:.3f} median={cr[len(cr)//2]:.3f} (calibrates whether exp is right)")
print("examples (chr arm cs ce got exp):")
for e in examples:
    print("  ", *e)
by = Counter((m[0].split('_')[0], m[1]) for m in missing)
print("missing by region/arm:")
for k, v in sorted(by.items()):
    print(f"  {k[0]}\t{k[1]}\t{v}")

out = os.path.join(NEW, "manifests/gapfill_tasks.tsv")
with open(out, "w") as fh:
    fh.write("task\tbundle\tzmat\tmode\tarm\tcstart\tcend\tregion\tmem\ttime\n")
    for i, (region, arm, cs, ce, bun, zm, got, exp) in enumerate(missing, 1):
        fh.write(f"{i}\t{bun}\t{zm}\thitbundle\t{arm}\t{cs}\t{ce}\t{region}\t6G\t24:00:00\n")
print(f"GAPFILL_BLOCKS={len(missing)} -> {out}")
