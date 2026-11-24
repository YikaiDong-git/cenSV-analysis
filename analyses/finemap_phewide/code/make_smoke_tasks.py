#!/usr/bin/env python3
"""smoke_svs.tsv -> smoke_tasks.tsv: one SLURM task per (region, arm)."""
import sys, re, collections
svs, out = sys.argv[1], sys.argv[2]
regions = collections.OrderedDict()
for ln in open(svs).read().splitlines()[1:]:
    f = ln.split("\t"); rb, col = f[3], f[4]
    base = rb.split("/")[-1].replace("_bundle.npz", "")
    m = re.match(r"chr(\d+)_super\d+_(\d+)_(\d+)", base)
    chrom, start, end = m.group(1), m.group(2), m.group(3)
    regions.setdefault(rb, [chrom, start, end, []])[3].append(col)
with open(out, "w") as fh:
    fh.write("task\tregion_bundle\tchrom\tstart\tend\tcols\tarm\n")
    t = 0
    for rb, (chrom, start, end, cols) in regions.items():
        for arm in ("matched", "impute"):
            t += 1
            fh.write(f"{t}\t{rb}\t{chrom}\t{start}\t{end}\t{','.join(cols)}\t{arm}\n")
print(f"{t} tasks ({len(regions)} regions x 2 arms) -> {out}")
