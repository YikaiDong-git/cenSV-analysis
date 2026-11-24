#!/usr/bin/env python3
"""keyinterest_svs.tsv -> keyinterest_tasks.tsv: one SLURM task per (region, arm)."""
import sys, re, collections
import numpy as np

svs = sys.argv[1] if len(sys.argv) > 1 else "manifests/keyinterest_svs.tsv"
out = sys.argv[2] if len(sys.argv) > 2 else "manifests/keyinterest_tasks.tsv"
GMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 600

regions = collections.OrderedDict()
for ln in open(svs).read().splitlines()[1:]:
    f = ln.split("\t"); rb, col = f[3], f[4]
    base = rb.split("/")[-1].replace("_bundle.npz", "")
    m = re.match(r"chr(\d+)_super\d+_(\d+)_(\d+)", base)
    chrom, start, end = m.group(1), m.group(2), m.group(3)
    regions.setdefault(rb, [chrom, start, end, []])[3].append(col)

with open(out, "w") as fh:
    fh.write("task\tregion_bundle\tchrom\tstart\tend\tcols\tarm\tmodes\n")
    t = 0
    for rb, (chrom, start, end, cols) in regions.items():
        b = np.load(rb, allow_pickle=True)
        p = int(len(b["snp_pos"])); G = int(len(b["sv_locus"]))
        modes = "hitbundle" if G > GMAX else "allwindow,hitbundle"
        cap = "  <-- HEAVY: all_window CAPPED" if G > GMAX else ""
        print(f"  {rb.split('/')[-1]:44s} p={p:5d} G={G:5d} cols={len(cols)} modes={modes}{cap}")
        for arm in ("matched", "impute"):
            t += 1
            fh.write(f"{t}\t{rb}\t{chrom}\t{start}\t{end}\t{','.join(cols)}\t{arm}\t{modes}\n")
print(f"[{t} tasks ({len(regions)} regions x 2 arms), GMAX={GMAX}] -> {out}")
