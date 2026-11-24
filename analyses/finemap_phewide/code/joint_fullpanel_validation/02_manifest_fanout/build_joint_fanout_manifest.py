#!/usr/bin/env python3
"""Build the JOINT-arm full-panel fanout manifest from the joint 464-scan survivor_gc cells."""
import re, sys, os
from collections import Counter
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]

M   = CENSV_ROOT + "/analyses/finemap_multi_sv"
NEW = CENSV_ROOT + "/analyses/finemap_phewide"
SURV = M + "/code/pipeline/10_block_multisv_464trait_scan/output/aggregate/survivors.tsv"
TRAITS = NEW + "/manifests/traits_all464.tsv"
SSV_MANIFEST = NEW + "/code/fullpanel_validation/manifest_fanout.tsv"
OUT = sys.argv[1] if len(sys.argv) > 1 else CENSV_SCRATCH + "/ssv_joint/manifest_joint_0p3.tsv"

traits = {}
for ln in open(TRAITS):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 2:
        traits[p[0]] = p[1]

ssv_cells = set()
if os.path.exists(SSV_MANIFEST):
    with open(SSV_MANIFEST) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(h)}
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            ssv_cells.add((c[ix["region"]], c[ix["sv_locus"]], c[ix["trait"]]))

rx = re.compile(r"^chr([\dXY]+)_super\d+_(\d+)_(\d+)$")
cells = []
seen = set()
n_missing_bgz = 0
n_parsefail = 0
with open(SURV) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    H = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if c[H["survivor_gc"]] != "True":
            continue
        region = c[H["region_id"]]
        sv_locus = c[H["sv_locus"]]
        sv_index = c[H["sv_index"]]
        pheno = c[H["pheno"]]
        m = rx.match(region)
        if not m:
            n_parsefail += 1
            continue
        chrom = "chr" + m.group(1); rs = m.group(2); re_ = m.group(3)
        key = (region, sv_locus, pheno)
        if key in seen:
            continue
        seen.add(key)
        bgz = traits.get(pheno, "")
        if not bgz:
            n_missing_bgz += 1
        cv = c[H["cv_r2"]]; marg = c[H["marg_z_SV"]]; jpip = c[H["sv_pip"]]; mxz = c[H["max_abs_z_snp"]]
        cells.append([region, rs, re_, chrom, sv_locus, sv_index, pheno, bgz, cv, mxz, marg, jpip])

cells.sort(key=lambda x: (x[0], x[4], x[6]))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("region\tregion_start\tregion_end\tchrom\tsv_locus\tsv_index\ttrait\tbgz\tcv\tsvSNP_r\tpip_matched\tpip_impute\n")
    for c in cells:
        fh.write("\t".join(map(str, c)) + "\n")

n_overlap = sum(1 for c in cells if (c[0], c[4], c[6]) in ssv_cells)
byreg = Counter(c[0] for c in cells)
n_svs = len(set(c[4] for c in cells))
print(f"survivor_gc cells written: {len(cells)}  distinct SVs: {n_svs}  regions: {len(byreg)}")
print(f"missing bgz: {n_missing_bgz}   region-parsefail: {n_parsefail}")
print(f"overlap with single-SV 473 fanout (region,sv_locus,trait): {n_overlap}  (reusable)")
print(f"NEW cells needing fanout: {len(cells) - n_overlap}")
print("cells per region:")
for reg, n in byreg.most_common():
    print(f"  {reg}\t{n}")
