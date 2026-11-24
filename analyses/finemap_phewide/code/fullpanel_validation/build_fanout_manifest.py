import re
from collections import Counter
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
NEW = CENSV_ROOT + "/analyses/finemap_phewide"
OUT = CENSV_SCRATCH + "/ssv_fanout"
import os; os.makedirs(OUT, exist_ok=True)

traits = {}
for ln in open(NEW + "/manifests/traits_all464.tsv"):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 2:
        traits[p[0]] = p[1]

rows = [l.rstrip("\n").split("\t") for l in open(NEW + "/output/full_v2_catalog.tsv")][1:]
cells = []
for r in rows:
    if len(r) < 13 or r[12] != "DISTINCT_CREDIBLE":
        continue
    region, sv, trait, cv, pm, pi, rr = r[1], r[0], r[2], r[3], r[4], r[5], r[6]
    mreg = re.match(r"chr[\dXY]+_super\d+_(\d+)_(\d+)$", region)
    msv = re.match(r"chr([\dXY]+)_sv(\d+)_", sv)
    if not mreg or not msv:
        print("PARSE_FAIL", region, sv); continue
    rs, re_ = mreg.group(1), mreg.group(2)
    chrom = "chr" + msv.group(1); svidx = msv.group(2)
    bgz = traits.get(trait, "")
    cells.append([region, rs, re_, chrom, sv, svidx, trait, bgz, cv, pm, pi, rr])

cells.sort(key=lambda c: (c[0], c[4], c[6]))
with open(OUT + "/manifest_fanout.tsv", "w") as fh:
    fh.write("region\tregion_start\tregion_end\tchrom\tsv_locus\tsv_index\ttrait\tbgz\tcv\tpip_matched\tpip_impute\tsvSNP_r\n")
    for c in cells:
        fh.write("\t".join(map(str, c)) + "\n")

byreg = Counter(c[0] for c in cells)
n_missing_bgz = sum(1 for c in cells if not c[7])
n_svs = len(set(c[4] for c in cells))
print(f"TOTAL distinct-credible cells: {len(cells)}  unique SVs: {n_svs}  regions: {len(byreg)}  missing bgz: {n_missing_bgz}")
print("cells per region:")
for reg, n in byreg.most_common():
    print(f"  {reg}\t{n}")
