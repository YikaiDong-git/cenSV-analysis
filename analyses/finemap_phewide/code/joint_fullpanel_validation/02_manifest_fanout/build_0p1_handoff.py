#!/usr/bin/env python3
import re, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
M   = CENSV_ROOT + "/analyses/finemap_multi_sv"
NEW = CENSV_ROOT + "/analyses/finemap_phewide"
FAN = CENSV_SCRATCH + "/ssv_joint"
SURV = f"{FAN}/aggregate_0p1/survivors.tsv"
TRAITS = NEW + "/manifests/traits_all464.tsv"
M0P3 = f"{FAN}/manifest_joint_0p3.tsv"

traits = {}
for ln in open(TRAITS):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 2: traits[p[0]] = p[1]

cells_0p3 = set()
with open(M0P3) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        cells_0p3.add((c[ix["region"]], c[ix["sv_locus"]], c[ix["trait"]]))

rx = re.compile(r"^chr([\dXY]+)_super\d+_(\d+)_(\d+)$")
full, new = [], []
seen = set(); miss_bgz = 0
with open(SURV) as fh:
    hdr = fh.readline().rstrip("\n").split("\t"); H = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if c[H["survivor_gc"]] != "True": continue
        region = c[H["region_id"]]; sv_locus = c[H["sv_locus"]]; sv_index = c[H["sv_index"]]; pheno = c[H["pheno"]]
        m = rx.match(region)
        if not m: continue
        key = (region, sv_locus, pheno)
        if key in seen: continue
        seen.add(key)
        chrom = "chr" + m.group(1); rs = m.group(2); re_ = m.group(3)
        bgz = traits.get(pheno, "")
        if not bgz: miss_bgz += 1
        row = [region, rs, re_, chrom, sv_locus, sv_index, pheno, bgz,
               c[H["cv_r2"]], c[H["max_abs_z_snp"]], c[H["marg_z_SV"]], c[H["sv_pip"]]]
        full.append(row)
        if key not in cells_0p3:
            new.append(row)

hdr_out = "region\tregion_start\tregion_end\tchrom\tsv_locus\tsv_index\ttrait\tbgz\tcv\tsvSNP_r\tpip_matched\tpip_impute\n"
full.sort(key=lambda r: (r[0], r[4], r[6])); new.sort(key=lambda r: (r[0], r[4], r[6]))
with open(f"{FAN}/manifest_joint_0p1_FULL.tsv", "w") as fh:
    fh.write(hdr_out); [fh.write("\t".join(map(str, r)) + "\n") for r in full]
with open(f"{FAN}/manifest_joint_0p1.tsv", "w") as fh:
    fh.write(hdr_out); [fh.write("\t".join(map(str, r)) + "\n") for r in new]

def cvbin(r):
    try: return "0.1-0.3" if float(r[8]) < 0.3 else ">=0.3"
    except: return "na"
from collections import Counter
print(f"0p1 survivor_gc cells (FULL): {len(full)}  distinct SVs: {len(set(r[4] for r in full))}  missing bgz: {miss_bgz}")
print(f"  reuse from 0p3 (already fanned): {len(full)-len(new)}")
print(f"  NEW cells to fanout: {len(new)}  distinct SVs: {len(set(r[4] for r in new))}")
print(f"  NEW cells cv2 split: {dict(Counter(cvbin(r) for r in new))}")
print(f"  FULL cells cv2 split: {dict(Counter(cvbin(r) for r in full))}")
