#!/usr/bin/env python3
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
MX = M + "/results/meqtl_cis/chm13_native_same_chrom/matrices"
UNIT = (CENSV_ROOT + "/analyses/"
        "meqtl_replication")
UNIV = UNIT + "/output/null_design/cpg_universe"
OUT = UNIT + "/output/headline_svs"

cpgs = []
for f in sorted(glob.glob(OUT + "/inputs/*_cpgs.tsv")):
    d = pd.read_csv(f, sep="\t")
    d["chrom"] = d.feature_id.astype(str).str.split(":").str[0]
    d["pos"] = d.feature_id.astype(str).str.split(":").str[1].str.split("-").str[0].astype(np.int64)
    cpgs.append(d)
C = pd.concat(cpgs, ignore_index=True)
print("headline CpGs: " + str(len(C)) + "  " + str(dict(C.sv_name.value_counts())))
print("by chromosome: " + str(dict(C.chrom.value_counts())))

SVS = [("chr12_DEL333", "chr12", 37183837), ("chr8_DEL113", "chr8", 44300221)]

avail = {}
for p in sorted(glob.glob(MX + "/*.methylation.tsv.gz")):
    b = os.path.basename(p).replace(".methylation.tsv.gz", "")
    parts = b.rsplit("_", 2)
    if len(parts) != 3:
        continue
    ch, s, e = parts[0], int(parts[1]), int(parts[2])
    avail.setdefault(ch, []).append((s, e, b, p))

need = set()
for ch, sub in C.groupby("chrom"):
    for s, e, b, p in avail.get(ch, []):
        if ((sub.pos >= s) & (sub.pos < e)).any():
            need.add((b, p))
for name, ch, pos in SVS:
    for s, e, b, p in avail.get(ch, []):
        if not (e < pos - 10_000_000 or s > pos + 10_000_000):
            need.add((b, p))

have = set(os.path.basename(f).replace(".cpg.tsv.gz", "")
           for f in glob.glob(UNIV + "/*/*.cpg.tsv.gz"))
print("\nchunks needed: " + str(len(need)) + "   already built: "
      + str(len([b for b, _ in need if b in have])))
missing = sorted([(b, p) for b, p in need if b not in have])
print("missing: " + str(len(missing)))
for b, p in missing:
    print("  " + b + "   " + format(os.path.getsize(p) / 1e6, ".0f") + " MB")

tl = UNIV + "/_tasklist_ext.txt"
with open(tl + ".tmp", "w") as fh:
    for b, p in missing:
        fh.write(p.replace(M + "/", "") + "\n")
os.replace(tl + ".tmp", tl)
print("\nwrote " + tl + "  (" + str(len(missing)) + " chunks)")
print("total input: " + format(sum(os.path.getsize(p) for _, p in missing) / 1e9, ".1f") + " GB")
