#!/usr/bin/env python3
"""Chunk a fanout manifest into region-contiguous cell ranges (~CELLS_PER_CHUNK each)."""
import sys
from collections import OrderedDict
MANIFEST = sys.argv[1]
OUT = sys.argv[2]
CPC = int(sys.argv[3]) if len(sys.argv) > 3 else 20

byreg = OrderedDict()
with open(MANIFEST) as fh:
    fh.readline()
    for ln in fh:
        reg = ln.split("\t", 1)[0]
        byreg[reg] = byreg.get(reg, 0) + 1

chunks = []
cid = 0
for reg, n in byreg.items():
    s = 1
    while s <= n:
        e = min(s + CPC - 1, n)
        cid += 1
        chunks.append((reg, s, e, f"{cid:04d}"))
        s = e + 1
with open(OUT, "w") as fh:
    for reg, s, e, c in chunks:
        fh.write(f"{reg}\t{s}\t{e}\t{c}\n")
print(f"regions: {len(byreg)}  total cells: {sum(byreg.values())}  chunks: {len(chunks)}  (<= {CPC} cells/chunk)")
