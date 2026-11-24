#!/usr/bin/env python3
import sys

inp = sys.argv[1] if len(sys.argv) > 1 else "/tmp/phewide/inputs/quant.pheno"
out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/phewide/inputs/quant_clean.pheno"

with open(inp) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    loincs = hdr[2:]; L = len(loincs)
    first = [None] * L
    varies = [False] * L
    n = 0
    for line in f:
        p = line.rstrip("\n").split("\t"); vals = p[2:]
        for i in range(L):
            v = vals[i]
            if v == "NA":
                continue
            if first[i] is None:
                first[i] = v
            elif not varies[i] and v != first[i]:
                varies[i] = True
        n += 1
        if n % 50000 == 0:
            print(f"  scanned {n} rows", flush=True)

drop = [i for i in range(L) if not varies[i]]
keep_cols = [0, 1] + [2 + i for i in range(L) if varies[i]]
print(f"scanned {n} rows, {L} LOINCs; dropping {len(drop)} zero-variance: "
      f"{[loincs[i] for i in drop][:20]}{' ...' if len(drop) > 20 else ''}", flush=True)

with open(inp) as f, open(out, "w") as fo:
    for line in f:
        p = line.rstrip("\n").split("\t")
        fo.write("\t".join(p[i] for i in keep_cols) + "\n")
print(f"wrote {out}: {len(keep_cols)-2} LOINCs kept", flush=True)
