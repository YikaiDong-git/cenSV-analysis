#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cross-check 01_build_gao_cdrs.py's merge against an INDEPENDENT implementation."""
import os
import sys
import subprocess
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

B = CENSV_WORKSPACE + "/tools/CDR-Finder"
PY = B + "/env/.pixi/envs/default/bin/python"
CF = (CENSV_ROOT + "/subproject_methylation/"
      "CDR-Finder/workflow/scripts/cdr_finder.py")
SRC = (CENSV_ROOT + "/analyses/"
       "cdr_calls/output/calls")
TMP = "/tmp/gao_xcheck"
GAP = 150_000


def call(binned, out, m):
    subprocess.run([PY, CF, "-i", binned, "-o", out,
                    "--thr_height_perc_valley", "0.34", "--thr_prom_perc_valley", "0.30",
                    "--bp_merge", str(m), "--bp_edge", "500000",
                    "--baseline_avg_methyl", "0.4", "--edge_height_heuristic", "min"],
                   check=True, stderr=subprocess.DEVNULL)


def load(p):
    by = {}
    with open(p) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            by.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    for k in by:
        by[k].sort()
    return by


def merge(iv, gap, strict):
    out = []
    for st, en in iv:
        if out:
            d = st - out[-1][1]
            if (d < gap) if strict else (d <= gap):
                out[-1][1] = max(out[-1][1], en)
                continue
        out.append([st, en])
    return out


os.makedirs(TMP, exist_ok=True)
fails = 0
for tag in sys.argv[1:]:
    binned = os.path.join(SRC, tag, "binned_aa.bed")
    if not os.path.exists(binned):
        print("[skip] %s: no binned_aa.bed" % tag)
        continue
    a, b = TMP + "/%s_m1.bed" % tag, TMP + "/%s_m75k.bed" % tag
    call(binned, a, 1)
    call(binned, b, 75000)
    dips, tool = load(a), load(b)

    inv = sum(1 for ctg in tool for s, e in tool[ctg] if e < s)
    assert inv == 0, "%s: newer CDR-Finder emitted %d inverted intervals" % (tag, inv)

    mine = {ctg: merge(iv, GAP, strict=True) for ctg, iv in dips.items()}
    nd = sum(len(v) for v in dips.values())
    nm = sum(len(v) for v in mine.values())
    nt = sum(len(v) for v in tool.values())
    ok = all([sorted(map(tuple, mine.get(c, []))) == sorted(tool.get(c, []))
              for c in set(mine) | set(tool)])
    print("%-16s dips=%-4d  my_merge=%-4d  tool_bp_merge=%-4d  %s"
          % (tag, nd, nm, nt, "MATCH" if ok else "MISMATCH"))
    if not ok:
        fails += 1
        for c in sorted(set(mine) | set(tool)):
            m_, t_ = sorted(map(tuple, mine.get(c, []))), sorted(tool.get(c, []))
            if m_ != t_:
                print("   %s\n     mine: %s\n     tool: %s" % (c, m_[:6], t_[:6]))

if fails:
    raise SystemExit("CROSS-CHECK FAILED on %d haplotype(s)" % fails)
print("\n[gate] merge implementation reproduces CDR-Finder's own --bp_merge exactly")
