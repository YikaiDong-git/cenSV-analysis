#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""How much of the Gao CDR span does our half-prominence call extent miss?"""
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
U = ROOT + "/analyses/cdr_definition"
CANON = ROOT + "/analyses/cdr_calls/output/calls"
EXT = U + "/output/callsets/ext_dips"
GAP = 150_000


def load(p):
    by = {}
    if not os.path.exists(p):
        return by
    for line in open(p):
        c = line.rstrip("\n").split("\t")
        if len(c) < 3:
            continue
        s, e = int(c[1]), int(c[2])
        if e > s:
            by.setdefault(c[0], []).append((s, e))
    for k in by:
        by[k].sort()
    return by


def merge(iv, gap):
    out = []
    for s, e in iv:
        if out and s - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


tags = sorted(os.path.basename(d) for d in glob.glob(EXT + "/*") if os.path.isdir(d))
print("haplotypes re-called with --extend_edges_std -1 : %d" % len(tags))
if len(tags) < 384:
    print("(partial: %d of 384 -- numbers below are over the completed subset)" % len(tags))

os.makedirs(U + "/output/callsets/ext_gao_cdrs", exist_ok=True)
rows = []
n_c = n_e = n_cm = n_em = 0
for tag in tags:
    C = load(os.path.join(CANON, tag, "cdrs.bed"))
    E = load(os.path.join(EXT, tag, "cdrs.bed"))
    if not E:
        continue
    os.makedirs(U + "/output/callsets/ext_gao_cdrs/%s" % tag, exist_ok=True)
    with open(U + "/output/callsets/ext_gao_cdrs/%s/cdrs.bed" % tag, "w") as f:
        for ctg, iv in sorted(E.items()):
            for s, e in merge(iv, GAP):
                f.write("%s\t%d\t%d\n" % (ctg, s, e))
    for ctg in set(C) | set(E):
        c, e = C.get(ctg, []), E.get(ctg, [])
        cm, em = merge(c, GAP) if c else [], merge(e, GAP) if e else []
        n_c += len(c); n_e += len(e); n_cm += len(cm); n_em += len(em)
        rows.append(dict(tag=tag, contig=ctg,
                         n_call_canon=len(c), n_call_ext=len(e),
                         bp_call_canon=sum(y - x for x, y in c),
                         bp_call_ext=sum(y - x for x, y in e),
                         n_cdr_canon=len(cm), n_cdr_ext=len(em),
                         bp_cdr_canon=sum(y - x for x, y in cm),
                         bp_cdr_ext=sum(y - x for x, y in em)))

D = pd.DataFrame(rows)
D.to_csv(U + "/output/callsets/extent_comparison.tsv", sep="\t", index=False)

print()
print("%-34s %14s %14s %10s" % ("", "half-prominence", "edge-extended", "ratio"))
print("-" * 76)
print("%-34s %14d %14d %10.2f" % ("calls before merging", n_c, n_e, n_e / float(n_c)))
print("%-34s %14d %14d %10.2f" % ("CDRs after the 150-kb merge", n_cm, n_em, n_em / float(n_cm)))
tot_c, tot_e = D.bp_call_canon.sum(), D.bp_call_ext.sum()
print("%-34s %14.1f %14.1f %10.2f" % ("total called Mb (calls)", tot_c / 1e6, tot_e / 1e6,
                                      tot_e / float(tot_c)))
tot_cm, tot_em = D.bp_cdr_canon.sum(), D.bp_cdr_ext.sum()
print("%-34s %14.1f %14.1f %10.2f" % ("total called Mb (merged CDRs)", tot_cm / 1e6, tot_em / 1e6,
                                      tot_em / float(tot_cm)))

print()
print("Per active array, median CDR bases after merging:")
print("   half-prominence %8.1f kb    edge-extended %8.1f kb    ratio %.3f"
      % (D.bp_cdr_canon.median() / 1e3, D.bp_cdr_ext.median() / 1e3,
         D.bp_cdr_ext.median() / float(D.bp_cdr_canon.median())))
same_n = float((D.n_cdr_canon == D.n_cdr_ext).mean())
print("   arrays where the CDR COUNT is unchanged by extension: %.1f%%" % (100 * same_n))

print()
print("[written] %s/output/callsets/extent_comparison.tsv" % U)
print("[written] %s/output/callsets/ext_gao_cdrs/  (%d haplotypes)" % (U, len(tags)))
print()
print("Read this as: extension is a large change to a CALL and a small one to a CDR, because a")
print("CDR span is set by the distance between calls, not by where each call's edge is drawn.")
