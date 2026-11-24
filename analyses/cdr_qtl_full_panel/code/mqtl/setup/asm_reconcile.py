#!/usr/bin/env python
"""Are the methylation-pipeline assembly and the SV-calling assembly the same, under two namings?"""
import gzip, os, sys
from collections import Counter
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "HG00099_hap1"
METH_FA = f"{P1}/input/HPRC/{SAMPLE}.fa.gz"
SV_FAI = f"{P1}/input/HPRC_bgzip/{SAMPLE}.fa.gz.fai"


def fasta_lengths(path):
    out, name, ln = {}, None, 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = ln
                name, ln = line[1:].split()[0], 0
            else:
                ln += len(line.strip())
    if name:
        out[name] = ln
    return out


print(f"sample: {SAMPLE}")
meth = fasta_lengths(METH_FA)
print(f"  methylation-pipeline assembly  {METH_FA}")
print(f"    contigs {len(meth)}   total {sum(meth.values()):,} bp")
for k in list(meth)[:3]:
    print(f"      {k}  {meth[k]:,}")

if not os.path.exists(SV_FAI):
    print(f"  SV-side .fai MISSING: {SV_FAI}")
    sys.exit(0)
sv = {}
for line in open(SV_FAI):
    f = line.split("\t")
    sv[f[0]] = int(f[1])
print(f"  SV-calling assembly            {SV_FAI}")
print(f"    contigs {len(sv)}   total {sum(sv.values()):,} bp")
for k in list(sv)[:3]:
    print(f"      {k}  {sv[k]:,}")

print("\n--- reconciliation by exact contig length ---")
mc, sc = Counter(meth.values()), Counter(sv.values())
shared = set(mc) & set(sc)
uniq_1to1 = [L for L in shared if mc[L] == 1 and sc[L] == 1]
print(f"  distinct lengths: methylation {len(mc)}   SV {len(sc)}   shared {len(shared)}")
print(f"  lengths matching UNAMBIGUOUSLY 1:1  {len(uniq_1to1)}")
print(f"  methylation contigs covered by a 1:1 length match: "
      f"{len(uniq_1to1)}/{len(meth)} ({100*len(uniq_1to1)/max(1,len(meth)):.1f}%)")
print(f"  identical total assembly size: {sum(meth.values()) == sum(sv.values())}")

if uniq_1to1:
    m_by = {v: k for k, v in meth.items() if mc[v] == 1}
    s_by = {v: k for k, v in sv.items() if sc[v] == 1}
    print("\n  example name map (methylation -> SV side):")
    for L in sorted(uniq_1to1, reverse=True)[:6]:
        print(f"    {m_by[L]:32s} -> {s_by[L]:26s} ({L:,} bp)")

print("\n--- the contigs the methylation beds actually use ---")
bed = f"{P1}/analyses/cdr_calls/output/calls/{SAMPLE}/binned_freq.bed"
if os.path.exists(bed):
    used = set()
    for line in open(bed):
        used.add(line.split("\t")[0])
    print(f"  {len(used)} contigs used in binned_freq.bed")
    inmeth = sum(1 for c in used if c in meth)
    print(f"  present in the methylation-pipeline assembly: {inmeth}/{len(used)}")
    mapped = sum(1 for c in used if c in meth and mc[meth[c]] == 1 and sc.get(meth[c], 0) == 1)
    print(f"  of those, liftable via a 1:1 length match:    {mapped}/{len(used)}")
    print("  VERDICT: " + ("SAME ASSEMBLY, pure renaming -> base-level liftover is SAFE"
                           if mapped == len(used) and sum(meth.values()) == sum(sv.values())
                           else "NOT a clean renaming -> base-level liftover NOT safe as-is"))
else:
    print(f"  bed missing: {bed}")
