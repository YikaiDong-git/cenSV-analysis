#!/usr/bin/env python
import gzip, os, sys
from collections import Counter, defaultdict
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "HG00099_hap1"
CHAIN = f"{P1}/subproject_methylation/raw/chains_chm13/{SAMPLE}_vs_CHM13.chain.gz"
FAI = f"{P1}/input/HPRC_bgzip/{SAMPLE}.fa.gz.fai"
BED = f"{P1}/analyses/cdr_calls/output/calls/{SAMPLE}/binned_freq.bed"

asm, chm, nblk = {}, {}, 0
with gzip.open(CHAIN, "rt") as fh:
    for line in fh:
        if line.startswith("chain"):
            f = line.split()
            asm[f[2]] = int(f[3])
            chm[f[7]] = int(f[8])
            nblk += 1
print(f"{SAMPLE}: {nblk} chain blocks")
print(f"  methylation-assembly contigs in chain : {len(asm)}")
print(f"  CHM13 sequences in chain              : {len(chm)}   e.g. {sorted(chm)[:4]}")

fai = {}
for line in open(FAI):
    f = line.split("\t")
    fai[f[0]] = int(f[1])
print(f"  hifiasm (SV-side) contigs in .fai     : {len(fai)}")

ac, fc = Counter(asm.values()), Counter(fai.values())
uniq = [L for L in set(ac) & set(fc) if ac[L] == 1 and fc[L] == 1]
print("\n--- same assembly under two namings? ---")
print(f"  chain-contig lengths matching a hifiasm contig 1:1 : {len(uniq)}/{len(asm)}")
a_by = {v: k for k, v in asm.items() if ac[v] == 1}
f_by = {v: k for k, v in fai.items() if fc[v] == 1}
for L in sorted(uniq, reverse=True)[:5]:
    print(f"    {a_by[L]:22s} == {f_by[L]:26s} ({L:,} bp)")

used = Counter()
for line in open(BED):
    f = line.split("\t")
    used[f[0].split("#")[-1]] += 1
cov = [c for c in used if c in asm]
print("\n--- do the chains cover the methylation windows? ---")
print(f"  contigs used by binned_freq.bed      : {len(used)}")
print(f"  of those, present in the chain file  : {len(cov)}")
print(f"  windows on chain-covered contigs     : {sum(used[c] for c in cov)}/{sum(used.values())} "
      f"({100*sum(used[c] for c in cov)/max(1,sum(used.values())):.1f}%)")
missing = [c for c in used if c not in asm]
if missing:
    print(f"  NOT in chain: {missing[:6]}")

ok = len(cov) == len(used) and len(uniq) >= 0.9 * len(asm)
print("\nVERDICT: " + ("chains give a COMPLETE CHM13 -> methylation-contig liftover; "
                       "no new alignment needed" if len(cov) == len(used)
                       else "chains do NOT cover all methylation contigs"))
print("         assembly identity: " + ("SAME assembly, two namings (lengths match 1:1)"
                                        if len(uniq) >= 0.9 * len(asm) else "NOT confirmed same"))
