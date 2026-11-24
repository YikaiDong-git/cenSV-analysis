#!/usr/bin/env python3
"""Step 1 of the AoU imputed-SV replication: variant match + per-SV coverage."""
import csv
import os
from collections import defaultdict

HOME = os.path.expanduser("~")
HANDOFF = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_handoff"
PGEN_DIR = f"{HOME}/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen"
OUT_DIR = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results/step1"
os.makedirs(OUT_DIR, exist_ok=True)

rows = list(csv.DictReader(open(f"{HANDOFF}/sv_weights_long.tsv"), delimiter="\t"))
by_chrom = defaultdict(list)
for r in rows:
    by_chrom[r["chrom"]].append(r)

sv_total = defaultdict(int)
sv_biallelic = defaultdict(int)
sv_recoverable = defaultdict(int)
sv_swapped = defaultdict(int)
sv_chrom = {}
for r in rows:
    sv_total[r["sv_locus"]] += 1
    sv_chrom[r["sv_locus"]] = r["chrom"]


def census_chrom(chrom, wrows):
    want = defaultdict(list)
    for r in wrows:
        want[r["pos"]].append((r["sv_locus"], r["ref"], r["alt"]))
    pvar = f"{PGEN_DIR}/acaf_threshold.{chrom}.pvar"
    seen = {}
    with open(pvar) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t", 5)
            pos = f[1]
            if pos not in want:
                continue
            ref, alt = f[3], f[4]
            seen[pos] = (ref, alt.split(","))
    for pos, items in want.items():
        rec = seen.get(pos)
        for sv, ref, alt in items:
            if rec is None:
                continue
            pref, palts = rec
            biallelic = len(palts) == 1
            fwd = (ref == pref and alt in palts)
            rev = (alt == pref and ref in palts)
            if fwd or rev:
                sv_recoverable[sv] += 1
                if rev and not fwd:
                    sv_swapped[sv] += 1
                if biallelic and (
                    (ref == pref and alt == palts[0]) or (alt == pref and ref == palts[0])
                ):
                    sv_biallelic[sv] += 1


for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
    census_chrom(chrom, by_chrom[chrom])
    print(f"  done {chrom}", flush=True)

chr13_bi = sum(v for k, v in sv_biallelic.items() if sv_chrom[k] == "chr13")
assert chr13_bi == 8, f"chr13 self-check FAILED: got {chr13_bi}, expected 8"
print("chr13 self-check OK (8 biallelic)")

out = f"{OUT_DIR}/step1_coverage.tsv"
with open(out, "w") as f:
    f.write("sv_locus\tchrom\tn_total\tn_biallelic\tn_recoverable\tn_multiallelic_only\tn_swapped\tpct_biallelic\n")
    for sv in sorted(sv_total, key=lambda s: (sv_chrom[s], s)):
        t = sv_total[sv]
        bi = sv_biallelic[sv]
        rec = sv_recoverable[sv]
        f.write(f"{sv}\t{sv_chrom[sv]}\t{t}\t{bi}\t{rec}\t{rec-bi}\t{sv_swapped[sv]}\t{bi/t:.3f}\n")

T = sum(sv_total.values())
BI = sum(sv_biallelic.values())
REC = sum(sv_recoverable.values())
print(f"\nSUMMARY over {len(sv_total)} SVs / {T} weight SNPs:")
print(f"  biallelic (scoreable now): {BI} ({BI/T:.1%})")
print(f"  recoverable if multiallelics split: {REC} ({REC/T:.1%})")
print(f"  multiallelic-only loss: {REC-BI} ({(REC-BI)/T:.1%})")
print(f"  absent entirely: {T-REC} ({(T-REC)/T:.1%})")
print(f"\nwrote {out}")
