#!/usr/bin/env python
"""decomp_dupe_audit.py — antinomy #4: how many INDEPENDENT genomic events are behind the handoff SVs?"""
import csv, re
from collections import defaultdict
from itertools import combinations
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P   = CENSV_WORKSPACE + ""
CT  = f"{P}/Project1_Centromere/analyses/disease_replication"
CAT = f"{P}/Project1_Centromere/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"

rows = {r["sv_locus"]: r for r in csv.DictReader(open(CAT), delimiter="\t")}
TOK = re.compile(r"svim_asm\.(?:INS|DEL)\.\d+")

def read_col1(path):
    return [r[0] for r in csv.reader(open(path), delimiter="\t")][1:]

def cluster(sv_locus):
    r = rows.get(sv_locus)
    if r is None: return None
    return frozenset(TOK.findall(r["sv_id"]))

def audit(name, svs, meta_ntraits=None):
    print(f"\n{'='*84}\n{name}: {len(svs)} handoff sv_loci\n{'='*84}")
    miss = [s for s in svs if s not in rows]
    if miss: print(f"  !! {len(miss)} NOT found in catalog: {miss}")
    clusters = {s: cluster(s) for s in svs if s in rows}
    groups = defaultdict(list)
    for s, c in clusters.items(): groups[c].append(s)
    dupes = {k: v for k, v in groups.items() if len(v) > 1}
    print(f"  distinct genomic events (exact svim-cluster): {len(groups)}  (from {len(clusters)} SVs)")
    print(f"  EXACT-duplicate event groups (>1 handoff SV = one event double-counted): {len(dupes)}")
    for k, v in sorted(dupes.items(), key=lambda kv: -len(kv[1])):
        cats = [rows[s]["sv_category"] for s in v]
        nt = [meta_ntraits.get(s, "?") for s in v] if meta_ntraits else None
        print(f"    * {len(v)} SVs, {len(k)} svim calls, chrom {rows[v[0]]['primary_chrom']}: "
              f"{[s.split('_')[0]+'_'+s.split('_')[1] for s in v]} cats={cats}" + (f" n_traits={nt}" if nt else ""))
    partials = []
    items = list(clusters.items())
    for (s1, c1), (s2, c2) in combinations(items, 2):
        if c1 and c2 and c1 != c2 and (c1 & c2):
            j = len(c1 & c2) / len(c1 | c2)
            partials.append((j, s1.split('_')[0]+'_'+s1.split('_')[1], s2.split('_')[0]+'_'+s2.split('_')[1], len(c1 & c2)))
    partials.sort(reverse=True)
    print(f"  PARTIAL-overlap SV pairs (share svim calls, not identical): {len(partials)}")
    for j, a, b, n in partials[:12]:
        print(f"    ~ {a} & {b}: Jaccard={j:.2f} ({n} shared calls)")
    return clusters

meta222 = {r["sv_locus"]: r["n_traits_hit"] for r in csv.DictReader(open(f"{CT}/phewide_222_handoff/sv_metadata_222.tsv"), delimiter="\t")}
sv33 = read_col1(f"{CT}/phewide_222_handoff/sv_metadata_222.tsv")
cl33 = audit("222-SET (33 SVs)", sv33, meta222)

sv50 = read_col1(f"{CT}/aou_handoff/sv_metadata.tsv")
cl50 = audit("HANDOFF-50", sv50)

pairs = list(csv.DictReader(open(f"{CT}/phewide_222_handoff/pairs_222.tsv"), delimiter="\t"))
ev_of = {s: cl33.get(s) for s in sv33}
uniq_pair = set()
for r in pairs:
    ev = ev_of.get(r["sv_locus"])
    uniq_pair.add((ev, r["trait"]))
print(f"\n{'='*84}\n222 PAIRS inflation\n{'='*84}")
print(f"  raw pairs = {len(pairs)}")
print(f"  unique (EVENT, trait) pairs after event-dedup = {len(uniq_pair)}")
print(f"  => decomposition double-counting inflates the pair count by {len(pairs)-len(uniq_pair)} ({(1-len(uniq_pair)/len(pairs))*100:.0f}%)")
print("\nDONE_DECOMP_AUDIT")
