#!/usr/bin/env python3
import csv, sys
from collections import defaultdict

MANIFEST, PIP, FG_DESC, NL_DESC, OUT = sys.argv[1:6]


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


rank = {}
for row in csv.DictReader(open(MANIFEST), delimiter="\t"):
    rank[(row["study"], row["pheno_code"])] = int(row["pheno_rank"])

sig, rob = defaultdict(set), defaultdict(set)
for row in csv.DictReader(open(PIP), delimiter="\t"):
    key = (row["study"], row["pheno"])
    sig[key].add(row["sv_locus"])
    if max(f(row["impute_pip"]), f(row["matched_pip"])) >= 0.5:
        rob[key].add(row["sv_locus"])

fg_name = {r["phenocode"]: (r.get("phenotype") or "").strip()
           for r in csv.DictReader(open(FG_DESC), delimiter="\t")}
nl_name = {r["phenotype_code"]: (r.get("description") or "").strip()
           for r in csv.DictReader(open(NL_DESC), delimiter="\t")}

rows = []
for (study, pheno), rk in rank.items():
    name = (fg_name.get(pheno) if study == "finngen" else nl_name.get(pheno)) or pheno
    rows.append((pheno, name, study, rk,
                 len(sig.get((study, pheno), ())), len(rob.get((study, pheno), ()))))
rows.sort(key=lambda x: (x[2], x[3], x[0]))

with open(OUT, "w") as o:
    o.write("trait_code\ttrait_name\tsource_code\tpheno_rank\tn_sv_pip_ge_0p01\tn_sv_pip_ge_0p50\n")
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")

miss = [r[0] for r in rows if r[1] == r[0]]
by = defaultdict(int)
for r in rows:
    by[r[2]] += 1
sys.stderr.write(f"traits: {len(rows)} {dict(by)}; missing trait_name (fell back to code): {len(miss)}\n")
if miss:
    sys.stderr.write("  missing: " + ", ".join(miss[:15]) + ("..." if len(miss) > 15 else "") + "\n")
