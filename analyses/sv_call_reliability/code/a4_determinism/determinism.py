#!/usr/bin/env python3
import csv, json, sys, os, random
from collections import defaultdict
import numpy as np
import os as _censv_os
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

J = (CENSV_COLLAB_ROOT + "/exploration/simulate_svim_asm_accuracy/"
     "results/cohort_frequency_v1/joint_analysis")
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10**9)
random.seed(1)
NH = 256

rows = {}
with open(J + "/merged_rows.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        r.pop("ref"); r.pop("alt")
        rows[(r["chrom"], int(r["row_index"]))] = r
car = defaultdict(set); car_status = defaultdict(list)
with open(J + "/source_membership.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["chrom"], int(r["row_index"]))
        car[k].add(int(r["haplotype_index"])); car_status[k].append(r["status"])
ev = defaultdict(set); evinfo = {}
with open(J + "/event_carriers.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        k = (r["chrom"], r["event_id"])
        ev[k].add(int(r["haplotype_index"]))
        evinfo[k] = (r["recipe"], r["true_type"], int(r["size"]), int(r["true_ac"]))

blocks = lambda s: {h // 32 for h in s}
ac32 = [s for k, s in ev.items() if evinfo[k][3] == 32]
one_block32 = sum(len(blocks(s)) == 1 and len(s) == 32 for s in ac32)
block_note = dict(events_ac32=len(ac32), events_ac32_exactly_one_block=one_block32)


def jac(a, b):
    return len(a & b) / len(a | b) if a or b else 0.0


def block_union_j(s):
    bs = set(b for b in range(NH // 32) if len(s & set(range(32 * b, 32 * b + 32))) > 16)
    u = set(h for b in bs for h in range(32 * b, 32 * b + 32))
    return jac(s, u), sorted(bs)


out = []
for k, r in rows.items():
    if r["passes_ac_af"] != "True" or r["dedup_retained"] != "True":
        continue
    maf = float(r["maf"])
    st = car_status[k]
    if not st or any(x != "fp" for x in st):
        continue
    grp = "common" if maf >= 0.05 else ("rare" if maf < 0.01 else None)
    if grp is None:
        continue
    C = car[k]
    best = max(((jac(C, T), e) for e, T in ev.items() if e[0] == k[0]), default=(0.0, None))
    bj, bset = block_union_j(C)
    null = [jac(set(random.sample(range(NH), len(C))), ev[best[1]]) for _ in range(200)] if best[1] else [0]
    out.append(dict(chrom=k[0], row_index=k[1], pos=r["pos"], svtype=r["svtype"], svlen=r["svlen"], maf=maf,
                    group=grp, n_carriers=len(C), best_event=best[1][1] if best[1] else "",
                    best_event_recipe=evinfo[best[1]][0] if best[1] else "",
                    best_event_true_ac=evinfo[best[1]][3] if best[1] else "",
                    jaccard_best_event=round(best[0], 3),
                    null_jaccard_mean=round(float(np.mean(null)), 3), null_jaccard_p99=round(float(np.percentile(null, 99)), 3),
                    jaccard_block_union=round(bj, 3), blocks=";".join(map(str, bset))))

with open(OUT + "/false_rows_carrier_alignment.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
summ = {"founder_block_check": block_note}
for g in ("common", "rare"):
    S = [o for o in out if o["group"] == g]
    if not S:
        continue
    jb = np.array([o["jaccard_best_event"] for o in S]); jn = np.array([o["null_jaccard_p99"] for o in S])
    jbk = np.array([o["jaccard_block_union"] for o in S])
    summ[g] = dict(n_all_false_rows=len(S),
                   median_carriers=float(np.median([o["n_carriers"] for o in S])),
                   median_jaccard_best_event=float(np.median(jb)),
                   frac_jaccard_ge_0_9=round(float(np.mean(jb >= 0.9)), 3),
                   frac_jaccard_ge_0_5=round(float(np.mean(jb >= 0.5)), 3),
                   frac_above_null_p99=round(float(np.mean(jb > jn)), 3),
                   frac_block_union_ge_0_9=round(float(np.mean(jbk >= 0.9)), 3),
                   by_chrom={c: sum(o["chrom"] == c for o in S) for c in sorted(set(o["chrom"] for o in S))})
json.dump(summ, open(OUT + "/summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
