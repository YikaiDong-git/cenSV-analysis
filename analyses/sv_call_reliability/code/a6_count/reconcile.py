#!/usr/bin/env python3
"""A6: reconcile the manuscript's 188,589 SVs with the 188,591 rows in Josh's composition summary."""
import csv, glob, os, sys, json
from collections import Counter, defaultdict
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

P = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
CUR = P + "/sv_extraction/decompose/" + SUF
ARCH = P + "/discovery/2026-09-11_regen_chr15_chr18/archive_pre_regen_2026-09-11"
COMP = (CENSV_COLLAB_ROOT + "/exploration/simulate_svim_asm_accuracy/"
        "results/panel_repeat_composition.tsv")
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10**9)
CH = ["chr%d" % i for i in range(1, 23)] + ["chrX"]


def table_keys(path):
    with open(path, errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ip, ie, it, il = (hdr.index(x) for x in ("pos", "end", "svtype", "svlen"))
        mx = max(ip, ie, it, il)
        keys = Counter()
        for line in fh:
            f = line.split("\t", mx + 1)
            if len(f) <= mx:
                continue
            keys[(f[ip], f[ie], f[it], f[il])] += 1
    return keys


snap = defaultdict(dict)
for r in csv.DictReader(open(P + "/paper/figures/data/catalog/per_chrom_svtype.tsv"), delimiter="\t"):
    snap[r["chrom"]][r["svtype"]] = int(r["n"])
cur = {}
for c in CH:
    cur[c] = table_keys("%s/%s/%s_variants.tsv" % (CUR, c, c))
res = dict(per_chrom={}, diffs={})
for c in CH:
    ct = Counter()
    for (p_, e_, t_, l_), n in cur[c].items():
        ct[t_] += n
    st = snap[c]
    res["per_chrom"][c] = dict(current=dict(ct), snapshot=st,
                               differs=any(ct.get(t, 0) != st.get(t, 0) for t in set(ct) | set(st)))
tot_cur = Counter(); tot_snap = Counter()
for c in CH:
    tot_cur.update(res["per_chrom"][c]["current"]); tot_snap.update(res["per_chrom"][c]["snapshot"])
res["totals"] = dict(current=dict(tot_cur), current_total=sum(tot_cur.values()),
                     snapshot=dict(tot_snap), snapshot_total=sum(tot_snap.values()))

archived = {}
for c in [c for c in CH if res["per_chrom"][c]["differs"]]:
    cand = glob.glob("%s/**/%s_variants.tsv" % (ARCH, c), recursive=True)
    res["diffs"][c] = dict(archive_files=cand)
    if cand:
        old = table_keys(cand[0])
        archived[c] = old
        added = [list(k) + [n] for k, n in (cur[c] - old).items()]
        removed = [list(k) + [n] for k, n in (old - cur[c]).items()]
        res["diffs"][c].update(n_old=sum(old.values()), n_new=sum(cur[c].values()),
                               added=added, removed=removed)

drop = set()
for c, d in res["diffs"].items():
    for p_, e_, t_, l_, n in d.get("added", []):
        for L in {l_, l_.lstrip("-"), str(abs(int(float(l_)))) if l_.lstrip("-").replace(".", "").isdigit() else l_}:
            drop.add("%s:%s:%s:%s:%s" % (c, p_, e_, t_, L))
n_rows = 0; comp = Counter(); lens = []; ndrop = 0
for r in csv.DictReader(open(COMP), delimiter="\t"):
    if r["catalog_row_id"] in drop:
        ndrop += 1; continue
    n_rows += 1
    comp["all_alpha_rich" if r["all_alleles_alpha_rich"] == "True" else "other"] += 1
    if int(r["n_TE_rich_alleles"]) > 0:
        comp["TE_rich_present"] += 1
    if int(r["n_TE_positive_alleles"]) > 0:
        comp["TE_any"] += 1
    if r["annotation_status"] == "unsupported_type":
        comp["symbolic_unclassified"] += 1
    if r["all_alleles_alpha_rich"] == "True" and r["svtype"] in ("INS", "DEL") and abs(int(r["svlen"])) >= 50:
        lens.append(abs(int(r["svlen"])))
lens = np.array(lens)
k = np.maximum(1, np.rint(lens / 171))
hit10 = float(np.mean(np.abs(lens - k * 171) <= 10))
n_added = sum(a[-1] for d in res["diffs"].values() for a in d.get("added", []))
n_removed = sum(a[-1] for d in res["diffs"].values() for a in d.get("removed", []))
res["manuscript_catalog"] = dict(rows=n_rows, rows_dropped=ndrop, rows_added_since_snapshot=n_added,
                                 rows_removed_since_snapshot=n_removed,
                                 drop_keys_matched_all_added=(ndrop == n_added), composition=dict(comp),
                                 alpha_rich_pct=round(100 * comp["all_alpha_rich"] / n_rows, 2),
                                 TE_rich_pct=round(100 * comp["TE_rich_present"] / n_rows, 3),
                                 n_alpha_indels_ge50=int(len(lens)), frac_171_pm10=round(hit10, 4))
json.dump(res, open(OUT + "/reconcile.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "per_chrom"}, indent=1))
