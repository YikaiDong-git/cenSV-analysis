#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Per-origin evidence tier for every multi-origin row the browser flags."""
import collections
import csv
import io
import json
import os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
DISC = os.path.join(P1, "analyses/cdr_qtl_full_panel/output")
DATA = os.path.join(P1, "cenSV_web/web/public/data")
TRUE = ("True", "true", "1")
TAB = "\t"


def rd(p):
    with io.open(p, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=TAB))


def num(x):
    x = (x or "").strip()
    return None if x in ("", "NA", "nan", "None") else float(x)


mc = {}
for r in rd(os.path.join(DISC, "monophyletic_components.tsv")):
    mc[(r["chrom"], int(float(r["pos"])), r["mapped_haplogroup"])] = r

cat = json.load(io.open(os.path.join(DATA, "catalog.json"), encoding="utf-8"))
ci = {k: i for i, k in enumerate(cat["columns"])}
G = lambda k, i: cat["data"][ci[k]][i]
pos_of = {}
cat_rows = set()
for i in range(cat["rows"]):
    pos_of[(G("chrom", i), G("sv_index", i))] = int(G("pos", i))
    cat_rows.add((G("chrom", i), int(G("pos", i)), G("mapped_haplogroup", i)))

idx = json.load(io.open(os.path.join(DATA, "hap_index.json"), encoding="utf-8"))
out, demo = {}, []
for chrom in sorted(idx["multi"]):
    blk = json.load(io.open(os.path.join(DATA, "hap_%s.json" % chrom), encoding="utf-8"))
    tips = blk["tree"]["n_tip"]
    for sv in idx["multi"][chrom]:
        origins = blk["origins"][sv]
        pos = pos_of.get((chrom, int(sv)))
        recs = []
        for branch, r2, _packed in origins:
            m = mc.get((chrom, pos, branch))
            modelled = bool(m) and m["modelled"] in TRUE
            recs.append({
                "branch": branch,
                "tier": "imputed" if modelled else "assembly only",
                "cv_r2": num(m["cv_r2"]) if m else None,
                "n_carriers": int(float(m["n_carriers"])) if m and m["n_carriers"] else None,
                "in_catalog": (chrom, pos, branch) in cat_rows,
            })
        out.setdefault(chrom, {})[sv] = recs
        if any(r["tier"] == "assembly only" for r in recs[1:]):
            demo.append((chrom, sv, pos, tips, [(r["branch"], r["tier"], r["n_carriers"])
                                                for r in recs]))

p = os.path.join(DISC, "origin_tiers.json")
with io.open(p, "w", encoding="utf-8") as fh:
    fh.write(json.dumps(out, separators=(",", ":"), sort_keys=True))
print("wrote %s  (%d chromosomes, %d rows)" % (p, len(out), sum(len(v) for v in out.values())))

flat = [r for v in out.values() for rs in v.values() for r in rs]
print("origins: %d | imputed %d | assembly only %d"
      % (len(flat), sum(1 for r in flat if r["tier"] == "imputed"),
         sum(1 for r in flat if r["tier"] == "assembly only")))
later = [rs[1:] for v in out.values() for rs in v.values()]
later = [r for g in later for r in g]
print("later origins: %d | imputed %d | assembly only %d | absent from catalog.json %d"
      % (len(later), sum(1 for r in later if r["tier"] == "imputed"),
         sum(1 for r in later if r["tier"] == "assembly only"),
         sum(1 for r in later if not r["in_catalog"])))
print()
print("rows where a later origin is assembly only: %d" % len(demo))
print("candidates (largest trees first):")
for c, sv, pos, tips, rs in sorted(demo, key=lambda x: -x[3])[:12]:
    print("  %-6s sv_index %-6s pos %-10s tips %-4d  %s"
          % (c, sv, pos, tips, "  ".join("%s/%s/%s" % r for r in rs)))
