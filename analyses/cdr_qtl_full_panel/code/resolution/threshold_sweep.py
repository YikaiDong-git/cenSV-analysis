#!/usr/bin/env python
"""What the resolution returns when its three thresholds are moved."""
import argparse
import importlib.util
import os

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
VCF = os.path.join(ROOT, "variant_calling/svim_asm", "cohort_svim_asm_%s.vcf.gz" % SUF)
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
OUT = os.path.join(WU, "output/resolution/sweep")
HC = os.path.join(ROOT, "workflow/scripts/haplogroup_clustering.py")

R2S = [0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
FLOORS = [5, 8, 10, 15, 20, 30]
SECOND = [3, 5, 10]
PROD = (0.3, 10, 5)

spec = importlib.util.spec_from_file_location("hc", HC)
hc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hc)

ap = argparse.ArgumentParser()
ap.add_argument("--chrom", required=True)
a = ap.parse_args()
C = a.chrom
os.makedirs(OUT, exist_ok=True)

matrix, variants, samples = hc.parse_vcf_region(VCF, C, 2, 1.0, None)
hg = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, C, C), sep="\t")
TAB = hg["sample"].astype(str).tolist()
at = {s: i for i, s in enumerate(samples)}
order = TAB + [s for s in samples if s not in set(TAB)]
matrix = matrix[[at[s] for s in order]]
samples = order

res = hc.perform_clustering(matrix, samples, (2, 30), 5, "cosine_nj")
assert res["samples"] == TAB, "%s: haplotype order does not match the haplogroup table" % C
bv = hc.extract_branch_vectors(res["nj_tree"], res["samples"], min_tips=5)
X = res["matrix"]
print("%s: %d haplotypes x %d SVs, %d branches" % (C, X.shape[0], X.shape[1], len(bv)))


def run(r2t, floor, second):
    recs, _, _, _ = hc.decompose_svs(X, bv, variants, res["samples"], r2_threshold=r2t,
                                     min_carriers=second, min_haplogroup_size=floor)
    m = [r for r in recs if r["type"] == "monophyletic"]
    per = pd.Series([r["sv_id"] for r in m]).value_counts()
    return m, per


shipped = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, C, C), sep="\t")
shipped = shipped[shipped.type == "monophyletic"]
shipped_car = shipped.carrier_indices.fillna("").apply(
    lambda s: frozenset(int(i) for i in s.split(",") if s))
m0, per0 = run(*PROD)
assert len(m0) == len(shipped), (
    "%s: at the production setting the sweep yields %d components against the shipped %d"
    % (C, len(m0), len(shipped)))
by = {}
for r in m0:
    by.setdefault(r["sv_id"], []).append(frozenset(r["carrier_indices"]))
same = sum(1 for sid, car in zip(shipped.sv_id, shipped_car) if car in by.get(sid, []))
got_r2 = np.sort([r["r2"] for r in m0])
want_r2 = np.sort(shipped.r2.values)
assert np.allclose(got_r2, want_r2, atol=1e-9), (
    "%s: at the production setting the component R2 values differ from the shipped ones "
    "(largest gap %.3g)" % (C, float(np.abs(got_r2 - want_r2).max())))
print("  production cell: %d components with the shipped R2 values, %d of %d with identical "
      "carriers" % (len(m0), same, len(shipped)))

rows = []
for r2t in R2S:
    for floor in FLOORS:
        for second in SECOND:
            if second > floor:
                continue
            m, per = run(r2t, floor, second)
            r2v = np.array([r["r2"] for r in m]) if m else np.array([0.0])
            nc = np.array([r["n_carriers"] for r in m]) if m else np.array([0])
            rows.append(dict(
                chrom=C, r2_threshold=r2t, first_floor=floor, second_floor=second,
                components=len(m), svs_resolved=int(per.size),
                multi_origin=int((per >= 2).sum()),
                later_origins=int(sum(1 for r in m if r["iteration"] >= 1)),
                components_over_10=int((nc >= 10).sum()),
                median_branch_r2=float(np.median(r2v)),
                median_carriers=float(np.median(nc)),
                total_carriers=int(nc.sum()),
                is_production=(r2t, floor, second) == PROD))
    print("  R2 >= %.2f done" % r2t)

pd.DataFrame(rows).to_csv("%s/%s_sweep.tsv" % (OUT, C), sep="\t", index=False)
print("wrote %s/%s_sweep.tsv  (%d cells)" % (OUT, C, len(rows)))
