#!/usr/bin/env python
"""Every table the resolution-section figures read, assembled once from the primary outputs."""
import ast
import os
import re
import sys

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
       "qgap200_rgap200")
DEC = os.path.join(ROOT, "sv_extraction/decompose", SUF)
EN = os.path.join(ROOT, "results/elastic_net/decompose/cosine_nj/edge1000_" + SUF)
B8 = os.path.join(ROOT, "results/imputation/decompose/cosine_nj/best8_" + SUF)
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
OUT = os.path.join(WU, "output/resolution")
BED = os.path.join(ROOT, "resources/chm13v2.0_censat_v2.1.active.filtered.bed")
TRACKS = os.path.join(ROOT, "cenSV_web/etl")
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
HAND = {"HG002": "EUR"}

sys.path.insert(0, TRACKS)
import tracks                                                             # noqa: E402

os.makedirs(OUT, exist_ok=True)
NUM = re.compile(r"-?\d+")

MCOLS = ["r2", "cv_r2", "max_coef", "n_nonzero_coef", "gt_counts", "l1_ratio", "svfreq",
         "svpos", "svend", "svtype", "svlen", "alpha", "sv_id", "sv_category",
         "original_sv_id", "mapped_haplogroup", "r2_score", "sv_index",
         "selected_sv_snp_dist"]


def snp_stats(s):
    if not isinstance(s, str) or "[" not in s:
        return (np.nan,) * 5
    d = np.array([int(x) for x in NUM.findall(s.split("[", 1)[1])], dtype=float)
    if not len(d):
        return (np.nan,) * 5
    a = np.abs(d)
    return (float(a.min()), float(np.median(a)), float(a.max()),
            float(d.max() - d.min()), float((d < 0).mean()))


STRATS = ["edge500", "edge1000", "edge2000", "fksize100k", "fksize200k", "fksize500k",
          "midpoint500", "midpoint1000"]

met, dec, brn, hap = [], [], [], []
for c in CHROMS:
    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=MCOLS).assign(chrom=c)

    b8 = pd.read_csv("%s/%s/%s_best_strategy_selection.tsv" % (B8, c, c), sep="\t")
    assert len(b8) == len(m), ("%s: best-of-eight holds %d rows against sv_metrics' %d"
                               % (c, len(b8), len(m)))
    assert (b8.sv_index.values == np.arange(len(b8))).all(), (
        "%s: best-of-eight sv_index is not the row index, so the positional join is unsafe" % c)
    assert np.allclose(b8.cv_r2_edge1000.values, m.cv_r2.values, atol=1e-9, equal_nan=True), (
        "%s: best-of-eight's own edge1000 column disagrees with the edge1000 table it indexes"
        % c)
    m["cv_r2_edge1000"] = m.cv_r2
    m["cv_r2"] = b8.best_cv_r2.values
    m["best_strategy"] = b8.best_strategy_by_cvr2.values
    m["n_strategies"] = b8.n_strategies_available.values
    for s in STRATS:
        m["cv_" + s] = b8["cv_r2_" + s].values

    st = m.selected_sv_snp_dist.apply(snp_stats)
    m["snp_dist_min"] = [x[0] for x in st]
    m["snp_dist_med"] = [x[1] for x in st]
    m["snp_dist_max"] = [x[2] for x in st]
    m["snp_span"] = [x[3] for x in st]
    m["snp_frac_upstream"] = [x[4] for x in st]
    met.append(m.drop(columns=["selected_sv_snp_dist"]))

    d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, c, c), sep="\t",
                    usecols=["sv_id", "iteration", "type", "branch_id", "r2", "n_carriers"])
    dec.append(d.assign(chrom=c))

    b = pd.read_csv("%s/%s/%s_branches.tsv" % (DEC, c, c), sep="\t")
    brn.append(b.assign(chrom=c))

    h = pd.read_csv("%s/%s/%s_haplogroups.tsv" % (DEC, c, c), sep="\t")
    hap.append(h.assign(chrom=c))
    print("  %-6s metrics %6d  decomp %7d  branches %4d  haplotypes %4d"
          % (c, len(m), len(d), len(b), len(h)))

met = pd.concat(met, ignore_index=True)
dec = pd.concat(dec, ignore_index=True)
brn = pd.concat(brn, ignore_index=True)
hap = pd.concat(hap, ignore_index=True)

hap["hap_name"] = hap["sample"].astype(str)
hap["sample_id"] = hap.sample_id.astype(str)
sup, basis, byhand, panel_pop, why = tracks.resolve_ancestry(
    [s for s in sorted(hap.sample_id.unique()) if s not in HAND])
assert why is None, "ancestry join failed: %s" % why
sup.update(HAND)
hap["ancestry"] = hap.sample_id.map(sup)
assert hap.ancestry.notna().all(), "%d haplotypes carry no ancestry" % int(
    hap.ancestry.isna().sum())

gt = met.gt_counts.fillna("").str.extract(
    r"0:(?P<n0>\d+)/1:(?P<n1>\d+)/2:(?P<n2>\d+)").astype(float)
met = pd.concat([met, gt], axis=1)
met["n_called"] = met[["n0", "n1", "n2"]].sum(axis=1)

parent = met[met.sv_category == "original_decomposed"][
    ["chrom", "original_sv_id", "cv_r2", "svfreq", "n_nonzero_coef", "snp_dist_med"]].rename(
    columns={"cv_r2": "parent_cv_r2", "svfreq": "parent_svfreq",
             "n_nonzero_coef": "parent_n_snp", "snp_dist_med": "parent_snp_dist_med"})
assert not parent.duplicated(["chrom", "original_sv_id"]).any(), (
    "a resolved SV has two pre-resolution records; the join below would multiply rows")
met = met.merge(parent, on=["chrom", "original_sv_id"], how="left")

arr = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "name", "score", "strand", "ts", "te", "rgb"])
arr = arr.groupby("chrom").agg(array_start=("start", "min"), array_end=("end", "max"),
                               n_block=("start", "size"),
                               covered_bp=("end", lambda s: 0)).reset_index()
cov = pd.read_csv(BED, sep="\t", header=None, usecols=[0, 1, 2],
                  names=["chrom", "start", "end"])
arr["covered_bp"] = arr.chrom.map((cov.end - cov.start).groupby(cov.chrom).sum())
arr["array_len"] = arr.array_end - arr.array_start

met.to_csv(OUT + "/metrics_all.tsv", sep="\t", index=False)
dec.to_csv(OUT + "/decomp_all.tsv", sep="\t", index=False)
brn.to_csv(OUT + "/branches_all.tsv", sep="\t", index=False)
hap.to_csv(OUT + "/haplotypes_all.tsv", sep="\t", index=False)
arr.to_csv(OUT + "/arrays.tsv", sep="\t", index=False)

n = met.sv_category.value_counts()
print("\n--- record counts by class, against the manuscript ---")
for k, want in (("non_decomposed", 16609), ("original_decomposed", 2032),
                ("monophyletic", 1967), ("residual", 359)):
    print("  %-20s %6d   (manuscript %d)" % (k, n.get(k, 0), want))
    assert n.get(k, 0) == want, ("%s holds %d records, the manuscript reports %d"
                                 % (k, n.get(k, 0), want))
mean = met.groupby("sv_category").cv_r2.mean()
print("\n--- mean CV R2 (best of eight), against Table 1 ---")
for k, want in (("monophyletic", 0.407), ("original_decomposed", 0.371),
                ("non_decomposed", 0.065), ("residual", 0.078)):
    print("  %-20s %.3f   (Table 1 %.3f)" % (k, mean[k], want))
    assert abs(mean[k] - want) < 5e-4, (
        "%s has mean CV R2 %.4f here and %.3f in Table 1; the figures would contradict the "
        "table they sit beside" % (k, mean[k], want))
for k, want in (("monophyletic", 90.4), ("original_decomposed", 89.4),
                ("non_decomposed", 25.0), ("residual", 32.3)):
    got = 100 * (met[met.sv_category == k].cv_r2 > 0.1).mean()
    assert abs(got - want) < 0.05, ("%s: %.1f%% above R2 0.1 here, %.1f%% in Table 1"
                                    % (k, got, want))
print("  per-class %% above R2 0.1 and the eight-strategy join both reproduce Table 1")
print("\n--- which of the eight windows wins, by class ---")
print(pd.crosstab(met.sv_category, met.best_strategy).to_string())

print("\n--- best-branch R2 recorded for records the search turned down ---")
print(met.groupby("sv_category").r2_score.agg(["size", "count", "median"]).to_string())
print("\n--- decomposition steps ---")
print(dec.type.value_counts().to_string())
print("\n--- haplotypes by ancestry ---")
print(hap.drop_duplicates("sample_id").ancestry.value_counts().to_string())
print("\nwrote five tables to %s" % OUT)
