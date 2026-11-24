#!/usr/bin/env python
"""Candidate multi-origin SVs for the ten tree examples, each ranked on a different property."""
import os

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
O = WU + "/output"

geo = pd.read_csv(O + "/multiorigin_tree_geometry.tsv", sep="\t")
det = pd.read_csv(O + "/multiorigin_components_detail.tsv", sep="\t")
nul = pd.read_csv(O + "/multiorigin_null.tsv", sep="\t")
met = pd.read_csv(O + "/resolution/metrics_all.tsv", sep="\t", low_memory=False)

mono = met[met.sv_category == "monophyletic"][
    ["chrom", "sv_id", "original_sv_id", "mapped_haplogroup", "cv_r2", "svfreq",
     "n_nonzero_coef", "best_strategy"]]
det["haplogroups"] = det.haplogroups.fillna("").astype(str)
det["ancestry"] = det.ancestry.fillna("").astype(str)
det = det.merge(mono, left_on=["chrom", "original_sv_id", "branch"],
                right_on=["chrom", "original_sv_id", "mapped_haplogroup"], how="left")
det["modelled"] = det.cv_r2.notna()

g = det.groupby(["chrom", "original_sv_id"]).agg(
    pos=("pos", "first"), svtype=("svtype", "first"), svlen=("svlen", "first"),
    origins=("branch", "size"), modelled=("modelled", "sum"),
    carriers=("n_carriers", "sum"), min_car=("n_carriers", "min"),
    max_car=("n_carriers", "max"), min_r2=("branch_r2", "min"),
    max_r2=("branch_r2", "max"), best_cv=("cv_r2", "max"),
    anc=("n_ancestries", "min"), anc_max=("n_ancestries", "max"),
    hgs=("haplogroups", lambda s: " | ".join(s)),
    ancs=("ancestry", lambda s: " | ".join(s))).reset_index()
g = g.merge(geo[["chrom", "original_sv_id", "residual_carriers", "max_sep"]],
            on=["chrom", "original_sv_id"], how="left")
g = g.merge(nul.groupby(["chrom", "original_sv_id"]).percentile.min().rename("pct_min"),
            on=["chrom", "original_sv_id"], how="left")
g = g.merge(nul.groupby(["chrom", "original_sv_id"]).percentile.max().rename("pct_max"),
            on=["chrom", "original_sv_id"], how="left")

SHOW = ["chrom", "pos", "svtype", "svlen", "origins", "modelled", "carriers", "min_car",
        "min_r2", "max_r2", "best_cv", "anc", "anc_max", "residual_carriers", "pct_min"]


def top(name, frame, n=6):
    print("\n=== %s ===" % name)
    print(frame[SHOW].head(n).to_string(index=False))


top("most origins", g.sort_values(["origins", "carriers"], ascending=False))
top("every origin modelled", g[g.modelled == g.origins].sort_values("carriers",
                                                                   ascending=False))
top("most carriers", g.sort_values("carriers", ascending=False))
top("both branch matches near perfect", g[g.min_r2 >= 0.8].sort_values("min_r2",
                                                                      ascending=False))
top("later origin far below the ten-carrier floor", g.sort_values("min_car"))
top("origins in the fewest ancestries (lineage-restricted)",
    g[g.origins == 2].sort_values(["anc_max", "carriers"], ascending=[True, False]))
top("origins spanning the most ancestries", g.sort_values("anc", ascending=False))
top("clades closest together against the null", g.sort_values("pct_min"))
top("clades farthest apart against the null", g.sort_values("pct_min", ascending=False))
top("largest variant", g.sort_values("svlen", ascending=False))
top("residual carriers left unexplained", g.sort_values("residual_carriers",
                                                        ascending=False))

print("\n=== the chromosome-9 hypomethylated-fraction leads (branch B561) ===")
b = det[det.branch == "B561"]
print(b[["chrom", "pos", "svtype", "svlen", "branch", "n_carriers", "branch_r2", "cv_r2",
         "ancestry", "haplogroups"]].to_string(index=False))
print("\nand their other origins:")
print(det[det.original_sv_id.isin(b.original_sv_id)][
    ["chrom", "pos", "branch", "iteration", "n_carriers", "branch_r2", "cv_r2", "ancestry"]
].to_string(index=False))

print("\n=== SV length spectrum of multi-origin SVs ===")
print((g.svlen // 100 * 100).value_counts().head(12).to_string())
print("\nlengths within 15 bp of a multiple of 171: %d of %d"
      % (int(((g.svlen % 171).apply(lambda x: min(x, 171 - x)) <= 15).sum()), len(g)))
g.to_csv(O + "/multiorigin_summary.tsv", sep="\t", index=False)
print("\nwrote %s" % (O + "/multiorigin_summary.tsv"))
