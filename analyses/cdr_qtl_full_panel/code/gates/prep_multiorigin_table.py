#!/usr/bin/env python
"""One tidy row per monophyletic component, carrying everything a supplementary display needs."""
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
DEC = ROOT + "/sv_extraction/decompose/" + SUF
EN = ROOT + "/results/elastic_net/decompose/cosine_nj/edge1000_" + SUF
WU = ROOT + "/analyses/cdr_qtl_full_panel"
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]

logs, var, met = [], [], []
for c in CHROMS:
    d = pd.read_csv("%s/%s/%s_decomposition.tsv" % (DEC, c, c), sep="\t")
    d = d[d.type == "monophyletic"].assign(chrom=c)
    logs.append(d[["chrom", "sv_id", "iteration", "branch_id", "r2", "n_carriers"]])
    v = pd.read_csv("%s/%s/%s_sv_variants.tsv" % (DEC, c, c), sep="\t", low_memory=False,
                    usecols=["chrom", "id", "pos", "svtype", "svlen", "ac", "af",
                             "original_sv_id", "mapped_haplogroup", "sv_category"])
    var.append(v[v.sv_category == "monophyletic"])
    m = pd.read_csv("%s/%s/sv_metrics.tsv" % (EN, c), sep="\t", low_memory=False,
                    usecols=["sv_id", "sv_category", "cv_r2", "r2", "n_nonzero_coef", "svfreq"])
    met.append(m[m.sv_category == "monophyletic"].assign(chrom=c))
logs = pd.concat(logs, ignore_index=True)
var = pd.concat(var, ignore_index=True)
met = pd.concat(met, ignore_index=True).rename(columns={"r2": "in_sample_r2"})

t = var.merge(logs, left_on=["chrom", "original_sv_id", "mapped_haplogroup"],
              right_on=["chrom", "sv_id", "branch_id"], how="left",
              suffixes=("", "_log")).drop(columns=["sv_id"])
print("components: %d | with an iteration matched: %d" % (len(t), int(t.iteration.notna().sum())))

t = t.merge(met[["chrom", "sv_id", "cv_r2", "in_sample_r2", "n_nonzero_coef"]],
            left_on=["chrom", "id"], right_on=["chrom", "sv_id"], how="left").drop(
    columns=["sv_id"])
t["modelled"] = t.cv_r2.notna()
t["origin"] = t.iteration.map(lambda i: "first" if i == 0 else ("later" if i >= 1 else "?"))

n_or = t.groupby(["chrom", "original_sv_id"]).size().rename("sv_n_origins")
t = t.merge(n_or, on=["chrom", "original_sv_id"], how="left")
t["sv_multi_origin"] = t.sv_n_origins >= 2

OUT = WU + "/output/monophyletic_components.tsv"
t.to_csv(OUT, sep="\t", index=False)
print("wrote %s  (%d rows, %d columns)" % (OUT, len(t), t.shape[1]))

print("\n--- sanity, against this script's own last run ---")
print("  components                    %5d   (expect 2247)" % len(t))
print("  modelled, edge1000            %5d   (expect 1954)" % int(t.modelled.sum()))
print("  SVs with >=1 component        %5d   (expect 2103)"
      % t.groupby(["chrom", "original_sv_id"]).ngroups)
print("  SVs with >1 component         %5d   (expect  135)" % int(n_or[n_or >= 2].size))
print("  later origins, joined         %5d   (expect  141)" % int((t.origin == "later").sum()))
print("  later origins modelled        %5d   (expect   30)"
      % int(((t.origin == "later") & t.modelled).sum()))
print("  components with no iteration  %5d   (expect   19, all chr10)"
      % int(t.iteration.isna().sum()))

print("\n--- the question no display can be drawn without ---")
print("cross-validated R2 of MODELLED components, first vs later origin:")
for nm, s in t[t.modelled].groupby("origin").cv_r2:
    print("  %-6s n=%4d  mean %.3f  median %.3f  %%>0.1 %.1f  %%>0.3 %.1f"
          % (nm, len(s), s.mean(), s.median(), 100 * (s > .1).mean(), 100 * (s > .3).mean()))
print("\nallele count, modelled vs not, by origin:")
print(t.groupby(["origin", "modelled"]).ac.agg(["size", "min", "median", "max"]).to_string())
print("\nbranch-match R2 (the tree fit, not imputation), by origin:")
print(t.groupby("origin").r2.agg(["size", "mean", "median"]).to_string())
