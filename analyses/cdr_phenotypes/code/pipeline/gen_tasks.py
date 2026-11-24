import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
RD = (CENSV_ROOT + "/discovery"
      "/cdr_phenotypes")
OUT = RD + "/output"
V1D = (CENSV_ROOT + "/discovery"
       "/cdr_qtl/output")

POSITION = ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_pos_iqr", "cdr_pos_var"]
OTHER = ["n_cdrs", "cdrs_per_mb", "cdr_fraction", "median_cdr_width",
         "mean_5mC", "mean_5mC_covwt"]

F = pd.read_csv(OUT + "/frame_comparison.tsv", sep="\t")
D = pd.read_csv(OUT + "/phenotype_diff_v1_v2.tsv", sep="\t")
v1 = pd.read_csv(V1D + "/per_hap_chrom_phenotypes.tsv", sep="\t", engine="python")
v2 = pd.read_csv(OUT + "/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")

lost = [c for c in sorted(set(v1.chm13_chrom) | set(v2.chm13_chrom))
        if set(v1.loc[v1.chm13_chrom == c, "sample"]) != set(v2.loc[v2.chm13_chrom == c, "sample"])]

multi = set(F.loc[F.pct_multiblock > 0, "chrom"])
allchrom = sorted(set(F.chrom))
tasks = set()

for c in allchrom:
    for ph in POSITION:
        for fr in (["span", "block", "splice"] if c in multi else ["block"]):
            tasks.add((c, ph, fr))

for _, r in D[(D.verdict == "CHANGED") & (D.phenotype_v1.isin(OTHER))].iterrows():
    tasks.add((r.chrom, r.phenotype_v1, "block"))

for c in lost:
    for ph in POSITION + OTHER:
        tasks.add((c, ph, "block"))

tasks = sorted(tasks, key=lambda t: (t[2], t[1], t[0]))
with open(OUT + "/rescan_tasks.txt", "w") as f:
    for c, ph, fr in tasks:
        f.write("%s\t%s\t%s\n" % (c, ph, fr))

print("chromosomes with any multi-block array: %d  %s" % (len(multi), sorted(multi)))
print("chromosomes that lost an individual   : %s" % lost)
print("TOTAL TASKS: %d" % len(tasks))
print(pd.Series([t[2] for t in tasks]).value_counts().to_string())
