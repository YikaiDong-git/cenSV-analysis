#!/usr/bin/env python
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
a5 = pd.read_csv(f"{WU}/output/all_tests_1e5.tsv.gz", sep="\t")
a7 = pd.read_csv(f"{WU}/output/all_tests_1e7.tsv.gz", sep="\t")
done = set(map(tuple, a7[["chrom", "phenotype"]].drop_duplicates().values))

THR = 1e-4
g = (a5.groupby(["chrom", "phenotype"])
     .agg(min_p=("p_emp_perm", "min"), n_loci=("pos", "size")).reset_index())
g["escalated_r1"] = [(c, p) in done for c, p in zip(g.chrom, g.phenotype)]
todo = g[(~g.escalated_r1) & (g.min_p <= THR)].sort_values("min_p")
print(f"round-1 cells: {g.escalated_r1.sum()};  candidate round-2 cells (min p <= {THR:g}): {len(todo)}")
print(todo.to_string(index=False))
print(f"\ntotal extra tests at 1e7: {int(todo.n_loci.sum())}")
print(f"largest single cell: {int(todo.n_loci.max()) if len(todo) else 0} loci")

out = f"{WU}/code/scan/escalate_tasks_round2.txt"
with open(out, "w") as fh:
    for r in todo.itertuples():
        fh.write(f"{r.chrom}\t{r.phenotype}\n")
print(f"\nwrote {out} ({len(todo)} tasks)")

r1 = g[g.escalated_r1]
print(f"\nround-1 cells ranged {int(r1.n_loci.min())}-{int(r1.n_loci.max())} loci, "
      f"{int(r1.n_loci.sum())} tests total")
