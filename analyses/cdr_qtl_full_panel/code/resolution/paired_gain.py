#!/usr/bin/env python
"""What resolution is worth, measured four ways on the same variants."""
import os

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = os.path.join(ROOT, "analyses/cdr_qtl_full_panel")
A = os.path.join(WU, "output/resolution")

met = pd.read_csv(A + "/metrics_all.tsv", sep="\t", low_memory=False)
c = met[met.sv_category == "monophyletic"].copy()
p = met[met.sv_category == "original_decomposed"].copy()

pair = c[c.parent_cv_r2.notna()].copy()
pair["d"] = pair.cv_r2 - pair.parent_cv_r2
print("components with a parent record: %d of %d" % (len(pair), len(c)))
print("parents: %d; components per parent: mean %.3f"
      % (len(p), len(c) / float(len(p))))


def describe(name, v):
    v = np.asarray(v, dtype=float)
    print("  %-46s n=%5d  mean %+.4f  median %+.4f  q1 %+.4f  q3 %+.4f  above 0 %.1f%%"
          % (name, len(v), v.mean(), np.median(v), np.quantile(v, .25), np.quantile(v, .75),
             100 * (v > 0).mean()))


print("\n--- 1. one row per component ---")
describe("component minus its own parent", pair.d)

print("\n--- 2. averaged within parent first ---")
g = pair.groupby(["chrom", "original_sv_id"]).agg(comp=("cv_r2", "mean"),
                                                  par=("parent_cv_r2", "first"),
                                                  k=("cv_r2", "size"))
describe("mean of a parent's components minus the parent", g.comp - g.par)
print("     parents splitting into one component:  %d" % int((g.k == 1).sum()))
print("     parents splitting into two or more:    %d" % int((g.k >= 2).sum()))

print("\n--- 3. the best component of each parent ---")
b = pair.groupby(["chrom", "original_sv_id"]).agg(best=("cv_r2", "max"),
                                                  par=("parent_cv_r2", "first"))
describe("best component minus the parent", b.best - b.par)
print("     parents whose best component clears R2 > 0.1 where the parent did not: %d"
      % int(((b.best > .1) & (b.par <= .1)).sum()))
print("     parents that cleared it and whose best component does not:             %d"
      % int(((b.best <= .1) & (b.par > .1)).sum()))

print("\n--- 4. the class means of Table 1, and what separates them ---")
print("  monophyletic components  n=%5d  mean %.4f" % (len(c), c.cv_r2.mean()))
print("  pre-resolution parents   n=%5d  mean %.4f" % (len(p), p.cv_r2.mean()))
print("  difference of the class means                     %+.4f"
      % (c.cv_r2.mean() - p.cv_r2.mean()))
print("  same difference, taken pairwise (framing 1)       %+.4f" % pair.d.mean())
print("  same difference, one row per parent (framing 2)   %+.4f" % (g.comp - g.par).mean())
one = g[g.k == 1]
multi = g[g.k >= 2]
print("  pairwise, parents giving ONE component  n=%4d    %+.4f"
      % (len(one), (one.comp - one.par).mean()))
print("  pairwise, parents giving TWO or more    n=%4d    %+.4f"
      % (len(multi), (multi.comp - multi.par).mean()))

print("\n--- and the comparison Figure 2F draws: components against UNRESOLVED SVs ---")
u = met[met.sv_category == "non_decomposed"]
for lo, hi, lab in ((0, .05, "<0.05"), (.05, .1, "0.05-0.10"), (.1, .2, "0.10-0.20"),
                    (.2, .3, "0.20-0.30"), (.3, 1.01, ">=0.30")):
    a1 = c[(c.svfreq >= lo) & (c.svfreq < hi)].cv_r2
    a2 = u[(u.svfreq >= lo) & (u.svfreq < hi)].cv_r2
    a3 = p[(p.svfreq >= lo) & (p.svfreq < hi)].cv_r2
    print("  %-10s component %.3f (n=%4d)   unresolved %.3f (n=%5d)   parent %.3f (n=%4d)"
          % (lab, a1.mean(), len(a1), a2.mean(), len(a2), a3.mean(), len(a3)))

print("\n--- how often resolution changed the genotype at all ---")
met["carriers"] = met.n1.fillna(0) + met.n2.fillna(0)
pcar = met[met.sv_category == "original_decomposed"][
    ["chrom", "original_sv_id", "carriers"]].rename(columns={"carriers": "par_car"})
pair = pair.merge(met[["chrom", "sv_id", "carriers"]], on=["chrom", "sv_id"], how="left")
pair = pair.merge(pcar, on=["chrom", "original_sv_id"], how="left")
same_geno = pair[pair.carriers >= pair.par_car]
print("  components holding every one of their parent's carriers: %d of %d (%.1f%%)"
      % (len(same_geno), len(pair), 100 * len(same_geno) / len(pair)))
describe("  their paired difference", same_geno.d)
print("     of those, CV R2 identical to the parent to 1e-6: %d"
      % int((same_geno.d.abs() < 1e-6).sum()))
moved = pair[pair.carriers < pair.par_car]
describe("components that did lose carriers", moved.d)
if len(moved):
    print("  of those, carriers removed: median %d of %d (%.0f%%)"
          % (int((moved.par_car - moved.carriers).median()), int(moved.par_car.median()),
             100 * float(((moved.par_car - moved.carriers) / moved.par_car).median())))

nopar = set(zip(p.chrom, p.original_sv_id)) - set(zip(pair.chrom, pair.original_sv_id))
print("\n  parents in Table 1 whose components never reached a model: %d of %d"
      % (len(nopar), len(p)))

print("\n--- how much of the class-mean gap is the counting rule ---")
w = pair.groupby(["chrom", "original_sv_id"]).size().rename("k")
pw = pair.merge(w, on=["chrom", "original_sv_id"])
print("  mean over components, weighting each parent once: %.4f"
      % (pw.assign(w=1 / pw.k).pipe(lambda f: (f.cv_r2 * f.w).sum() / f.w.sum())))
print("  mean over components as Table 1 counts them:      %.4f" % c.cv_r2.mean())
