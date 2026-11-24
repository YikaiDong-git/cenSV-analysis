#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Emit the per-chromosome resolution table as LaTeX, for sec/supp_tables.tex."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resfig import CH, all_components, resolved_svs  # noqa: E402

ac = all_components()
svs = resolved_svs()

g = svs.assign(m=svs.n_origins >= 2).groupby("chrom").agg(svs=("m", "size"), multi=("m", "sum"))
comp = ac.groupby("chrom").size().rename("components")
later = ac[ac.origin == "later"].groupby("chrom").size().rename("later")
lmod = ac[(ac.origin == "later") & ac.modelled].groupby("chrom").size().rename("lmod")
t = pd.concat([g, comp, later, lmod], axis=1).reindex(CH).fillna(0).astype(int)


def n(v):
    return "{:,}".format(int(v)).replace(",", "{,}")


print(r"\begin{table}[htbp]")
print(r"\centering")
print(r"\caption{\textbf{Haplogroup resolution by chromosome} (one row per centromere). "
      r"``Resolved SVs'' are the SVs whose carriers were explained by at least one branch of "
      r"that chromosome's haplotype tree; ``components'' counts the monophyletic components "
      r"they were separated into, so a chromosome with more components than resolved SVs "
      r"carries SVs of more than one origin. ``Later origins'' are the components found after "
      r"the first branch match was masked, and ``imputed'' how many of those reached an "
      r"elastic-net model. The last column is the share of that chromosome's resolved SVs "
      r"with more than one origin, the quantity drawn in Supplementary "
      r"Figure~\ref{fig:supp_resolution}e.}")
print(r"\label{tab:resolution_chrom}")
print(r"{\footnotesize")
print(r"\begin{tabular}{lrrrrrr}")
print(r"\toprule")
print(r"Chromosome & Resolved SVs & $>1$ origin & Components & Later origins & Imputed & "
      r"$>1$ origin (\%) \\")
print(r"\midrule")
for c, r in zip(CH, t.itertuples()):
    pct = ("%.1f" % (100.0 * r.multi / r.svs)) if r.svs else "---"
    print(r"%s & %s & %d & %s & %d & %d & %s \\"
          % (c.replace("chr", ""), n(r.svs), r.multi, n(r.components), r.later, r.lmod, pct))
print(r"\midrule")
print(r"All & %s & %d & %s & %d & %d & %.1f \\"
      % (n(t.svs.sum()), t.multi.sum(), n(t.components.sum()), t.later.sum(), t.lmod.sum(),
         100.0 * t.multi.sum() / t.svs.sum()))
print(r"\bottomrule")
print(r"\end{tabular}")
print(r"}")
print(r"\end{table}")
