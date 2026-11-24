#!/usr/bin/env python
"""Are any of the 38 gate leads carried by an SV that has more than one origin?"""
import re

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"

comp = pd.read_csv(WU + "/output/multi_origin_svs.tsv", sep="\t")
multi = comp.groupby(["chrom", "sv_id"]).size().rename("n_origins").reset_index()
key = dict(zip(zip(multi.chrom, multi.sv_id), multi.n_origins))

vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
vt["orig_id"] = vt.sv_id.astype(str).str.replace(r"_(B\d+|residual|original)$", "", regex=True)
vt["branch"] = vt.sv_id.astype(str).str.extract(r"_(B\d+|residual|original)$")[0]
loc2 = vt.set_index("sv_locus")

gate = pd.read_csv(WU + "/output/gates/active_len_gate_all_hits.tsv", sep="\t")
al = pd.read_csv(WU + "/output/gates/full_hits_gate_1e7.tsv", sep="\t")
gate = gate.merge(al[["chrom", "sv_locus", "phenotype", "boundary_verdict"]],
                  on=["chrom", "sv_locus", "phenotype"], how="left", suffixes=("", "_g"))

rows = []
for r in gate.itertuples():
    v = loc2.loc[r.sv_locus]
    n = key.get((r.chrom, v.orig_id), 1)
    rows.append(dict(chrom=r.chrom, sv=r.sv, phenotype=r.phenotype, tier=r.tier,
                     sv_locus=r.sv_locus, orig_id=v.orig_id, branch=v.branch,
                     n_origins=n, prior=r.prior_verdict,
                     length_confounded=bool(r.p_active_len < 0.05)))
h = pd.DataFrame(rows)

print("=== the %d gate leads, by how many origins their SV has ===" % len(h))
print(h.groupby("n_origins").size().rename("leads").to_string())
mo = h[h.n_origins >= 2]
print("\nleads carried by a MULTI-ORIGIN SV: %d, on %d chromosomes, %d distinct SVs"
      % (len(mo), mo.chrom.nunique() if len(mo) else 0, mo.orig_id.nunique() if len(mo) else 0))
if not len(mo):
    print("  none -- every hit comes from an SV the tree explained with a single branch")
else:
    print()
    print("%-6s %-22s %-24s %-12s %-9s %8s  %s"
          % ("chrom", "SV", "phenotype", "tier", "component", "origins", "length-conf"))
    for r in mo.sort_values(["chrom", "phenotype"]).itertuples():
        print("%-6s %-22s %-24s %-12s %-9s %8d  %s"
              % (r.chrom, r.sv, r.phenotype, r.tier, r.branch, r.n_origins,
                 "YES" if r.length_confounded else "no"))

    print("\n--- the component structure of each such SV ---")
    for oid in mo.orig_id.unique():
        c = comp[comp.sv_id == oid].sort_values("iteration")
        lead_branch = set(mo[mo.orig_id == oid].branch)
        print("\n  %s  (%s)" % (c.chrom.iloc[0], oid if len(oid) < 70 else oid[:67] + "..."))
        for x in c.itertuples():
            mark = "  <-- the lead" if ("B%s" % x.branch_id) in lead_branch or \
                                       str(x.branch_id) in lead_branch else ""
            print("      iter %d  branch %-6s  R2 %.2f  %3d carriers%s"
                  % (x.iteration, x.branch_id, x.r2, x.n_carriers, mark))

print("\n=== how often is a lead a second-or-later component? ===")
print(h.branch.value_counts().to_string())
