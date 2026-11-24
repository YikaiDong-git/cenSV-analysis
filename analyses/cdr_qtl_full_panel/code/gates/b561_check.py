#!/usr/bin/env python
"""Three chr9 cdr_fraction leads all lead on branch B561."""
import itertools

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = ROOT + "/analyses/cdr_qtl_full_panel"

comp = pd.read_csv(WU + "/output/multi_origin_svs.tsv", sep="\t")
b = comp[(comp.chrom == "chr9") & (comp.branch_id == "B561")]
print("chr9 components landing on branch B561: %d" % len(b))
sets = {}
for r in b.itertuples():
    sets[r.sv_id] = set(map(int, str(r.carrier_indices).split(",")))
    print("  %-3d carriers  R2 %.2f  %s" % (r.n_carriers, r.r2, r.sv_id[:60] + "..."))

print("\npairwise carrier-set overlap (Jaccard):")
for a, c in itertools.combinations(sets, 2):
    j = len(sets[a] & sets[c]) / len(sets[a] | sets[c])
    print("  %.3f  (%d shared of %d union)" % (j, len(sets[a] & sets[c]), len(sets[a] | sets[c])))

G = np.load(WU + "/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(WU + "/geno/full_variants.tsv", sep="\t")
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}
gate = pd.read_csv(WU + "/output/gates/active_len_gate_all_hits.tsv", sep="\t")
leads = gate[(gate.chrom == "chr9") & (gate.phenotype == "cdr_fraction")]
print("\nchr9 cdr_fraction leads in the gate table: %d" % len(leads))
names, mats = [], []
for r in leads.itertuples():
    d = G[loc2row[r.sv_locus], :].copy()
    d[d == -1] = np.nan
    names.append("%s %s" % (r.sv, r.tier[:3]))
    mats.append(d)
M = np.vstack(mats)
print("\npairwise dosage correlation between the leads the scan tested:")
print("%-26s %s" % ("", "  ".join("%8s" % n.split()[0][:8] for n in names)))
for i, n in enumerate(names):
    cells = []
    for j in range(len(names)):
        ok = np.isfinite(M[i]) & np.isfinite(M[j])
        cells.append("%8.2f" % np.corrcoef(M[i][ok], M[j][ok])[0, 1])
    print("%-26s %s" % (n, "  ".join(cells)))
