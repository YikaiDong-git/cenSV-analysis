#!/usr/bin/env python3
import glob
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

M = CENSV_COLLAB_ROOT + "/exploration/meQTL"
SUB = (CENSV_ROOT + "/analyses/"
       "cdr_qtl_full_panel")

TARGETS = [("chr12 DEL-333 (SV-specific meQTL)", "chr12", 37183837, 333),
           ("chr8 DEL-113 (strongest cis meQTL)", "chr8", 44300221, 113)]

print("=" * 72)
print("A. IS THE SV IN THE HPRC ASSEMBLY CATALOG (this unit's cohort B)?")
print("=" * 72)
var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
gs = np.array(open(SUB + "/geno/full_samples.txt").read().split())
for name, ch, pos, svlen in TARGETS:
    sub = var[var.chrom == ch]
    exact = sub[(sub.pos == pos)]
    near = sub[(sub.pos - pos).abs() <= 500]
    print("\n" + name)
    print("  catalog SVs on " + ch + ": " + str(len(sub))
          + "   exact pos match: " + str(len(exact))
          + "   within 500 bp: " + str(len(near)))
    if len(near):
        n = near.assign(d=(near.pos - pos).abs())
        cols = ["sv_locus", "pos", "end", "svtype", "svlen", "af", "geno_group",
                "status", "d"]
        print(n.sort_values("d")[cols].head(8).to_string(index=False))
        gg = sorted(near.geno_group.unique())
        print("  distinct genotype-groups: " + str(len(gg)))
        for g in gg[:4]:
            ri = int(near[near.geno_group == g].index[0])
            v = G[ri, :]
            print("    group " + str(g) + ": carriers " + str(int((v > 0).sum()))
                  + "/" + str(int((v >= 0).sum())) + "  het " + str(int((v == 1).sum()))
                  + "  hom " + str(int((v == 2).sum())))
    else:
        print("  -> NOT in the HPRC catalog at this position")

print("\n" + "=" * 72)
print("B. IS THE SV IN COHORT A'S RESULTS -- BROAD PC30 vs CORRECTED MAC8/PC10?")
print("=" * 72)

SCANS = {
    "MAC8 / PC10 (corrected, post-audit)":
        M + "/results/meqtl_cis/chm13_observed_sv_robust_mac8_pc10/summary/association_results.tsv.gz",
}
for p in sorted(glob.glob(M + "/results/meqtl_cis/chm13_observed_sv_same_chrom/summary/**/*assoc*",
                          recursive=True)):
    SCANS["broad PC30: " + os.path.basename(os.path.dirname(p)) + "/" + os.path.basename(p)] = p

for label, path in SCANS.items():
    print("\n--- " + label)
    if not os.path.exists(path):
        print("    (missing) " + path)
        continue
    try:
        d = pd.read_csv(path, sep="\t")
    except Exception as e:
        print("    unreadable: " + str(e))
        continue
    poscol = "sv_pos" if "sv_pos" in d.columns else ("source_pos" if "source_pos" in d.columns else None)
    chrcol = "sv_chrom" if "sv_chrom" in d.columns else ("source_chrom" if "source_chrom" in d.columns else None)
    if poscol is None or chrcol is None:
        print("    columns: " + str(list(d.columns)[:12]))
        continue
    print("    rows " + str(len(d)) + "   distinct SVs " + str(d[poscol].nunique()))
    for name, ch, pos, svlen in TARGETS:
        hit = d[(d[chrcol] == ch) & ((d[poscol] - pos).abs() <= 500)]
        print("    " + name.split(" (")[0].ljust(14) + " rows: " + str(len(hit)))
        if len(hit):
            qc = [c for c in hit.columns if "bh_q" in c and "primary" in c and "_t" in c]
            bc = [c for c in hit.columns if c in ("rint_beta", "expected_rint_beta")]
            pc = [c for c in hit.columns if c in ("rint_hc3_p_t", "expected_rint_p")]
            show = ["feature_id"] + bc + pc + qc[:1]
            show = [c for c in show if c in hit.columns]
            if qc:
                sig = hit[hit[qc[0]] < 0.01]
                print("      at q<0.01: " + str(len(sig)))
            print(hit[show].head(4).to_string(index=False))
