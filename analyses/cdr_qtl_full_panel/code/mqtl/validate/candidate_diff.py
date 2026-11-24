#!/usr/bin/env python
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
S = f"{WU}/output/mqtl/scans"
rows = []
for W in (10, 50, 100):
    for lab, d in [("pre-fix", "bam_mirrorbug"), ("corrected", "bam")]:
        f = f"{S}/{d}/escalated_{W}kb_bam_10e6.tsv"
        try:
            E = pd.read_csv(f, sep="\t")
        except FileNotFoundError:
            continue
        n_tested = len(pd.read_csv(f"{S}/{d}/scan_{W}kb_bam.tsv", sep="\t"))
        alpha = 0.05 / n_tested
        pcol = "p_meth_1e6" if "p_meth_1e6" in E.columns else "p_meth"
        tcol = "t_meth_1e6" if "t_meth_1e6" in E.columns else "t_meth"
        ccol = "p_cov_1e6" if "p_cov_1e6" in E.columns else "p_cov"
        P = E[(E[pcol] <= alpha) & (E[ccol] > 0.05)]
        for r in P.itertuples():
            rows.append(dict(window_kb=W, variant=lab, alpha=alpha, sv_locus=r.sv_locus,
                             chrom=r.chrom, pos=r.pos, svtype=r.svtype, af=r.af, n=r.n,
                             t=getattr(r, tcol), p=getattr(r, pcol), p_cov=getattr(r, ccol)))
D = pd.DataFrame(rows)
D.to_csv(f"{WU}/output/mqtl/validation/candidates_pre_vs_post.tsv", sep="\t", index=False,
         float_format="%.6g")

print("COVERAGE-CLEAN CANDIDATES PASSING WITHIN-TIER alpha, by window and substrate")
for W in (10, 50, 100):
    g = D[D.window_kb == W]
    if not len(g):
        print(f"\n  {W} kb: none"); continue
    print(f"\n  {W} kb  (alpha = {g.alpha.iloc[0]:.3g})")
    for lab in ("pre-fix", "corrected"):
        h = g[g.variant == lab]
        loci = sorted(set(zip(h.chrom, h.pos)))
        print(f"    {lab:10s} rows {len(h):2d}  distinct loci {len(loci)}: "
              + "; ".join(f"{c}:{p:,}" for c, p in loci))

print("\nPER-LOCUS, pre-fix vs corrected (t and p at 1e7):")
key = ["window_kb", "chrom", "pos", "sv_locus"]
a = D[D.variant == "pre-fix"].set_index(key)[["t", "p", "n"]]
b = D[D.variant == "corrected"].set_index(key)[["t", "p", "n"]]
J = a.join(b, how="outer", lsuffix="_pre", rsuffix="_post").reset_index()
J["dt"] = J.t_post - J.t_pre
J["status"] = np.where(J.t_pre.isna(), "NEW", np.where(J.t_post.isna(), "LOST", "kept"))
print(J.sort_values(["window_kb", "p_post"]).to_string(index=False))
print(f"\nkept {int((J.status=='kept').sum())}   lost {int((J.status=='LOST').sum())}   "
      f"new {int((J.status=='NEW').sum())}   max |dt| among kept "
      f"{J.loc[J.status=='kept','dt'].abs().max():.4g}")
print("\ndone")
