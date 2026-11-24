#!/usr/bin/env python
import glob
import numpy as np, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")
T = f"{WU}/output/trans"
bad = []


def chk(label, ok, detail=""):
    print(f"  {'OK  ' if ok else 'FAIL'} {label}{('   ' + detail) if detail else ''}")
    if not ok:
        bad.append(label)


print("=== 1. two independent scanners must agree on shared cells ===")
F = pd.read_csv(f"{T}/focused_trans_results.tsv", sep="\t")
G = pd.read_csv(f"{T}/gw_results.tsv", sep="\t")
cis = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{T}/gw/*/*_cis.tsv")
                 if pd.read_csv(f, sep="\t").shape[0]], ignore_index=True)
GW = pd.concat([G, cis], ignore_index=True)
k = ["sv_locus", "phenotype", "target_chrom"]
J = F[k + ["t", "n", "relation"]].merge(GW[k + ["t", "n"]], on=k, suffixes=("_foc", "_gw"))
print(f"  shared cells found: {len(J)}")
if len(J):
    dt = (J.t_foc - J.t_gw).abs()
    dn = (J.n_foc - J.n_gw).abs()
    chk("t agrees to 1e-6", bool((dt < 1e-6).all()), f"max |dt| = {dt.max():.3g}")
    chk("n identical", bool((dn == 0).all()), f"max |dn| = {int(dn.max())}")
    print(f"  (worst row: {J.loc[dt.idxmax(), 'sv_locus']} {J.loc[dt.idxmax(),'phenotype']} "
          f"{J.loc[dt.idxmax(),'target_chrom']})")

print("\n=== 2. the chr18 association, as it will be quoted ===")
h = F[(F.sv_locus == "chr8_sv1_INS_non_best8_midpoint1000")
      & (F.phenotype == "cdr_pos_offset_deepest") & (F.target_chrom == "chr18")]
r = h.iloc[0]
print(f"  focused: n={r.n} t={r.t:.4f} p_perm_1e6={r.p_perm_1e6:.3g} "
      f"BH={r.bh_q_trans_escalated:.4f} FWER(960)={r.fwer_global:.4f}")
X = pd.read_csv(f"{T}/gw_context_for_8j.tsv", sep="\t")
chk("all 8 focused hits demoted to gw FWER >= 0.99", bool((X.gw_fwer >= 0.99).all()),
    f"min gw_fwer = {X.gw_fwer.min():.3f}")

gmax = np.max(np.vstack([np.load(f) for f in sorted(glob.glob(f"{T}/gw/*/*_maxT.npy"))]), axis=0)
chk("chr18 |t| is below the 5th percentile of the genome-wide null max",
    abs(r.t) < np.quantile(gmax, 0.05),
    f"|t|={abs(r.t):.3f} vs null-max p05={np.quantile(gmax,0.05):.3f}")

print("\n=== 3. genome-wide calibration and result ===")
C = pd.read_csv(f"{T}/gw_cells.tsv", sep="\t")
chk("192 cells", len(C) == 192)
chk("lambda median < 1.1", C.lambda_trans.median() < 1.1, f"{C.lambda_trans.median():.4f}")
chk("candidates <= expectation", C.n_cand_p01.sum() <= C.exp_cand_p01.sum(),
    f"{int(C.n_cand_p01.sum()):,} vs {C.exp_cand_p01.sum():,.0f} "
    f"({C.n_cand_p01.sum()/C.exp_cand_p01.sum():.3f}x)")
chk("0 at FWER<0.05", int((G.fwer < 0.05).sum()) == 0)
chk("0 at BH q<0.05", int((G.bh_q < 0.05).sum()) == 0)
print(f"  family size: {int(C.n_trans.sum()):,} trans tests, {int(C.n_cis.sum()):,} cis")
top = G.nsmallest(1, "p_param").iloc[0]
print(f"  strongest trans: {top.chrom}:{top.pos:,} {top.svtype} -> {top.target_chrom} "
      f"{top.phenotype}  t={top.t:.3f} p={top.p_param:.3g} FWER={top.fwer:.3f} BH={top.bh_q:.3f}")
print(f"  cells with lambda > 1.2: {int((C.lambda_trans > 1.2).sum())} -> "
      f"{sorted(zip(C[C.lambda_trans>1.2].phenotype, C[C.lambda_trans>1.2].target_chrom))}")

print("\n=== 4. cis positive control ===")
c8 = cis[(cis.chrom == "chr8") & cis.pos.between(44_200_000, 44_500_000)]
chk("known chr8 leads are the strongest cis signal",
    float(cis.abs_t.max()) == float(c8.abs_t.max()),
    f"max |t| overall {cis.abs_t.max():.3f}, chr8 window {c8.abs_t.max():.3f}")

print("\n=== 5. liftover fix, as it will be quoted ===")
M = pd.read_csv(f"{WU}/output/mqtl/validation/fix_impact_matrix.tsv", sep="\t")
S = pd.read_csv(f"{WU}/output/mqtl/validation/fix_impact_scan.tsv", sep="\t")
RT = pd.read_csv(f"{WU}/output/mqtl/validation/roundtrip_mirror_comparison.tsv", sep="\t")
CD = pd.read_csv(f"{WU}/output/mqtl/validation/candidates_pre_vs_post.tsv", sep="\t")
print(M.to_string(index=False))
print(S.to_string(index=False))
print(RT[["variant", "pct_exact", "pct_within_1kb", "pct_width_exact"]].to_string(index=False))
pre = set(zip(CD[CD.variant == "pre-fix"].chrom, CD[CD.variant == "pre-fix"].pos))
post = set(zip(CD[CD.variant == "corrected"].chrom, CD[CD.variant == "corrected"].pos))
chk("candidate loci identical pre/post", pre == post, f"{sorted(post)}")
chk("max |dt| from the fix <= 0.26", float(S.max_abs_dt.max()) <= 0.26,
    f"{S.max_abs_dt.max():.3f}")

print("\n=== 6. chr18 anomaly ===")
A = pd.read_csv(f"{T}/chr18_anomaly_by_chrom.tsv", sep="\t")
P = pd.read_csv(f"{T}/chr18_phenotype_availability.tsv", sep="\t")
r18 = A[A.chrom == "chr18"].iloc[0]
chk("chr18 has the lowest alignment coverage", r18.pct_blocks == A.pct_blocks.min(),
    f"{r18.pct_blocks:.2f}% vs next {sorted(A.pct_blocks)[1]:.2f}%")
chk("chr18's gap is alignment, not methylation", r18.pct_aligned_but_no_meth < 1.0,
    f"aligned_but_no_meth = {r18.pct_aligned_but_no_meth:.4f}%")
p18 = P[P.chrom == "chr18"].iloc[0]
chk("chr18 individuals with a call == those with a complete phenotype set",
    p18.indiv_any_row == p18.indiv_deepest == p18.indiv_mean5mC,
    f"{p18.indiv_any_row:.0f} rows, {p18.indiv_deepest:.0f} deepest, {p18.indiv_mean5mC:.0f} 5mC")

print(f"\n{'ALL CHECKS PASSED' if not bad else 'FAILED: ' + str(bad)}")
