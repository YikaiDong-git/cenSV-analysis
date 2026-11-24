#!/usr/bin/env python3
"""A0 scope probe: can the jwein22 MAC8 observed-SV meQTL discoveries be re-tested in the HPRC assembly cohort?"""
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

PROJ = CENSV_ROOT + ""
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
AUD = PROJ + "/analyses/cdr_calls"
REPORT = (CENSV_COLLAB_ROOT + "/exploration/meQTL/"
          "reports/chm13_observed_sv_meqtl_mac8_files")
OUT = PROJ + "/analyses/meqtl_replication/output/arms/scope"

TOLERANCES = [0, 10, 50, 200]
LIFT_HALFWIDTH = 100000


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


os.makedirs(OUT, exist_ok=True)

pairs = pd.read_csv(REPORT + "/all_discovery_conditioning_consensus.tsv", sep="\t")
print("report pairs: " + str(len(pairs)))
pairs["source_type"] = pairs.sv_id.str.split("_").str[3]

src = (pairs.groupby(["sv_id", "source_chrom", "source_pos", "source_end", "source_type"])
       .size().reset_index(name="n_pairs"))
print("report distinct source SVs: " + str(len(src)))

var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t")
G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
samples = [l.strip() for l in open(SUB + "/geno/full_samples.txt") if l.strip()]
assert G.shape == (len(var), len(samples)), (G.shape, len(var), len(samples))
print("HPRC substrate: " + str(G.shape[0]) + " loci x " + str(G.shape[1]) + " samples")

coh = pd.read_csv(AUD + "/output/cohort_corrected_final.tsv", sep="\t")
coh = coh[coh.chemistry.isin(["R941", "R1041"])]
chem = (coh.assign(is_r1041=(coh.chemistry == "R1041").astype(int))
        .groupby("sample").is_r1041.max())
meth_samples = set(chem.index)
r1041_samples = set(chem[chem == 1].index)
print("HPRC methylation cohort: " + str(len(meth_samples)) + " samples ("
      + str(len(r1041_samples)) + " R10.4.1, "
      + str(len(meth_samples) - len(r1041_samples)) + " R9.4.1)")

geno_has_meth = np.array([s in meth_samples for s in samples])
geno_is_r1041 = np.array([s in r1041_samples for s in samples])
print("genotype x methylation overlap: " + str(int(geno_has_meth.sum()))
      + " samples (" + str(int(geno_is_r1041.sum())) + " of them R10.4.1)")

q1_rows = []
match_map = []
for tol in TOLERANCES:
    n_matched = 0
    n_groups = []
    for _, s in src.iterrows():
        hit = var[(var.chrom == s.source_chrom)
                  & (var.svtype == s.source_type)
                  & ((var.pos - s.source_pos).abs() <= tol)
                  & ((var.end - s.source_end).abs() <= tol)]
        if len(hit):
            n_matched += 1
            n_groups.append(hit.geno_group.nunique())
            if tol == 0:
                for _, h in hit.iterrows():
                    match_map.append(dict(sv_id=s.sv_id, n_pairs=s.n_pairs,
                                          sv_locus=h.sv_locus, geno_group=h.geno_group,
                                          hprc_af=h.af, hprc_cv_r2=h.cv_r2,
                                          chrom=h.chrom, pos=h.pos, end=h.end,
                                          svtype=h.svtype, svlen=h.svlen))
    q1_rows.append(dict(tolerance_bp=tol, report_svs=len(src), matched=n_matched,
                        matched_frac=n_matched / float(len(src)),
                        median_geno_groups_per_sv=(float(np.median(n_groups))
                                                   if n_groups else 0.0)))
q1 = pd.DataFrame(q1_rows)
print("\nQ1 SV matching:")
print(q1.to_string(index=False))
write_atomic(q1, OUT + "/q1_sv_match_by_tolerance.tsv")

MM = pd.DataFrame(match_map)
if len(MM) == 0:
    print("\nNO EXACT SV MATCHES -- stop here and inspect coordinate conventions "
          "(the report's source_pos/source_end vs the HPRC catalog pos/end).")
    write_atomic(pd.DataFrame([dict(exact_matches=0)]), OUT + "/scope_summary.tsv")
    raise SystemExit(0)

q2_rows = []
for gg, blk in MM.groupby("geno_group"):
    i = var.index[var.geno_group == gg][0]
    g = G[i, :].astype(float)
    g_meth = g[geno_has_meth]
    g_r1041 = g[geno_is_r1041]
    q2_rows.append(dict(
        geno_group=gg,
        report_sv_ids=";".join(sorted(set(blk.sv_id))),
        n_report_pairs=int(blk.n_pairs.sum()),
        chrom=blk.chrom.iloc[0], pos=int(blk.pos.iloc[0]),
        svtype=blk.svtype.iloc[0], svlen=int(blk.svlen.iloc[0]),
        hprc_af=float(blk.hprc_af.iloc[0]),
        n_meth=int(geno_has_meth.sum()),
        n_carrier_meth=int((g_meth > 0).sum()),
        n_het_meth=int((g_meth == 1).sum()),
        n_homalt_meth=int((g_meth == 2).sum()),
        n_r1041=int(geno_is_r1041.sum()),
        n_carrier_r1041=int((g_r1041 > 0).sum()),
        n_het_r1041=int((g_r1041 == 1).sum()),
    ))
q2 = pd.DataFrame(q2_rows).sort_values(["chrom", "pos", "geno_group"]).reset_index(drop=True)
print("\nQ2 matched genotype-groups: " + str(len(q2)))
print("  testable at n_carrier_meth >= 5 : " + str(int((q2.n_carrier_meth >= 5).sum())))
print("  ASM-testable at n_het_meth  >= 10: " + str(int((q2.n_het_meth >= 10).sum())))
print("  ASM-testable at n_het_r1041 >=  5: " + str(int((q2.n_het_r1041 >= 5).sum())))
write_atomic(q2, OUT + "/q2_matched_sv_power.tsv")

loci_by_chrom = {}
for c, b in var.groupby("chrom"):
    loci_by_chrom[c] = np.sort(b.pos.values)

tgt = (pairs.groupby(["feature_id", "target_chrom", "target_start"])
       .size().reset_index(name="n_pairs"))


def covered(chrom, p):
    if chrom not in loci_by_chrom:
        return False
    a = loci_by_chrom[chrom]
    j = int(np.searchsorted(a, p))
    for k in (j - 1, j):
        if 0 <= k < len(a) and abs(int(a[k]) - int(p)) <= LIFT_HALFWIDTH:
            return True
    return False


tgt["already_lifted"] = [covered(c, p) for c, p in
                         zip(tgt.target_chrom, tgt.target_start)]
q3 = (tgt.groupby("already_lifted")
      .agg(n_cpgs=("feature_id", "size"), n_pairs=("n_pairs", "sum")).reset_index())
print("\nQ3 target CpGs: " + str(len(tgt)) + " distinct")
print(q3.to_string(index=False))
write_atomic(tgt.sort_values(["target_chrom", "target_start"]), OUT + "/q3_target_cpgs.tsv")

matched_ids = set(MM.sv_id)
pairs["sv_matched_hprc"] = pairs.sv_id.isin(matched_ids)
pairs["dist_class"] = pd.cut(pairs.source_target_distance_bp.fillna(-1),
                             [-2, -0.5, 1e5, 1e6, 5e6, 2.5e7, np.inf],
                             labels=["trans", "<=100kb", "0.1-1Mb", "1-5Mb",
                                     "5-25Mb", ">25Mb"])
q4 = (pairs.groupby(["dist_class", "sv_matched_hprc"], observed=True)
      .agg(n_pairs=("pair_id", "size"),
           neg_beta_frac=("expected_rint_beta", lambda x: float((x < 0).mean())),
           median_beta=("expected_rint_beta", "median")).reset_index())
print("\nQ4 composition of the replication set:")
print(q4.to_string(index=False))
write_atomic(q4, OUT + "/q4_replication_set_composition.tsv")

summary = pd.DataFrame([dict(
    report_pairs=len(pairs),
    report_source_svs=len(src),
    hprc_loci=len(var),
    hprc_geno_samples=len(samples),
    hprc_meth_samples=len(meth_samples),
    hprc_r1041_samples=len(r1041_samples),
    geno_meth_overlap=int(geno_has_meth.sum()),
    geno_r1041_overlap=int(geno_is_r1041.sum()),
    svs_matched_exact=int(q1.loc[q1.tolerance_bp == 0, "matched"].iloc[0]),
    svs_matched_200bp=int(q1.loc[q1.tolerance_bp == 200, "matched"].iloc[0]),
    matched_geno_groups=len(q2),
    groups_carrier_ge5=int((q2.n_carrier_meth >= 5).sum()),
    groups_het_ge10=int((q2.n_het_meth >= 10).sum()),
    groups_het_r1041_ge5=int((q2.n_het_r1041 >= 5).sum()),
    target_cpgs=len(tgt),
    target_cpgs_already_lifted=int(tgt.already_lifted.sum()),
    pairs_with_matched_sv=int(pairs.sv_matched_hprc.sum()),
)])
write_atomic(summary, OUT + "/scope_summary.tsv")
print("\n===== SCOPE SUMMARY =====")
print(summary.T.to_string(header=False))
