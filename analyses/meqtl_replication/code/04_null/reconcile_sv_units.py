#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = CENSV_ROOT + ""
UNIT = PROJ + "/analyses/meqtl_replication"
SUB = PROJ + "/analyses/cdr_qtl_full_panel"
OUT = UNIT + "/output/null_design"


def write_atomic(df, path):
    tmp = path + ".tmp." + str(os.getpid())
    df.to_csv(tmp, sep="\t", index=False, float_format="%.10g")
    os.replace(tmp, path)
    print("  wrote " + path + "  (" + str(len(df)) + " rows)")


os.makedirs(OUT, exist_ok=True)

pairs = pd.read_csv(UNIT + "/output/targets/targets_pair.tsv", sep="\t")
svgg = pd.read_csv(UNIT + "/output/targets/targets_sv_genogroup.tsv", sep="\t")
var = pd.read_csv(SUB + "/geno/full_variants.tsv", sep="\t", engine="python")
G = np.load(SUB + "/geno/full_genotypes.npz")["genotypes"]
samples = np.array(open(SUB + "/geno/full_samples.txt").read().split())
assert G.shape == (len(var), len(samples)), (G.shape, len(var), len(samples))

tags = [t.strip() for t in open(UNIT + "/output/a1_lift/all_tags.txt") if t.strip()]
meth_indiv = sorted({t.rsplit("_", 1)[0] for t in tags})
print("methylation cohort: " + str(len(tags)) + " haplotypes, "
      + str(len(meth_indiv)) + " individuals")
in_geno = [s for s in meth_indiv if s in set(samples)]
print("of those, present in the genotype matrix: " + str(len(in_geno)))
midx = np.array([int(np.where(samples == s)[0][0]) for s in in_geno])

print("\nreport SVs: " + str(pairs.sv_id.nunique())
      + "   genotype-groups: " + str(svgg.geno_group.nunique())
      + "   (sv, group) rows: " + str(len(svgg)))

rows = []
for sv_id, blk in pairs.groupby("sv_id"):
    s = blk.iloc[0]
    hit = var[(var.chrom == s.source_chrom) & (var.svtype == s.source_type)
              & (var.pos == s.source_pos) & (var.end == s.source_end)]
    assert len(hit), sv_id
    for i, r in hit.iterrows():
        g = G[i, midx]
        rows.append(dict(
            sv_id=sv_id, sv_locus=r.sv_locus, geno_group=int(r.geno_group),
            svtype=r.svtype, svlen=int(r.svlen), af=float(r.af),
            row_index=int(i), strategy=r.strategy, arm=r.arm, status=r.status,
            fidelity=r.fidelity, cv_r2=float(r.cv_r2),
            n_called=int((g >= 0).sum()),
            n_carrier=int((g > 0).sum()), n_het=int((g == 1).sum()),
            n_hom=int((g == 2).sum()), ac=int(g[g > 0].sum())))
det = pd.DataFrame(rows)
print("catalog rows matched at 0 bp: " + str(len(det))
      + "  over " + str(det.sv_id.nunique()) + " SVs, "
      + str(det.geno_group.nunique()) + " groups")

agree = []
for sv_id, blk in det.groupby("sv_id"):
    gg = sorted(blk.geno_group.unique())
    vecs = {int(g): G[blk[blk.geno_group == g].row_index.iloc[0], midx] for g in gg}
    if len(gg) == 1:
        agree.append(dict(sv_id=sv_id, n_groups=1, min_pair_agreement=1.0,
                          min_pair_r=1.0, identical=True))
        continue
    ags, rs = [], []
    for a in range(len(gg)):
        for b in range(a + 1, len(gg)):
            x, y = vecs[gg[a]].astype(float), vecs[gg[b]].astype(float)
            ok = (x >= 0) & (y >= 0)
            ags.append(float((x[ok] == y[ok]).mean()))
            rs.append(float(np.corrcoef(x[ok], y[ok])[0, 1])
                      if x[ok].std() > 0 and y[ok].std() > 0 else np.nan)
    agree.append(dict(sv_id=sv_id, n_groups=len(gg),
                      min_pair_agreement=float(np.min(ags)),
                      min_pair_r=float(np.nanmin(rs)),
                      identical=bool(np.min(ags) == 1.0)))
ag = pd.DataFrame(agree)
det = det.merge(ag, on="sv_id", how="left")

multi = ag[ag.n_groups > 1]
print("\nSVs mapping to >1 genotype-group: " + str(len(multi)) + " / " + str(len(ag)))
if len(multi):
    print("  genotype agreement between competing groups of the SAME locus:")
    print("    identical vectors : " + str(int(multi.identical.sum())) + " / " + str(len(multi)))
    print("    min pair agreement: median " + str(round(float(multi.min_pair_agreement.median()), 4))
          + ", min " + str(round(float(multi.min_pair_agreement.min()), 4)))
    print("    min pair Pearson r: median " + str(round(float(multi.min_pair_r.median()), 4))
          + ", min " + str(round(float(multi.min_pair_r.min()), 4)))

print("\nstatus values present: " + str(sorted(det.status.astype(str).unique())))
det["is_admit"] = (det.status.astype(str) == "ADMIT")
det["status_rank"] = (~det.is_admit).astype(int)
det = det.sort_values(["sv_id", "status_rank", "cv_r2", "geno_group"],
                      ascending=[True, True, False, True])
canon = det.groupby("sv_id", as_index=False).first()
canon["canonical_rule"] = "ADMIT_then_max_cv_r2_then_min_geno_group"

alt = (det.sort_values(["sv_id", "cv_r2", "geno_group"], ascending=[True, False, True])
       .groupby("sv_id", as_index=False).first()[["sv_id", "geno_group"]]
       .rename(columns={"geno_group": "geno_group_cv_r2_only"}))
canon = canon.merge(alt, on="sv_id", how="left")
canon["pick_changed_by_status"] = canon.geno_group != canon.geno_group_cv_r2_only

n_admit_avail = det.groupby("sv_id").is_admit.any()
canon["has_admit_row"] = canon.sv_id.map(n_admit_avail)
print("canonical units: " + str(len(canon)))
print("  canonical pick is ADMIT      : " + str(int(canon.is_admit.sum())) + " / " + str(len(canon)))
print("  loci with NO ADMIT row at all: " + str(int((~canon.has_admit_row).sum()))
      + "   <- genotype unreliable, must be flagged not silently used")
print("  picks changed by the status term vs cv_r2 alone: "
      + str(int(canon.pick_changed_by_status.sum())))
print("  cv_r2 of canonical pick, min/median: "
      + str(round(float(canon.cv_r2.min()), 4)) + " / "
      + str(round(float(canon.cv_r2.median()), 4)))

th = pd.read_csv(UNIT + "/output/a1_lift/pilot_qc/target_hit_summary.tsv", sep="\t")
rec = set(th[th.n_hit > 0].target_id)
pairs["recoverable_pilot8"] = pairs.target_id.isin(rec)
recp = pairs[pairs.recoverable_pilot8]
print("\nrecoverable pairs (pilot, 8 haplotypes): " + str(len(recp))
      + "   distinct source SVs: " + str(recp.sv_id.nunique()))
canon["contributes_recoverable"] = canon.sv_id.isin(set(recp.sv_id))
canon["n_pairs_all"] = canon.sv_id.map(pairs.groupby("sv_id").size()).fillna(0).astype(int)
canon["n_pairs_recoverable"] = canon.sv_id.map(
    recp.groupby("sv_id").size()).fillna(0).astype(int)

acct = pd.DataFrame([
    dict(quantity="report source SVs", n=int(pairs.sv_id.nunique()),
         population="the 465 pre-registered pairs"),
    dict(quantity="HPRC catalog rows at 0 bp", n=int(len(det)),
         population="all imputation strategies x arms at those coordinates"),
    dict(quantity="distinct genotype-groups", n=int(det.geno_group.nunique()),
         population="same, deduplicated on identical genotype vector (188 HPRC)"),
    dict(quantity="canonical genotype units", n=int(len(canon)),
         population="one per physical locus after the canonical rule"),
    dict(quantity="units contributing >=1 recoverable pair", n=int(canon.contributes_recoverable.sum()),
         population="pilot-8 recoverability; PROVISIONAL, recomputed after v2"),
    dict(quantity="units with n_carrier >= 5 in the methylation cohort",
         n=int((canon.n_carrier >= 5).sum()),
         population=str(len(in_geno)) + " individuals with both genotype and methylation"),
    dict(quantity="units with n_het >= 10 in the methylation cohort",
         n=int((canon.n_het >= 10).sum()), population="same"),
    dict(quantity="canonical picks with ADMIT status", n=int(canon.is_admit.sum()),
         population="the 47 canonical units"),
    dict(quantity="units with NO ADMIT row available", n=int((~canon.has_admit_row).sum()),
         population="genotype call unreliable; flagged, not dropped"),
    dict(quantity="multi-group loci (pick is a tie-break)", n=int((canon.n_groups > 1).sum()),
         population="require a canonical-choice sensitivity analysis"),
    dict(quantity="contributing units with n_het >= 10",
         n=int(((canon.n_het >= 10) & canon.contributes_recoverable).sum()),
         population="the haplotype DiD arm's actual cluster count"),
])
print("\n" + acct.to_string(index=False))

write_atomic(det, OUT + "/b0_sv_units/b0_genogroup_detail.tsv")
write_atomic(canon, OUT + "/b0_sv_units/b0_sv_unit_map.tsv")
write_atomic(acct, OUT + "/b0_sv_units/b0_reconciliation.tsv")

print("\ncanonical unit composition (contributing only):")
c = canon[canon.contributes_recoverable]
print("  type: " + str(dict(c.svtype.value_counts())))
print("  n_carrier  min/median/max: " + str(int(c.n_carrier.min())) + " / "
      + str(int(c.n_carrier.median())) + " / " + str(int(c.n_carrier.max())))
print("  n_het      min/median/max: " + str(int(c.n_het.min())) + " / "
      + str(int(c.n_het.median())) + " / " + str(int(c.n_het.max())))
print("  pairs/unit min/median/max: " + str(int(c.n_pairs_recoverable.min())) + " / "
      + str(int(c.n_pairs_recoverable.median())) + " / "
      + str(int(c.n_pairs_recoverable.max())))
