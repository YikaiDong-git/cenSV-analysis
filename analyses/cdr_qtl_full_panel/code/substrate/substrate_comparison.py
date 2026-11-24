#!/usr/bin/env python
"""Head-to-head: the 6,519-ADMIT substrate vs the full 11,129 substrate."""
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
NPH = 10


def memset(s):
    return frozenset(str(s).split(";"))


vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
adm = vt[vt.status == "ADMIT"]
gate = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
surv = gate[gate.overall == "SURVIVES"]
frames = []
for tag in ["1e5", "1e7"]:
    p = f"{WU}/output/all_tests_{tag}.tsv.gz"
    if os.path.exists(p):
        d = pd.read_csv(p, sep="\t"); d["prio"] = 1 if tag == "1e7" else 0; frames.append(d)
allt = pd.concat(frames, ignore_index=True)
best = allt.sort_values(["prio", "p_emp_perm"], ascending=[False, True]) \
           .drop_duplicates(["sv_locus", "phenotype"], keep="first")

rows = []
def add(k, a, b, note=""):
    rows.append({"quantity": k, "6,519 ADMIT": a, "11,129 FULL": b, "note": note})

add("defining filter", "+ UKBB LD-panel match", "no panel filter",
    "both share af>=0.0266 and cv_r2>=0.1")
add("rows (sv_locus)", f"{len(adm):,}", f"{len(vt):,}", f"+{len(vt)-len(adm):,}")
add("distinct genotype vectors", f"{adm.geno_group.nunique():,}", f"{vt.geno_group.nunique():,}",
    f"+{vt.geno_group.nunique()-adm.geno_group.nunique():,}")
add("distinct svim-asm member sets", f"{adm.sv_id.map(memset).nunique():,}",
    f"{vt.sv_id.map(memset).nunique():,}", "")
add("chromosomes", f"{adm.chrom.nunique()}", f"{vt.chrom.nunique()}", "+chr14, +chrX")
add("cv_r2 median", f"{adm.cv_r2.median():.3f}", f"{vt.cv_r2.median():.3f}",
    "excluded set 0.241 -- indistinguishable")
add("maf median", f"{adm.maf.median():.4f}", f"{vt.maf.median():.4f}", "full set is rarer")
add("maf < 0.05", f"{100*(adm.maf<0.05).mean():.1f}%", f"{100*(vt.maf<0.05).mean():.1f}%", "")
add("genotypes reconstructed exactly", "100%", "100%", "0 drift, 0 unrecovered in both")
add("mean missing cells / locus", f"{adm.n_miss.mean():.3f}", f"{vt.n_miss.mean():.3f}",
    "0.001% vs 0.005% of cells")
a_row, a_uni = 0.05/(len(adm)*NPH), 0.05/(adm.geno_group.nunique()*NPH)
f_row, f_uni = 0.05/(len(vt)*NPH), 0.05/(vt.geno_group.nunique()*NPH)
add("alpha, row-level (x10)", f"{a_row:.2e}", f"{f_row:.2e}", "full is 1.7x stricter")
add("alpha, unique-genotype (x10)", f"{a_uni:.2e}", f"{f_uni:.2e}", "")

sub_a = best[(best.status == "ADMIT") & (best.p_emp_perm <= a_uni)]
sub_f = best[best.p_emp_perm <= f_uni]
add("raw hits (own uniq-geno alpha)", f"{len(sub_a)} tests / {sub_a.geno_group.nunique()} groups",
    f"{len(sub_f)} tests / {sub_f.geno_group.nunique()} groups", "")
add("chroms with raw hits",
    ",".join(sorted(sub_a.chrom.unique(), key=lambda c: (len(c), c))),
    ",".join(sorted(sub_f.chrom.unique(), key=lambda c: (len(c), c))), "")

s_adm = surv[surv.status == "ADMIT"]
s_new = surv[surv.status != "ADMIT"]
add("leads surviving ALL gates", f"{len(s_adm)}", f"{len(surv)}",
    f"+{len(s_new)} (both chr9, cdr_fraction)")
add("chroms surviving", ",".join(sorted(set(s_adm.chrom))), ",".join(sorted(set(surv.chrom))), "")
add("phenotypes surviving", ",".join(sorted(set(s_adm.phenotype))),
    ",".join(sorted(set(surv.phenotype))), "")

cmp = pd.DataFrame(rows)
print(cmp.to_string(index=False))
out = f"{WU}/output/tables/substrate_comparison_6519_vs_11129.tsv"
cmp.to_csv(out, sep="\t", index=False)

print("\n=== what the extra 4,610 loci bought, per chromosome (loci added) ===")
t = pd.DataFrame({"admit": adm.groupby("chrom").size(), "full": vt.groupby("chrom").size()}).fillna(0).astype(int)
t["added"] = t.full - t.admit
t = t.sort_values("added", ascending=False)
print(t.head(10).to_string())

print("\n=== survivors, side by side ===")
cols = ["chrom", "pos", "svtype", "svlen", "af", "phenotype", "p_best", "status"]
print(surv[cols].sort_values(["chrom", "p_best"]).to_string(index=False))
print(f"\nwrote {out}")
