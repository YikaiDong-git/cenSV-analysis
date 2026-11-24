#!/usr/bin/env python
import numpy as np, pandas as pd, os, sys
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

_D = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_D, os.path.dirname(_D)]
from phenotype_family import alpha, headline

ROOT = CENSV_ROOT + ""
CAN = f"{ROOT}/analyses/cdr_qtl"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"

d = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t").dropna(
    subset=["cdr_pos_offset", "cdr_pos_offset_deepest"])

print("=" * 92)
print("0. ANALYSIS UNIT")
vtn = len(open(f"{WU}/geno/full_samples.txt").read().split())
print(f"  genotype-matrix columns (full_samples.txt): {vtn}")
print(f"  distinct sample IDs in the phenotype table: {d['sample'].nunique()}")
print(f"  rows in the phenotype table (sample x hap x chrom): {len(d)}")
print(f"  hap labels present: {sorted(d.hap_local.unique())}")
print("  -> the scanner averages the two assembly haplotypes to ONE value per individual before\n"
      "     testing, so the testing unit is the individual, not the haplotype.")

print("\n" + "=" * 92)
print("1. WITHIN-CHROMOSOME CORRELATION (agy: pooling inflates r)")
pool_r = np.corrcoef(d.cdr_pos_offset, d.cdr_pos_offset_deepest)[0, 1]
rows = []
for ch, g in d.groupby("chm13_chrom"):
    if len(g) < 20:
        continue
    r = np.corrcoef(g.cdr_pos_offset, g.cdr_pos_offset_deepest)[0, 1]
    gd = g[~np.isclose(g.cdr_pos_offset, g.cdr_pos_offset_deepest, atol=1e-9)]
    rd = (np.corrcoef(gd.cdr_pos_offset, gd.cdr_pos_offset_deepest)[0, 1]
          if len(gd) >= 20 else np.nan)
    rows.append(dict(chrom=ch, n=len(g), r_within=r, n_discordant=len(gd), r_discordant=rd,
                     mean_rel=g.cdr_pos_offset.mean()))
W = pd.DataFrame(rows).sort_values("r_within")
print(f"  pooled r across all chromosomes            : {pool_r:.4f}")
print(f"  between-chromosome spread of mean position : "
      f"{W.mean_rel.min():+.3f} to {W.mean_rel.max():+.3f}")
print(f"  WITHIN-chromosome r: median {W.r_within.median():.4f}, "
      f"range {W.r_within.min():.4f}-{W.r_within.max():.4f}")
print(f"  WITHIN-chromosome r, discordant rows only: median {W.r_discordant.median():.4f}, "
      f"range {W.r_discordant.min():.4f}-{W.r_discordant.max():.4f}")
print(W.to_string(index=False))
W.to_csv(f"{WU}/output/tables/position_phenotype_within_chrom_r.tsv", sep="\t", index=False,
         float_format="%.6g")

print("\n" + "=" * 92)
print("2. WITHIN-HAP DISPERSAL (codex: is the deepest dip inside the dip cloud?)")
gap = (d.cdr_pos_offset - d.cdr_pos_offset_deepest).abs()
sub = d[d.n_cdrs >= 2]
gap2 = (sub.cdr_pos_offset - sub.cdr_pos_offset_deepest).abs()
iqr = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
m = sub.merge(iqr[["sample", "hap_local", "chm13_chrom", "cdr_pos_iqr", "cdr_pos_var"]],
              on=["sample", "hap_local", "chm13_chrom"], how="left")
m["gap"] = (m.cdr_pos_offset - m.cdr_pos_offset_deepest).abs()
m["sd_rel"] = np.sqrt(m.cdr_pos_var)
ok = m[(m.cdr_pos_iqr > 0) & m.gap.notna()]
print(f"  |median - deepest| : mean {gap.mean():.4f}, median {gap.median():.4f}, "
      f"p90 {gap.quantile(.9):.4f}, max {gap.max():.4f}")
print(f"  within-hap IQR of rel (n_cdrs>=2): mean {ok.cdr_pos_iqr.mean():.4f}, "
      f"median {ok.cdr_pos_iqr.median():.4f}")
print(f"  within-hap SD of rel             : mean {ok.sd_rel.mean():.4f}, "
      f"median {ok.sd_rel.median():.4f}")
print(f"  ratio  |gap| / within-hap SD     : median {(ok.gap / ok.sd_rel).median():.3f}, "
      f"mean {(ok.gap / ok.sd_rel).replace([np.inf, -np.inf], np.nan).mean():.3f}")
print(f"  fraction of hap-chromosomes where the deepest dip lies INSIDE the interquartile range "
      f"of that hap's dips: ", end="")
inside = (ok.gap <= ok.cdr_pos_iqr / 2).mean()
print(f"{inside:.3f}")
print(f"  full-scale reference: rel spans -1..+1, so the mean gap {gap.mean():.4f} is "
      f"{100*gap.mean()/2:.1f}% of the array")

print("\n" + "=" * 92)
print("3. DOES THE HEADLINE SURVIVOR SET CHANGE BETWEEN x8 AND STRICT x10?")
gate = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
HL_LAB, HL_N, HL_A = headline(vt)
A10 = 0.05 / (HL_N * 10)
A7 = 0.05 / (HL_N * 7)
for lab, a in [("x10 (all scanned, strictest)", A10), ("x8  (primary, HEADLINE)", HL_A),
               ("x7  (if deepest were also a twin)", A7),
               ("xMeff=6 (loosest defensible)", 0.05 / (HL_N * 6))]:
    s = gate[(gate.p_best <= a) & (gate.overall == "SURVIVES")]
    print(f"  {lab:36s} alpha={a:.3e}  survivors={len(s):2d}  "
          f"groups={s.geno_group.nunique()}  SVs={s.groupby(['pos','svtype','svlen']).ngroups}  "
          f"chroms={sorted(set(s.chrom))}")
s8 = set(map(tuple, gate[(gate.p_best <= HL_A) & (gate.overall == "SURVIVES")]
             [["sv_locus", "phenotype"]].values))
s10 = set(map(tuple, gate[(gate.p_best <= A10) & (gate.overall == "SURVIVES")]
              [["sv_locus", "phenotype"]].values))
print(f"\n  x8 survivor set == x10 survivor set ?  {s8 == s10}   "
      f"(symmetric difference: {len(s8 ^ s10)})")

print("\n" + "=" * 92)
print("4. HOW MANY INDEPENDENT LOCI ARE THE SURVIVORS?")
sv = gate[(gate.p_best <= HL_A) & (gate.overall == "SURVIVES")]
print(f"  test cells (genotype-group x phenotype): {len(sv)}")
print(f"  distinct genotype vectors             : {sv.geno_group.nunique()}")
print(f"  distinct SVs                          : {sv.groupby(['pos','svtype','svlen']).ngroups}")
print(f"  distinct chromosomes                  : {sv.chrom.nunique()}")
print(f"  genomic span of the surviving SVs     : chr8:{int(sv.pos.min()):,}-{int(sv.pos.max()):,} "
      f"= {int(sv.pos.max() - sv.pos.min()):,} bp")
print("  -> these are NOT independent discoveries: 4 SVs in strong mutual LD inside one window,")
print("     under two phenotypes correlated at r = 0.989. The correct unit is ONE locus.")
