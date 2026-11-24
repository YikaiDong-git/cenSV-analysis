#!/usr/bin/env python3
"""RM-ALR spot-check (humility #6) WITHOUT RepeatMasker."""
import os, re, glob, pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT = CENSV_ROOT + ""
SUB = f"{ROOT}/subproject_methylation"
PHENO = f"{ROOT}/analyses/cdr_qtl/output/per_hap_chrom_phenotypes.tsv"
CENSAT = f"{SUB}/raw/censat"

def is_alpha(label):
    l = str(label).strip().lower()
    return ("hor" in l) or l.startswith("mon")

def overlap_len(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))

ph = pd.read_csv(PHENO, sep="\t")
c8 = ph[ph.chm13_chrom == "chr8"].copy()
print(f"chr8 haps in phenotype table: {len(c8)}")

def hap_of(contig):
    m = re.search(r"#(\d)#", str(contig)); return f"hap{m.group(1)}" if m else None

rows = []
for _, r in c8.iterrows():
    contig = r["asm_contig"]; hap = hap_of(contig)
    if hap is None: continue
    a0, a1 = int(r["active_start"]), int(r["active_end"])
    cands = glob.glob(f"{CENSAT}/{r['sample']}_{hap}_*.cenSat.bed")
    if not cands: continue
    cs = pd.read_csv(cands[0], sep="\t", header=None, usecols=[0,1,2,3],
                     names=["contig","s","e","lab"], dtype={"contig":str}, skiprows=1)
    cs = cs[cs.s.astype(str).str.fullmatch(r"\d+")]
    cs["s"] = cs["s"].astype(int); cs["e"] = cs["e"].astype(int)
    cs = cs[cs.contig == contig]
    if len(cs) == 0: continue
    alpha_bp = sum(overlap_len(a0,a1,int(s),int(e)) for s,e,l in zip(cs.s,cs.e,cs.lab) if is_alpha(l))
    nonalpha_bp = sum(overlap_len(a0,a1,int(s),int(e)) for s,e,l in zip(cs.s,cs.e,cs.lab) if not is_alpha(l))
    active_bp = a1 - a0
    annot_bp = alpha_bp + nonalpha_bp
    rows.append(dict(sample=r["sample"], hap=hap, active_mb=round(active_bp/1e6,2),
                     alpha_frac_of_annot=round(alpha_bp/annot_bp,4) if annot_bp else float("nan"),
                     alpha_frac_of_active=round(alpha_bp/active_bp,4) if active_bp else float("nan"),
                     nonalpha_kb=round(nonalpha_bp/1e3,1), n_cdrs=int(r["n_cdrs"])))

R = pd.DataFrame(rows)
print(f"\nresolved {len(R)} chr8 haps with cenSat annotation")
if len(R):
    print(f"\n=== alpha-satellite (ALR) fraction of the active region, chr8 ===")
    print(f"  alpha_frac_of_annotated: median={R.alpha_frac_of_annot.median():.4f}  min={R.alpha_frac_of_annot.min():.4f}  "
          f"q05={R.alpha_frac_of_annot.quantile(0.05):.4f}")
    print(f"  non-alpha satellite in active region (kb): median={R.nonalpha_kb.median():.1f}  max={R.nonalpha_kb.max():.1f}")
    nonalpha_heavy = R[R.alpha_frac_of_annot < 0.95]
    print(f"  haps with <95% alpha (ALR restriction MIGHT shift median): {len(nonalpha_heavy)}/{len(R)}")
    print(R.sort_values("alpha_frac_of_annot").head(6)[["sample","hap","active_mb","alpha_frac_of_annot","nonalpha_kb","n_cdrs"]].to_string(index=False))
    print(f"\nVERDICT: " + ("active region is ~all alpha-satellite => restrict_alr would NOT change the region median => CDR calls validated (active-HOR is a stricter subset, as claimed)" if R.alpha_frac_of_annot.median() >= 0.95 else "non-trivial non-alpha satellite in active region => restrict_alr could shift the median => needs a full re-call"))
    R.to_csv(f"{ROOT}/analyses/cdr_qtl/results/canonical_admit/rm_alr_spotcheck_chr8.tsv", sep="\t", index=False)
    print(f"wrote rm_alr_spotcheck_chr8.tsv")
