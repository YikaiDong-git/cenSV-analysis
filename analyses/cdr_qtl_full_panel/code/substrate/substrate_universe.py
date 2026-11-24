#!/usr/bin/env python
"""What IS the 11,129 set, exactly?"""
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
ADM = f"{ROOT}/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
SUF = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
EN = f"{ROOT}/results/elastic_net/decompose/cosine_nj"
SVD = f"{ROOT}/sv_extraction/decompose/{SUF}"

a = pd.read_csv(ADM, sep="\t")
print("=== (a) the 11,129 table ===")
print(f"rows: {len(a)}   unique sv_locus: {a.sv_locus.nunique()}")
print("statuses present:")
print(a.status.value_counts().to_string())
print(f"\ncv_r2: min={a.cv_r2.min():.4f}  max={a.cv_r2.max():.4f}  "
      f"n below 0.1 = {(a.cv_r2 < 0.1).sum()}")
print(f"FAIL_cv_r2_below_floor rows present: {(a.status == 'FAIL_cv_r2_below_floor').sum()}")
print(f"arms: {a.arm.value_counts().to_dict()}")
print(f"strategies: {a.strategy.nunique()} -> {sorted(a.strategy.unique())[:10]}")


def memset(s):
    return frozenset(str(s).split(";"))


print("\n=== (b) the DECOMPOSED SV universe actually present in the assemblies ===")
rows = []
for ch in sorted([f"chr{i}" for i in range(1, 23)] + ["chrX"], key=lambda c: (len(c), c)):
    p = f"{SVD}/{ch}/{ch}_sv_variants.tsv"
    n_decomp = None
    if os.path.exists(p):
        try:
            n_decomp = sum(1 for _ in open(p)) - 1
        except Exception:
            n_decomp = None
    q = None
    for s in ["midpoint1000", "midpoint500", "edge1000"]:
        cand = f"{EN}/{s}_{SUF}/{ch}/sv_metrics.tsv"
        if os.path.exists(cand):
            q = cand
            break
    n_metrics = (sum(1 for _ in open(q)) - 1) if q else None
    n_in_11129 = int((("chr" + a.primary_chrom.astype(str)) == ch).sum())
    rows.append(dict(chrom=ch, sv_variants_decomposed=n_decomp, sv_metrics_rows=n_metrics,
                     in_11129=n_in_11129))
t = pd.DataFrame(rows)
print(t.to_string(index=False))
print(f"\nTOTALS  decomposed sv_variants rows: {t.sv_variants_decomposed.sum():,}   "
      f"sv_metrics rows: {t.sv_metrics_rows.sum():,}   in the 11,129: {t.in_11129.sum():,}")

tot_uniq = 0
for ch in t.chrom:
    for s in ["midpoint1000", "midpoint500", "edge1000"]:
        p = f"{EN}/{s}_{SUF}/{ch}/sv_metrics.tsv"
        if os.path.exists(p):
            m = pd.read_csv(p, sep="\t", usecols=["sv_id"])
            tot_uniq += m.sv_id.map(memset).nunique()
            break
print(f"\nDISTINCT genotyped variants (unique svim-asm member sets in sv_metrics): {tot_uniq:,}")
print(f"Covered by the 11,129 weight-model table: {a.sv_id.map(memset).nunique():,}")
print(f"=> the CDR-QTL substrate covers {100*a.sv_id.map(memset).nunique()/max(tot_uniq,1):.1f}% "
      f"of the decomposed, genotyped SV catalogue.")
