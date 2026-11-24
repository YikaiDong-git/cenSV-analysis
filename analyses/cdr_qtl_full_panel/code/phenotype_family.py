#!/usr/bin/env python
"""SINGLE SOURCE OF TRUTH for the CDR-QTL phenotype family and the multiple-testing threshold."""
import os
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = (CENSV_ROOT + "/"
      "analyses/cdr_qtl_full_panel")

PRIMARY = [
    "cdr_pos_offset",
    "cdr_pos_offset_deepest",
    "cdr_pos_iqr",
    "n_cdrs",
    "cdrs_per_mb",
    "cdr_fraction",
    "median_cdr_width",
    "mean_5mC",
]

SUPPLEMENTARY = {
    "cdr_pos_var": ("cdr_pos_iqr",
                    "spread of dip positions as a moment (variance) rather than an order "
                    "statistic (IQR)"),
    "mean_5mC_covwt": ("mean_5mC",
                       "array methylation level weighted by read coverage rather than unweighted"),
}

ALL_SCANNED = PRIMARY + list(SUPPLEMENTARY)
NPH = len(PRIMARY)
assert NPH == 8, NPH
assert len(ALL_SCANNED) == 10, ALL_SCANNED

GROUPS = [
    ("Position / geometry -- WHERE the dip sits",
     ["cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_pos_iqr", "cdr_pos_var"]),
    ("Abundance / extent -- HOW MUCH of the array is hypomethylated",
     ["n_cdrs", "cdrs_per_mb", "cdr_fraction", "median_cdr_width"]),
    ("Level -- HOW METHYLATED the array is overall",
     ["mean_5mC", "mean_5mC_covwt"]),
]


def denominators(variant_table=None):
    vt = variant_table if variant_table is not None else pd.read_csv(
        f"{WU}/geno/full_variants.tsv", sep="\t")
    adm = vt[vt.status == "ADMIT"]
    return [
        ("row-level, FULL (HEADLINE)", len(vt)),
        ("unique-genotype, FULL", int(vt.geno_group.nunique())),
        ("row-level, ADMIT (predecessor)", len(adm)),
        ("unique-genotype, ADMIT", int(adm.geno_group.nunique())),
    ]


def alpha(n, nph=NPH):
    return 0.05 / (n * nph)


def headline(variant_table=None):
    lab, n = denominators(variant_table)[0]
    return lab, n, alpha(n)


def is_primary(ph):
    return ph in PRIMARY


def label(ph):
    if ph in SUPPLEMENTARY:
        return f"{ph}  [supp: twin of {SUPPLEMENTARY[ph][0]}]"
    return ph


if __name__ == "__main__":
    vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
    print(f"PRIMARY ({NPH}):        {PRIMARY}")
    print(f"SUPPLEMENTARY ({len(SUPPLEMENTARY)}):  "
          f"{ {k: v[0] for k, v in SUPPLEMENTARY.items()} }")
    print(f"\nalpha = 0.05 / (n x {NPH})")
    for nm, n in denominators(vt):
        print(f"  {nm:34s} n={n:6d}   alpha={alpha(n):.3e}")
    lab, n, a = headline(vt)
    print(f"\nHEADLINE: {lab}  n={n}  alpha={a:.3e}")
