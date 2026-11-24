#!/usr/bin/env python3
import glob
import numpy as np
import pandas as pd
from pathlib import Path
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

DISC = CENSV_ROOT + "/analyses/cdr_qtl"
PERCHROM = f"{DISC}/results/cdr_qtl_corrected/perchrom"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"
N_PERM = 100000
FLOOR_P = 1.0 / (N_PERM + 1)
PHENOS = ["mean_5mC", "n_cdrs", "cdrs_per_mb", "cdr_fraction",
          "median_cdr_width", "cdr_pos_offset", "cdr_pos_var", "cdr_pos_iqr"]


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = p[i] * m / (rank + 1)
        prev = min(prev, val)
        q[i] = prev
    return np.clip(q, 0, 1)


def main():
    files = sorted(glob.glob(f"{PERCHROM}/*_perm.tsv"))
    print(f"[load] {len(files)} perchrom files")
    dfs = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        if d.empty:
            continue
        d["stratum"] = Path(f).stem.replace("_perm", "")
        dfs.append(d)
    alld = pd.concat(dfs, ignore_index=True)
    print(f"[load] total SV x pheno tests: {len(alld)}")
    print(f"[load] strata: {alld['stratum'].nunique()}  phenotypes: {sorted(alld.phenotype.unique())}")

    alld["q_within"] = np.nan
    for s, g in alld.groupby("stratum"):
        alld.loc[g.index, "q_within"] = bh(g["p_emp_perm"].values)

    alld["is_floor"] = alld["count_perm_ge_obs"] == 0
    n_floor = int(alld["is_floor"].sum())
    print(f"\n=== FLOOR HITS (count_perm_ge_obs==0, p_emp={FLOOR_P:.2e}) ===")
    print(f"  total floor hits: {n_floor}")
    by_ph = alld[alld.is_floor].groupby("phenotype").size().reindex(PHENOS).fillna(0).astype(int)
    print("  by phenotype:")
    for ph in PHENOS:
        nf = by_ph.get(ph, 0)
        nstrat = alld[alld.is_floor & (alld.phenotype == ph)]["stratum"].nunique()
        print(f"    {ph:18s}: {nf:4d} floor hits across {nstrat} chrom-strata")
    by_chrom = alld[alld.is_floor].copy()
    by_chrom["chrom_only"] = by_chrom["stratum"].str.split("_").str[0]
    print("  by chromosome (top):")
    print(by_chrom.groupby("chrom_only").size().sort_values(ascending=False).head(10).to_string())

    strat_min = alld.groupby("stratum")["p_emp_perm"].min().reset_index()
    strat_min["q_stage2"] = bh(strat_min["p_emp_perm"].values)
    n_sig_strata = int((strat_min["q_stage2"] < 0.05).sum())
    print(f"\n=== STAGE-2 BH on {len(strat_min)} stratum-min p ===")
    print(f"  strata with q_stage2<0.05: {n_sig_strata}")
    print(strat_min.sort_values("p_emp_perm").head(12).to_string(index=False))

    ph = pd.read_csv(PHENO, sep="\t")
    indiv = ph.groupby(["sample", "chm13_chrom"])[PHENOS].mean().reset_index()
    corr = indiv[PHENOS].corr().values
    corr = np.nan_to_num(corr, nan=0.0)
    eig = np.linalg.eigvalsh(corr)
    eig = np.clip(eig, 0, None)
    Meff = 1 + (len(PHENOS) - 1) * (1 - np.var(eig, ddof=1) / len(PHENOS))
    print(f"\n=== Li & Ji effective #phenotypes: Meff = {Meff:.2f} (of {len(PHENOS)}) ===")

    sig = alld[alld.q_within < 0.10].sort_values("p_emp_perm")
    print(f"\n=== top SVs by p_emp (q_within<0.10: {len(sig)}) ===")
    cols = ["stratum", "chrom", "pos", "svtype", "svlen", "af", "phenotype",
            "r_obs", "t_obs", "count_perm_ge_obs", "p_emp_perm", "q_within", "n_indiv"]
    cols = [c for c in cols if c in alld.columns]
    print(sig.head(25)[cols].to_string(index=False))

    print("\n=== chr8 DEL-342 cdr_pos_offset survivor (permutation) ===")
    s8 = alld[(alld.phenotype == "cdr_pos_offset") & (alld.chrom == "chr8") &
              (alld.svtype == "DEL") & (alld.svlen.between(340, 344)) &
              (alld.pos.between(44456000, 44457000))]
    if len(s8):
        print(s8[cols].to_string(index=False))
    else:
        print("  (DEL-342 @44456532 not in kept SVs; nearby chr8 cdr_pos_offset floor hits:)")
        print(alld[(alld.phenotype=="cdr_pos_offset")&(alld.chrom=="chr8")&alld.is_floor]
              .sort_values("p_emp_perm").head(10)[cols].to_string(index=False))

    out = f"{DISC}/output/scan/cohort_floor_hits_corrected.tsv"
    alld[alld.is_floor].sort_values("p_emp_perm").to_csv(out, sep="\t", index=False)
    alld.to_csv(f"{DISC}/output/scan/all_sv_pheno_tests_corrected.tsv.gz", sep="\t",
                index=False, compression="gzip")
    print(f"\n[written] {out} ({n_floor} floor hits)")
    print(f"[written] {DISC}/output/scan/all_sv_pheno_tests_corrected.tsv.gz ({len(alld)} tests)")
    print(f"\n>>> CORRECTED: {len(alld)} tests, {n_floor} floor hits @1e5, "
          f"{n_sig_strata} strata q<0.05.  PRIOR(over-called): 126 floor @1e5, 55 @1e6.")


if __name__ == "__main__":
    main()
