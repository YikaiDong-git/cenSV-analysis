#!/usr/bin/env python3
import argparse
import struct
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT_DEFAULT = CENSV_ROOT + ""
PHENO_DEFAULT = (ROOT_DEFAULT +
    "/analyses/cdr_qtl/output/per_hap_chrom_phenotypes.tsv")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chrom", required=True)
    p.add_argument("--phenotype", required=True,
                   choices=["mean_5mC", "mean_5mC_covwt", "n_cdrs", "cdrs_per_mb",
                            "cdr_fraction", "median_cdr_width", "cdr_pos_offset",
                            "cdr_pos_var", "cdr_pos_iqr", "cdr_pos_offset_deepest"])
    p.add_argument("--n-perm", type=int, default=100_000)
    p.add_argument("--out-dir", default="results/cdr_qtl_corrected")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--root", default=ROOT_DEFAULT)
    p.add_argument("--pheno-file", default=PHENO_DEFAULT)
    p.add_argument("--cohort", choices=["all", "r941"], default="all")
    p.add_argument("--grm-prefix", default=None,
                   help="GRM prefix; default pericen_1m_all (all) / pericen_1m_r941 (r941)")
    p.add_argument("--extra-cols", default="",
                   help="Comma-separated extra columns from <chrom>_sv_variants.tsv to carry into the output "
                        "(e.g. sv_locus,geno_group,status). Default empty = unchanged output schema.")
    p.add_argument("--sv-source", default=None,
                   help="Dir holding <chrom>/<chrom>_sv_{genotypes.npz,samples.txt,variants.tsv}. "
                        "Default = legacy decompose. Set to the canonical 6519 perchrom dir for the canonical scan.")
    return p.parse_args()


def read_gcta_grm_bin(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), dtype=np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v
                K[j, i] = v
    ids = [line.strip().split() for line in open(f"{prefix}.grm.id")]
    return K, [idx[1] for idx in ids]


def main():
    args = parse_args()
    ROOT = Path(args.root)
    SUB = ROOT / "subproject_methylation"
    chrom = args.chrom
    pheno_name = args.phenotype
    N_PERM = args.n_perm
    SEED = args.seed
    OUT_DIR = Path(args.out_dir)
    if not OUT_DIR.is_absolute():
        OUT_DIR = ROOT / "analyses/cdr_qtl" / OUT_DIR
    perchrom_dir = OUT_DIR / "perchrom"
    perchrom_dir.mkdir(parents=True, exist_ok=True)
    grm_prefix = args.grm_prefix or str(
        SUB / ("results/greml/pericen_1m_all" if args.cohort == "all"
               else "results/greml/pericen_1m_r941"))

    print(f"[CONFIG] chrom={chrom} pheno={pheno_name} n_perm={N_PERM:,} seed={SEED} "
          f"cohort={args.cohort} grm={grm_prefix}")
    print(f"[CONFIG] pheno_file={args.pheno_file}")
    print(f"[CONFIG] out_dir={OUT_DIR}")

    if args.sv_source:
        SV_DIR = Path(args.sv_source) / chrom
    else:
        SV_DIR = (ROOT / "sv_extraction/decompose"
                  / "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
                  / chrom)
    npz_path = SV_DIR / f"{chrom}_sv_genotypes.npz"
    print(f"[1/6] Loading SVs from {npz_path}")
    npz_data = np.load(npz_path)
    key = "arr_0" if "arr_0" in npz_data else list(npz_data.keys())[0]
    sv_dosages = npz_data[key]
    sv_samples = open(SV_DIR / f"{chrom}_sv_samples.txt").read().split()
    sv_variants = pd.read_csv(SV_DIR / f"{chrom}_sv_variants.tsv", sep="\t")
    if sv_dosages.shape[0] == len(sv_samples):
        pass
    elif sv_dosages.shape[1] == len(sv_samples):
        sv_dosages = sv_dosages.T
    else:
        sys.exit(f"[FATAL] SV matrix {sv_dosages.shape} matches neither n_samples nor T")
    print(f"  oriented SV matrix: {sv_dosages.shape}  n_variants={len(sv_variants)}")

    print("[2/6] Loading CORRECTED per-hap-chrom phenotypes")
    pc = pd.read_csv(args.pheno_file, sep="\t")
    pc_chrom = pc[pc["chm13_chrom"] == chrom]
    print(f"  {chrom} rows: {len(pc_chrom)}")
    pheno_ind = pc_chrom.groupby("sample")[pheno_name].mean().reset_index()
    bc_map = pc_chrom.groupby("sample")["basecaller"].first().to_dict()
    print(f"  per-individual rows: {len(pheno_ind)}  non-NaN {pheno_name}: "
          f"{pheno_ind[pheno_name].notna().sum()}")

    if args.cohort == "r941":
        r941_set = {s for s, b in bc_map.items() if b == "R941"}
        cohort_set = r941_set
        print(f"[3/6] cohort=r941: {len(cohort_set)} samples")
    else:
        cohort_set = set(bc_map.keys())
        print(f"[3/6] cohort=all: {len(cohort_set)} samples "
              f"(R1041 frac={np.mean([bc_map[s]=='R1041' for s in cohort_set]):.2f})")

    print(f"[4/6] PCs from GRM {grm_prefix}")
    K, grm_samples = read_gcta_grm_bin(grm_prefix)
    eigvals, eigvecs = np.linalg.eigh(K)
    order = np.argsort(eigvals)[::-1]
    PCs = eigvecs[:, order[:10]]
    pc_df = pd.DataFrame(PCs, columns=[f"PC{i+1}" for i in range(10)])
    pc_df["sample"] = grm_samples

    print("[5/6] Aligning sample sets")
    keep = sorted(cohort_set & set(pheno_ind["sample"]) & set(grm_samples) & set(sv_samples))
    print(f"  intersection: n={len(keep)}")
    if len(keep) < 30:
        sys.exit(f"[FATAL] intersection too small: n={len(keep)}")

    pheno_map = pheno_ind.set_index("sample")[pheno_name].to_dict()
    pc_map = pc_df.set_index("sample")[[f"PC{i+1}" for i in range(10)]].to_dict("index")
    sv_idx_map = {s: i for i, s in enumerate(sv_samples)}

    y = np.array([pheno_map[s] for s in keep], dtype=float)
    PC_mat = np.array([[pc_map[s][f"PC{i+1}"] for i in range(10)] for s in keep], dtype=float)
    bc_vec = np.array([1.0 if bc_map.get(s) == "R1041" else 0.0 for s in keep], dtype=float)
    sv_subset = sv_dosages[[sv_idx_map[s] for s in keep], :].astype(float)

    y_mask = ~np.isnan(y)
    y = y[y_mask]; PC_mat = PC_mat[y_mask]; bc_vec = bc_vec[y_mask]
    sv_subset = sv_subset[y_mask]
    keep = [k for i, k in enumerate(keep) if y_mask[i]]
    print(f"  after NaN drop: n={len(y)}")
    if len(y) < 30:
        sys.exit(f"[FATAL] post-NaN n={len(y)} too small")

    def rint(a):
        r = stats.rankdata(a)
        return stats.norm.ppf((r - 0.5) / len(r))

    y_rint = rint(y)

    cov_cols = [np.ones(len(y_rint)), PC_mat]
    n_cov = 1 + 10
    if args.cohort == "all" and bc_vec.std() > 1e-10:
        cov_cols = [np.ones(len(y_rint)), PC_mat, bc_vec[:, None]]
        n_cov = 1 + 10 + 1
    X_pc = np.column_stack(cov_cols)
    df_resid = len(y_rint) - X_pc.shape[1] - 1
    print(f"  covariates={X_pc.shape[1]} (df_resid=n-{X_pc.shape[1]}-1={df_resid})")

    beta_pc = np.linalg.lstsq(X_pc, y_rint, rcond=None)[0]
    y_resid = y_rint - X_pc @ beta_pc

    print(f"[6/6] {N_PERM:,} within-chrom perms on {chrom} (n={len(y)})")
    rng = np.random.default_rng(SEED)
    sv_resid = sv_subset.copy()
    for j in range(sv_resid.shape[1]):
        bj = np.linalg.lstsq(X_pc, sv_resid[:, j], rcond=None)[0]
        sv_resid[:, j] = sv_resid[:, j] - X_pc @ bj
    sd = sv_resid.std(axis=0)
    mask = sd > 1e-10
    sv_resid = sv_resid[:, mask]
    sv_variants_kept = sv_variants.iloc[mask.nonzero()[0]].reset_index(drop=True)
    sv_resid = (sv_resid - sv_resid.mean(axis=0)) / sv_resid.std(axis=0)
    n = len(y_resid)
    y_resid_std = (y_resid - y_resid.mean()) / y_resid.std()
    print(f"  kept n_svs: {sv_resid.shape[1]}")

    r_obs = (sv_resid.T @ y_resid_std) / n
    t_obs = r_obs * np.sqrt(df_resid) / np.sqrt(np.maximum(1 - r_obs ** 2, 1e-30))
    abs_t_obs = np.abs(t_obs)

    counts = np.zeros(sv_resid.shape[1], dtype=np.int64)
    batch = 1000
    n_batches = N_PERM // batch
    for b in range(n_batches):
        Y_perm = np.empty((batch, n))
        for k in range(batch):
            Y_perm[k] = y_resid_std[rng.permutation(n)]
        R_perm = (Y_perm @ sv_resid) / n
        T_perm = R_perm * np.sqrt(df_resid) / np.sqrt(np.maximum(1 - R_perm ** 2, 1e-30))
        counts += (np.abs(T_perm) >= abs_t_obs[None, :]).sum(axis=0)
    p_emp = (counts + 1) / (N_PERM + 1)

    KEEP = ["chrom", "pos", "end", "svtype", "svlen", "ac", "an", "af", "maf", "n_miss"]
    KEEP = KEEP + [c.strip() for c in args.extra_cols.split(",") if c.strip()]
    keep_cols = [c for c in KEEP if c in sv_variants_kept.columns]
    out = sv_variants_kept[keep_cols].copy()
    if "chrom" not in out.columns:
        out["chrom"] = chrom
    out["t_obs"] = t_obs
    out["abs_t_obs"] = abs_t_obs
    out["r_obs"] = r_obs
    out["count_perm_ge_obs"] = counts
    out["p_emp_perm"] = p_emp
    out["phenotype"] = pheno_name
    out["n_indiv"] = n
    out["cohort"] = args.cohort
    out = out.sort_values("p_emp_perm").reset_index(drop=True)

    out_path = perchrom_dir / f"{chrom}_{pheno_name}_perm.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"WROTE {out_path}: {len(out)} SVs")
    print(out.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
