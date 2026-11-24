#!/usr/bin/env python3
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

warnings.filterwarnings("ignore")

ROOT = Path(CENSV_ROOT + "")
SUFFIX = "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
OUT = ROOT / "analyses" / "imputability" / "output"
FA = ROOT / "resources" / "chm13v2.0.fa.gz"
CENSAT = ROOT / "resources" / "chm13v2.0_censat_v2.1.bed"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def family_class(name):
    return name.split("_")[0]


def gc_in_window(fa, chrom, pos, half):
    start = max(0, pos - half); end = pos + half
    seq = fa.fetch(chrom, start, end).upper()
    if len(seq) == 0: return np.nan
    a = seq.count("A"); c = seq.count("C")
    g = seq.count("G"); t = seq.count("T")
    denom = a + c + g + t
    return (g + c) / denom if denom > 0 else np.nan


def main():
    full_rows = []
    for chrom in CHROMS:
        var = SV_BASE / chrom / f"{chrom}_sv_variants.tsv"
        if not var.exists():
            continue
        v = pd.read_csv(var, sep="\t", low_memory=False)
        v["chrom"] = chrom
        v = v.reset_index().rename(columns={"index": "sv_index"})
        full_rows.append(v)
    full = pd.concat(full_rows, ignore_index=True)
    print(f"[INFO] full SV catalog: {len(full)} SVs across "
          f"{full['chrom'].nunique()} chroms", flush=True)

    imp = pd.read_csv(OUT / "per_sv_features_v3.tsv", sep="\t", low_memory=False,
                      usecols=["chrom", "sv_index", "cv_r2"])
    imp["sv_index"] = imp["sv_index"].astype(int)
    imp_keys = set(zip(imp["chrom"], imp["sv_index"]))
    full["sv_index"] = full["sv_index"].astype(int)
    full["imputable"] = full.apply(
        lambda r: 1 if (r["chrom"], r["sv_index"]) in imp_keys else 0, axis=1)
    print(f"[INFO] imputable: {full['imputable'].sum()} / {len(full)}  "
          f"({full['imputable'].mean()*100:.1f}%)", flush=True)

    if "svlen" in full.columns:
        full["log_svlen"] = np.log10(np.maximum(full["svlen"].astype(float), 1))
    full["svtype_INS"] = (full["svtype"] == "INS").astype(int)
    full["svtype_DEL"] = (full["svtype"] == "DEL").astype(int)
    full["svtype_INV"] = (full["svtype"] == "INV").astype(int)
    full["svtype_DUP"] = (full["svtype"].str.startswith("DUP")).astype(int)

    anc = pd.read_csv(OUT / "sv_ancestry_features.tsv", sep="\t", low_memory=False)
    anc["sv_index"] = anc["sv_index"].astype(int)
    full = full.merge(anc, on=["chrom", "sv_index"], how="left",
                       suffixes=("", "_anc"))
    print(f"[INFO] ancestry merge coverage: "
          f"{full['af_skew_sv'].notna().mean()*100:.1f}%", flush=True)

    chrom_feats = imp.copy()
    chrom_lookup = {}
    base_imp = pd.read_csv(OUT / "per_sv_features_v3.tsv", sep="\t",
                           low_memory=False)
    for chrom in CHROMS:
        chr_sub = base_imp[base_imp["chrom"] == chrom]
        if len(chr_sub) == 0:
            continue
        chrom_lookup[chrom] = {
            "log_active_HOR_span": chr_sub["log_active_HOR_span"].iloc[0],
            "n_HOR_intervals": chr_sub["n_HOR_intervals"].iloc[0],
            "n_train_haps": chr_sub["n_train_haps"].iloc[0],
            "acrocentric": chr_sub["acrocentric"].iloc[0],
            "pericen_len_500000": chr_sub["pericen_len_500000"].iloc[0],
        }
    for col in ["log_active_HOR_span", "n_HOR_intervals", "n_train_haps",
                 "acrocentric", "pericen_len_500000"]:
        full[col] = full["chrom"].map(lambda c: chrom_lookup.get(c, {}).get(col, np.nan))

    print("[INFO] computing GC + sat features for full catalog...", flush=True)
    fa = pysam.FastaFile(str(FA))
    cs = pd.read_csv(CENSAT, sep="\t", header=None,
                     names=["chrom", "start", "end", "name",
                            "score", "strand", "ts", "te", "rgb"])
    cs["family"] = cs["name"].apply(family_class)

    pos_col = ("pos" if "pos" in full.columns else
               ("svpos" if "svpos" in full.columns else "start"))
    print(f"[INFO] using {pos_col} as SV pos column", flush=True)
    if pos_col not in full.columns:
        print("[FAIL] no pos column", flush=True); sys.exit(1)

    gc_10kb = []; gc_100kb = []; gc_1Mb = []
    sat_local = []
    for i, row in full.iterrows():
        chrom = row["chrom"]; pos = int(row[pos_col]) if pd.notna(row[pos_col]) else 0
        gc_10kb.append(gc_in_window(fa, chrom, pos, 10_000))
        gc_100kb.append(gc_in_window(fa, chrom, pos, 100_000))
        gc_1Mb.append(gc_in_window(fa, chrom, pos, 1_000_000))
        cs_chr = cs[cs["chrom"] == chrom]
        local = cs_chr[(cs_chr["end"] > pos - 100_000) &
                        (cs_chr["start"] < pos + 100_000)]
        sat_local.append(local["family"].nunique())
        if (i + 1) % 10000 == 0:
            print(f"  ...{i+1}/{len(full)}", flush=True)
    full["gc_10kb"] = gc_10kb
    full["gc_100kb"] = gc_100kb
    full["gc_1Mb"] = gc_1Mb
    full["sat_local_complexity"] = sat_local

    feats = [
        "log_svlen", "svtype_INS", "svtype_DEL", "svtype_INV", "svtype_DUP",
        "log_active_HOR_span", "n_HOR_intervals", "n_train_haps",
        "acrocentric", "pericen_len_500000",
        "gc_10kb", "gc_100kb", "gc_1Mb", "sat_local_complexity",
        "af_AFR_sv", "af_AMR_sv", "af_EAS_sv", "af_EUR_sv", "af_SAS_sv",
        "af_skew_sv", "max_FST_pair_sv", "total_af_sv",
    ]
    feats = [f for f in feats if f in full.columns]
    X = full[feats].astype(float).values
    y = full["imputable"].astype(int).values

    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    fold_aucs = []
    fold_perm = []
    for fold_i, (tr, te) in enumerate(kf.split(X)):
        m = HistGradientBoostingClassifier(
            max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
            min_samples_leaf=50, l2_regularization=1.0, random_state=42)
        m.fit(X[tr], y[tr])
        auc = m.score(X[te], y[te])
        from sklearn.metrics import roc_auc_score
        proba = m.predict_proba(X[te])[:, 1]
        roc = roc_auc_score(y[te], proba)
        fold_aucs.append(roc)
        pi = permutation_importance(m, X[te], y[te], n_repeats=10,
                                     random_state=42, n_jobs=2,
                                     scoring="roc_auc")
        fold_perm.append(pi.importances_mean)
        print(f"  fold{fold_i}: ROC-AUC = {roc:.4f}", flush=True)

    auc_arr = np.array(fold_aucs)
    perm_mean = np.mean(np.vstack(fold_perm), axis=0)
    perm_std = np.std(np.vstack(fold_perm), axis=0)

    summary = pd.DataFrame([{
        "n_total": len(full),
        "n_imputable": int(full["imputable"].sum()),
        "frac_imputable": float(full["imputable"].mean()),
        "roc_auc_mean_5fold": float(auc_arr.mean()),
        "roc_auc_std_5fold": float(auc_arr.std()),
        "roc_auc_min_5fold": float(auc_arr.min()),
        "roc_auc_max_5fold": float(auc_arr.max()),
    }])
    summary.to_csv(OUT / "v3_survivor_bias_summary.tsv", sep="\t", index=False)
    print(f"\n=== Survivor-bias classifier (imputable vs unimputable) ===",
          flush=True)
    print(summary.to_string(index=False), flush=True)

    perm = pd.DataFrame({
        "feature": feats,
        "perm_imp_mean": perm_mean,
        "perm_imp_std": perm_std,
    }).sort_values("perm_imp_mean", ascending=False)
    perm.to_csv(OUT / "v3_survivor_bias_perm.tsv", sep="\t", index=False)
    print(f"\nTop-10 features for imputable-vs-unimputable:", flush=True)
    print(perm.head(10).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
