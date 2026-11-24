#!/usr/bin/env python
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BASE = CENSV_ROOT + "/discovery"
V1D = BASE + "/cdr_qtl/output"
RD = BASE + "/cdr_phenotypes"
OUT = RD + "/output"

PAIRS = [
    ("cdr_pos_offset", "cdr_pos_offset_span"),
    ("cdr_pos_offset_deepest", "cdr_pos_offset_deepest_span"),
    ("cdr_pos_iqr", "cdr_pos_iqr_span"),
    ("cdr_pos_var", "cdr_pos_var_span"),
    ("cdr_fraction", "cdr_fraction_sum"),
    ("cdrs_per_mb", "cdrs_per_mb_sum"),
    ("n_cdrs", "n_cdrs"),
    ("median_cdr_width", "median_cdr_width"),
    ("mean_5mC", "mean_5mC"),
    ("mean_5mC_covwt", "mean_5mC_covwt"),
]
FRAMES = ["cdr_pos_offset_span", "cdr_pos_offset_block", "cdr_pos_offset_splice"]
NEW_PHENOS = ["cdr_pos_offset_block", "cdr_pos_offset_splice",
              "cdr_pos_offset_deepest_block", "cdr_pos_offset_deepest_splice"]


def rint(x):
    x = np.asarray(x, dtype=float)
    ok = ~np.isnan(x)
    out = np.full(x.shape, np.nan)
    if ok.sum() < 3:
        return out
    out[ok] = norm.ppf((rankdata(x[ok]) - 0.5) / ok.sum())
    return out


def per_individual(df, col, newname):
    if col not in df.columns:
        return None
    g = (df[["sample", "chm13_chrom", col]]
         .groupby(["sample", "chm13_chrom"], as_index=False)[col].mean())
    return g.rename(columns={col: newname})


def main():
    v1 = pd.read_csv(V1D + "/per_hap_chrom_phenotypes.tsv", sep="\t", engine="python")
    deep = pd.read_csv(V1D + "/per_hap_chrom_phenotypes_deepest.tsv", sep="\t", engine="python")
    key = ["sample", "hap_local", "chm13_chrom"]
    if "cdr_pos_offset_deepest" in deep.columns:
        v1 = v1.merge(deep[key + ["cdr_pos_offset_deepest"]], on=key, how="left")
    v2 = pd.read_csv(OUT + "/per_hap_chrom_phenotypes_frames.tsv", sep="\t", engine="python")

    print("v1 rows %d (%d samples) | v2 rows %d (%d samples)"
          % (len(v1), v1["sample"].nunique(), len(v2), v2["sample"].nunique()))

    rows, tasks = [], []
    for c1, c2 in PAIRS:
        a = per_individual(v1, c1, "V1")
        b = per_individual(v2, c2, "V2")
        if a is None or b is None:
            print("  SKIP %s -> %s (column missing)" % (c1, c2))
            continue
        m = a.merge(b, on=["sample", "chm13_chrom"], how="inner")
        for chrom, sub in m.groupby("chm13_chrom"):
            x, y = rint(sub["V1"].values), rint(sub["V2"].values)
            ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.sum() < 3:
                continue
            d = np.abs(x[ok] - y[ok])
            r = np.corrcoef(x[ok], y[ok])[0, 1] if np.std(x[ok]) > 0 and np.std(y[ok]) > 0 else np.nan
            changed = bool(d.max() > 1e-9)
            rows.append(dict(chrom=chrom, phenotype_v1=c1, phenotype_v2=c2, n=int(ok.sum()),
                             r=r, max_abs_delta=float(d.max()), mean_abs_delta=float(d.mean()),
                             verdict="CHANGED" if changed else "identical"))
            if changed:
                tasks.append((chrom, c2))

    R = pd.DataFrame(rows)
    R.to_csv(OUT + "/phenotype_diff_v1_v2.tsv", sep="\t", index=False, float_format="%.6g")

    agg = {}
    for f in FRAMES:
        agg[f] = per_individual(v2, f, f)
    fr = []
    for chrom, sub in v2.groupby("chm13_chrom"):
        rec = dict(chrom=chrom, n_hap=int(len(sub)),
                   pct_multiblock=100.0 * float((sub.n_blocks > 1).mean()))
        for i in range(len(FRAMES)):
            for j in range(i + 1, len(FRAMES)):
                p = agg[FRAMES[i]].merge(agg[FRAMES[j]], on=["sample", "chm13_chrom"])
                p = p[p.chm13_chrom == chrom]
                x, y = rint(p[FRAMES[i]].values), rint(p[FRAMES[j]].values)
                ok = ~np.isnan(x) & ~np.isnan(y)
                nm = "r_%s_vs_%s" % (FRAMES[i].split("_")[-1], FRAMES[j].split("_")[-1])
                rec[nm] = (np.corrcoef(x[ok], y[ok])[0, 1]
                           if ok.sum() > 2 and np.std(x[ok]) > 0 and np.std(y[ok]) > 0 else np.nan)
        fr.append(rec)
    F = pd.DataFrame(fr).sort_values("pct_multiblock", ascending=False)
    F.to_csv(OUT + "/frame_comparison.tsv", sep="\t", index=False, float_format="%.6g")

    for chrom in sorted(v2.chm13_chrom.unique()):
        for f in NEW_PHENOS:
            tasks.append((chrom, f))
    tasks = sorted(set(tasks), key=lambda t: (t[1], t[0]))
    with open(OUT + "/rescan_tasks.txt", "w") as fh:
        for chrom, ph in tasks:
            fh.write("%s\t%s\n" % (chrom, ph))

    ch = R[R.verdict == "CHANGED"]
    print("\n%d of %d (chrom, phenotype) cells CHANGED, %d identical"
          % (len(ch), len(R), len(R) - len(ch)))
    if len(ch):
        print("\nchanged cells per chromosome:")
        print(ch.groupby("chrom").size().to_string())
        print("\nlargest movers (rank-inverse-normal units):")
        print(ch.nlargest(15, "max_abs_delta")[
            ["chrom", "phenotype_v1", "n", "r", "max_abs_delta"]].to_string(index=False))
    print("\nrescan task list: %d tasks -> %s/rescan_tasks.txt" % (len(tasks), OUT))
    print("\nframe agreement per chromosome (top 12 by %% multi-block):")
    print(F.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
