#!/usr/bin/env python3
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

sys.path.insert(0, CENSV_ROOT + "/"
                   "cenSV_web/etl")
import tracks                                                              # noqa: E402

ROOT = Path(CENSV_ROOT + "")
SUFFIX = ("mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_"
          "qgap200_rgap200")
SV_BASE = ROOT / "sv_extraction" / "decompose" / SUFFIX
OUT = ROOT / "analyses" / "imputability" / "output"
PANEL = ROOT / "resources" / "1000G" / "20130606_g1k_3202_samples_ped_population.txt"
CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

pan = pd.read_csv(PANEL, sep=r"\s+", engine="python")
panel_map = dict(zip(pan.SampleID, pan.Superpopulation))

everyone = set()
for c in CHROMS:
    f = SV_BASE / c / ("%s_sv_samples.txt" % c)
    if f.exists():
        everyone |= set(f.read_text().split())
NO_CODE = {"HG002"}
resolved, _basis, _hand, _pp, why = tracks.resolve_ancestry(sorted(everyone - NO_CODE))
if resolved is None:
    raise SystemExit("could not resolve the SV cohort: %s" % why)
print("SV cohort across all chromosomes: %d individuals" % len(everyone))
print("  in the 3,202 panel   : %d" % sum(1 for s in everyone if s in panel_map))
print("  resolvable elsewhere : %d  (%s)"
      % (sum(1 for s in everyone if s not in panel_map and s in resolved),
         ", ".join(sorted(s for s in everyone if s not in panel_map and s in resolved))))
print("  no population code   : %s" % ", ".join(sorted(everyone & NO_CODE)))


def features(chrom, pop_of):
    npz, sam = (SV_BASE / chrom / ("%s_sv_genotypes.npz" % chrom),
                SV_BASE / chrom / ("%s_sv_samples.txt" % chrom))
    if not npz.exists() or not sam.exists():
        return None
    g = np.load(npz)
    arr = g[g.files[0]].astype(np.float32)
    cols = sam.read_text().splitlines()
    if len(cols) != arr.shape[1]:
        raise SystemExit("%s: %d sample labels for %d matrix columns"
                         % (chrom, len(cols), arr.shape[1]))
    pop_col = np.array([pop_of.get(s, "OTHER") for s in cols])
    idx = {p: np.where(pop_col == p)[0] for p in POPS}
    if min(len(v) for v in idx.values()) < 5:
        return None
    arr_m = np.where(np.isfinite(arr) & (arr >= 0), arr, np.nan)
    norm = 2.0 if np.nanmax(arr_m) > 1.05 else 1.0
    p = {}
    for pop in POPS:
        sub = arr_m[:, idx[pop]]
        n_obs = np.sum(~np.isnan(sub), axis=1)
        p[pop] = np.where(n_obs > 0, np.nansum(sub, axis=1) / (norm * n_obs), np.nan)
    pair = {}
    for i, j in combinations(POPS, 2):
        pi, pj = p[i], p[j]
        ni, nj = max(len(idx[i]), 2), max(len(idx[j]), 2)
        num = (pi - pj) ** 2 - pi * (1 - pi) / (ni - 1) - pj * (1 - pj) / (nj - 1)
        den = pi * (1 - pj) + pj * (1 - pi)
        pair["%s_%s" % (i, j)] = np.where(den > 0, num / den, np.nan)
    with np.errstate(all="ignore"):
        mx = np.nanmax(np.vstack(list(pair.values())), axis=0)
    return pd.DataFrame({"chrom": chrom, "sv_index": np.arange(arr.shape[0]),
                         "max_FST_pair_sv": mx,
                         "n_used": sum(len(v) for v in idx.values())})


def run(pop_of, label):
    frames = [f for f in (features(c, pop_of) for c in CHROMS) if f is not None]
    df = pd.concat(frames, ignore_index=True)
    print("  %-12s %d SVs, %d columns used on chr1"
          % (label, len(df), int(df.loc[df.chrom == "chr1", "n_used"].iloc[0])))
    return df.drop(columns="n_used")


print()
old = run(panel_map, "panel-only")
new = run(resolved, "resolved")

ship = pd.read_csv(OUT / "sv_ancestry_features.tsv", sep="\t",
                   usecols=["chrom", "sv_index", "max_FST_pair_sv"])
j = ship.merge(old, on=["chrom", "sv_index"], suffixes=("_ship", "_mine"))
both = j.dropna(subset=["max_FST_pair_sv_ship", "max_FST_pair_sv_mine"])
agree = np.allclose(both.max_FST_pair_sv_ship, both.max_FST_pair_sv_mine, atol=1e-9)
print("\nCONTROL  panel-only recomputation vs the shipped sv_ancestry_features.tsv: "
      "%d rows, max |diff| %.2e -> %s"
      % (len(both), float((both.max_FST_pair_sv_ship - both.max_FST_pair_sv_mine)
                          .abs().max()), "MATCH" if agree else "DIFFERS"))
if not agree:
    raise SystemExit("the recomputation does not reproduce the shipped file, so any delta "
                     "below would be measuring my reimplementation, not the cohort")

cat = pd.read_csv(OUT / "v3" / "per_sv_features_with_genome.tsv", sep="\t",
                  usecols=["chrom", "sv_index"])
print("catalog rows: %d" % len(cat))


def stats(df, name):
    v = cat.merge(df, on=["chrom", "sv_index"]).max_FST_pair_sv.dropna()
    print("  %-12s n=%d  mean %.4f  %%>0.15 %.2f%%  max %.3f"
          % (name, len(v), v.mean(), 100.0 * (v > 0.15).mean(), v.max()))
    return v.mean(), 100.0 * (v > 0.15).mean(), v.max()


print("\nR9 catalog figures (published: mean 0.0038, 0.44%% above 0.15, max 0.574)")
o = stats(old, "panel-only")
n = stats(new, "resolved")
print("  %-12s mean %+.4f  %%>0.15 %+.2f pp  max %+.3f"
      % ("delta", n[0] - o[0], n[1] - o[1], n[2] - o[2]))
