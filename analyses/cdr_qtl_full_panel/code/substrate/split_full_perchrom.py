#!/usr/bin/env python
import pandas as pd, numpy as np, os
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

WU = CENSV_ROOT + "/analyses/cdr_qtl_full_panel"
GEN = f"{WU}/geno"
PC = f"{GEN}/perchrom"
PH_CORR = CENSV_ROOT + "/analyses/cdr_qtl/output/per_hap_chrom_phenotypes.tsv"
PH_DEEP = CENSV_ROOT + "/analyses/cdr_qtl/output/per_hap_chrom_phenotypes_deepest.tsv"
MIN_HAPS = 20

G = np.load(f"{GEN}/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{GEN}/full_variants.tsv", sep="\t")
samples = open(f"{GEN}/full_samples.txt").read().split()
assert G.shape == (len(vt), len(samples)), (G.shape, len(vt), len(samples))

nmiss = int((G < 0).sum())
for i in range(G.shape[0]):
    row = G[i, :]
    m = row < 0
    if m.any():
        obs = row[~m]
        row[m] = obs.mean() if len(obs) else 0.0
        G[i, :] = row
print(f"imputed {nmiss} missing cells (locus-mean) over {G.shape[0]} loci")

for ch, sub in vt.groupby("chrom"):
    idx = sub.index.values
    d = f"{PC}/{ch}"
    os.makedirs(d, exist_ok=True)
    np.savez_compressed(f"{d}/{ch}_sv_genotypes.npz", genotypes=G[idx, :])
    open(f"{d}/{ch}_sv_samples.txt", "w").write("\n".join(samples) + "\n")
    sub.reset_index(drop=True).to_csv(f"{d}/{ch}_sv_variants.tsv", sep="\t", index=False)
print(f"wrote per-chrom files for {vt.chrom.nunique()} chroms under {PC}")

pc = pd.read_csv(PH_CORR, sep="\t")
pd_ = pd.read_csv(PH_DEEP, sep="\t")
avail = {"cdr_pos_offset": pc.groupby("chm13_chrom").cdr_pos_offset.count().to_dict(),
         "cdr_pos_offset_deepest": pd_.groupby("chm13_chrom").cdr_pos_offset_deepest.count().to_dict()}
tasks, skipped = [], []
for ch in sorted(vt.chrom.unique(), key=lambda c: (len(c), c)):
    for ph in ["cdr_pos_offset", "cdr_pos_offset_deepest"]:
        n = avail[ph].get(ch, 0)
        (tasks if n >= MIN_HAPS else skipped).append((ch, ph, n))
with open(f"{WU}/code/scan/scan_tasks.txt", "w") as fh:
    for ch, ph, n in tasks:
        fh.write(f"{ch}\t{ph}\n")
print(f"\nscan tasks: {len(tasks)}  (skipped {len(skipped)} for <{MIN_HAPS} phenotyped haps: {skipped})")
print("per-chrom loci in the scan:")
print(vt.chrom.value_counts().reindex(sorted(vt.chrom.unique(), key=lambda c: (len(c), c))).to_string())
