#!/usr/bin/env python
"""build_phewide_null_ext.py -- the discovery-NULL EXTENSION of the phewide handoff."""
import numpy as np, os, csv, glob
from collections import defaultdict
import os as _censv_os
CENSV_WORKSPACE = _censv_os.environ["CENSV_WORKSPACE"]

P = CENSV_WORKSPACE + ""
SSIMP = f"{P}/Project1_Centromere/analyses/finemap_multi_sv"
AGG = f"{SSIMP}/validation/ssimp_impg/work/phewide_agg"
BUND = f"{SSIMP}/output/bundles"
OUT = f"{P}/Project1_Centromere/analyses/disease_replication/phewide_all_handoff/discovery_null_ext"
os.makedirs(OUT, exist_ok=True)
FIT_PANEL = "all_1kGP_n185"
MAXROWS = 45000

signal = set(r["sv_locus"] for r in csv.DictReader(open(f"{AGG}/phewide_pip_ge_0.01.tsv"), delimiter="\t"))
print(f"signal sv_locus (excluded from ext): {len(signal)}", flush=True)

wrows = defaultdict(list)
mafd = defaultdict(list)
meta_rows = []
seen = set()
for npz in sorted(glob.glob(f"{BUND}/*/*_bundle.npz")):
    region = os.path.basename(npz).replace("_bundle.npz", "")
    chrom = region.split("_")[0].replace("chr", "")
    d = np.load(npz, allow_pickle=True)
    svn = np.array([str(x) for x in d["sv_locus"]])
    pos = d["snp_pos"].astype(int)
    ref = np.array([str(x).upper() for x in d["snp_ref"]])
    alt = np.array([str(x).upper() for x in d["snp_alt"]])
    mafp = d["snp_maf_panel"].astype(float)
    schr = d["snp_chrom"]
    cvr2 = d["sv_cv_r2"]
    metf = npz.replace("_bundle.npz", "_meta.tsv")
    metmap = {}
    if os.path.exists(metf):
        for r in csv.DictReader(open(metf), delimiter="\t"):
            metmap[r["sv_locus"]] = r
    for ti, sv in enumerate(svn):
        if sv in signal or sv in seen:
            continue
        seen.add(sv)
        w = d["B"][:, ti].astype(float)
        sup = np.where(w != 0)[0]
        for j in sup:
            cc = str(schr[j])
            cc = cc if cc.startswith("chr") else f"chr{cc}"
            wrows[cc].append((sv, cc, int(pos[j]), ref[j], alt[j], alt[j], f"{float(w[j]):.10g}"))
            mafd[cc].append((sv, cc, int(pos[j]), ref[j], alt[j], f"{float(mafp[j]):.6f}"))
        m = metmap.get(sv, {})
        meta_rows.append([sv, region, chrom, f"{float(cvr2[ti]):.3f}", len(sup),
                          m.get("sv_category", ""), m.get("strategy", ""), m.get("arm", ""), "none", FIT_PANEL])


def write_split(prefix, header, by):
    files = []
    for cc in sorted(by, key=lambda c: (len(c), c)):
        rr = sorted(by[cc], key=lambda x: (x[0], x[2]))
        nparts = (len(rr) + MAXROWS - 1) // MAXROWS
        for k in range(nparts):
            ch = rr[k * MAXROWS:(k + 1) * MAXROWS]
            suf = "" if nparts == 1 else f"_part{k + 1}"
            fn = f"{OUT}/{prefix}_{cc}{suf}.tsv"
            with open(fn, "w") as fh:
                fh.write("\t".join(header) + "\n")
                for row in ch:
                    fh.write("\t".join(map(str, row)) + "\n")
            files.append((os.path.basename(fn), len(ch), os.path.getsize(fn)))
    return files


wf = write_split("sv_weights_long_ext", ["sv_locus", "chrom", "pos", "ref", "alt", "a1_effect", "beta_std"], wrows)
mf = write_split("train_maf_ext", ["sv_locus", "chrom", "pos", "ref", "alt", "train_maf_panel"], mafd)
meta_rows.sort(key=lambda r: (r[2], r[0]))
with open(f"{OUT}/sv_metadata_ext.tsv", "w") as fh:
    fh.write("sv_locus\tregion\tchrom\tcv_r2\tn_weight_snps\tsv_category\tsv_strategy\tsv_arm\tdiscovery_prior\tfit_panel\n")
    for r in meta_rows:
        fh.write("\t".join(map(str, r)) + "\n")

maxsz = max([s for _, _, s in wf + mf], default=0)
print(f"null_ext sv_locus: {len(meta_rows)}  weight_files: {len(wf)}  maf_files: {len(mf)}  "
      f"max_file_bytes: {maxsz}  under_4.9MB: {maxsz < 4.9e6}", flush=True)
print("DONE_NULL_EXT", flush=True)
