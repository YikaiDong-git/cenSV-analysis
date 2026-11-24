#!/usr/bin/env python3
import sys, glob, os
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BUNDLES = CENSV_ROOT + "/analyses/finemap_multi_sv/output/bundles"
OUT = sys.argv[1] if len(sys.argv) > 1 else "manifests/smoke_svs.tsv"

ANCHORS = ["chr4_sv349", "chr11_sv263", "chr3_sv20", "chr5_sv324", "chr8_sv5", "chr16_sv1029"]
CHROMS = ["chr1", "chr3", "chr4", "chr5", "chr8", "chr11", "chr16", "chr20"]


def load_region(npz):
    b = np.load(npz, allow_pickle=True)
    sl = [str(x) for x in b["sv_locus"]]
    cv = np.asarray(b["sv_cv_r2"], float)
    nsup = (np.asarray(b["B"], float) != 0).sum(axis=0)
    p = int(len(b["snp_pos"]))
    return sl, cv, nsup, p


def physical(name):
    parts = name.split("_")
    return "_".join(parts[:2])


def main():
    rows = []
    seen_phys = set()
    picks = []
    region_cache = {}
    for chrom in CHROMS:
        npzs = sorted(glob.glob(f"{BUNDLES}/{chrom}/*_bundle.npz"))
        if not npzs:
            continue
        for npz in npzs:
            sl, cv, nsup, p = load_region(npz)
            region_cache[npz] = (sl, cv, nsup, p)
            for j, name in enumerate(sl):
                rows.append((name, physical(name), chrom, npz, j, float(cv[j]), int(nsup[j])))
    if not rows:
        sys.exit("no bundles found")

    for anc in ANCHORS:
        cands = [r for r in rows if r[1] == anc]
        if cands:
            best = max(cands, key=lambda r: r[6])
            if best[1] not in seen_phys:
                picks.append(best); seen_phys.add(best[1])

    remaining = [r for r in rows if r[1] not in seen_phys and r[6] >= 5]
    remaining.sort(key=lambda r: r[5])
    if remaining:
        idxs = np.linspace(0, len(remaining) - 1, num=max(0, 10 - len(picks)) + 4).astype(int)
        for i in idxs:
            r = remaining[i]
            if r[1] in seen_phys:
                continue
            picks.append(r); seen_phys.add(r[1])
            if len(picks) >= 10:
                break

    picks = picks[:10]
    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("sv_locus\tphysical\tchrom\tregion_bundle\tcol\tcv_r2\tn_support\n")
        for name, phys, chrom, npz, j, cv, nsup in picks:
            fh.write(f"{name}\t{phys}\t{chrom}\t{npz}\t{j}\t{cv:.4f}\t{nsup}\n")
    print(f"[picked {len(picks)} smoke SVs] -> {OUT}")
    for name, phys, chrom, npz, j, cv, nsup in picks:
        print(f"  {phys:14s} cv_r2={cv:.3f} nsup={nsup:4d} col={j:4d} {os.path.basename(npz)}")


if __name__ == "__main__":
    main()
