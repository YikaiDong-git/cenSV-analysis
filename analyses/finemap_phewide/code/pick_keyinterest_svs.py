#!/usr/bin/env python3
import sys, glob, os
import numpy as np
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

BUNDLES = CENSV_ROOT + "/analyses/finemap_multi_sv/output/bundles"
OUT = sys.argv[1] if len(sys.argv) > 1 else "manifests/keyinterest_svs.tsv"

TARGETS = [
    "chr4_sv398",
    "chr5_sv213", "chr8_sv2",
    "chr16_sv465", "chr16_sv37", "chr16_sv73", "chr16_sv1017", "chr16_sv559", "chr16_sv396",
    "chr10_sv697", "chr10_sv698", "chr11_sv200", "chr11_sv400", "chr11_sv148", "chr11_sv389",
    "chr17_sv317", "chr17_sv331", "chr1_sv1325",
    "chr11_sv271", "chr12_sv721", "chr1_sv650", "chr19_sv1949", "chr11_sv121", "chr3_sv40",
]
CHROMS = sorted({t.split("_")[0] for t in TARGETS}, key=lambda c: int(c[3:]))


def load_region(npz):
    b = np.load(npz, allow_pickle=True)
    sl = [str(x) for x in b["sv_locus"]]
    cv = np.asarray(b["sv_cv_r2"], float)
    nsup = (np.asarray(b["B"], float) != 0).sum(axis=0)
    return sl, cv, nsup


def physical(name):
    return "_".join(name.split("_")[:2])


def main():
    rows = []
    for chrom in CHROMS:
        for npz in sorted(glob.glob(f"{BUNDLES}/{chrom}/*_bundle.npz")):
            sl, cv, nsup = load_region(npz)
            for j, name in enumerate(sl):
                rows.append((name, physical(name), chrom, npz, j, float(cv[j]), int(nsup[j])))
    if not rows:
        sys.exit("no bundles found")

    picks, missing, multi = [], [], []
    for t in TARGETS:
        cands = [r for r in rows if r[1] == t]
        if not cands:
            missing.append(t)
            continue
        best = max(cands, key=lambda r: r[6])
        picks.append(best)
        variants = sorted({r[0] for r in cands})
        if len(variants) > 1:
            multi.append((t, len(variants), best[0]))

    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("sv_locus\tphysical\tchrom\tregion_bundle\tcol\tcv_r2\tn_support\n")
        for name, phys, chrom, npz, j, cv, nsup in picks:
            fh.write(f"{name}\t{phys}\t{chrom}\t{npz}\t{j}\t{cv:.4f}\t{nsup}\n")

    print(f"[picked {len(picks)}/{len(TARGETS)} key-interest SVs] -> {OUT}")
    for name, phys, chrom, npz, j, cv, nsup in picks:
        print(f"  {phys:14s} cv_r2={cv:.3f} nsup={nsup:4d} col={j:4d} {os.path.basename(npz):42s} [{name}]")
    if multi:
        print(f"\n[multi-variant physical SVs (picked max-support); {len(multi)}]:")
        for t, n, chosen in multi:
            print(f"  {t}: {n} strategy variants -> chose {chosen}")
    if missing:
        print(f"\n[NOT FOUND in 6,519-SV catalog; {len(missing)}]: {', '.join(missing)}")


if __name__ == "__main__":
    main()
