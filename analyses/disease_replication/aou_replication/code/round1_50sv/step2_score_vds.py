#!/usr/bin/env python3
"""Step 2 (VDS re-score): per-person imputed SV dosage from Josh's VDS-extracted pgens."""
import numpy as np, csv, os, gzip

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
H = f"{BASE}/aou_handoff"
BED = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
OUT = f"{BASE}/aou_replication/results_vds/step2"
os.makedirs(OUT, exist_ok=True)

WEIGHTS = f"{H}/sv_weights_long.tsv"

_CODE2A1 = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)
_LUT = np.empty((256, 4), dtype=np.float64)
for b in range(256):
    for k in range(4):
        _LUT[b, k] = _CODE2A1[(b >> (2 * k)) & 3]


def read_fam_iids(path):
    iids = []
    with open(path) as f:
        for line in f:
            p = line.split()
            iids.append(p[1])
    return iids


def index_bim(path, want_pos):
    pos_map = {}
    vidx = 0
    with open(path) as f:
        for line in f:
            c = line.split()
            if c[3] in want_pos:
                pos_map.setdefault(c[3], []).append((vidx, c[4], c[5]))
            vidx += 1
    return pos_map, vidx


def main():
    weights = list(csv.DictReader(open(WEIGHTS), delimiter="\t"))
    by_chrom = {}
    for r in weights:
        by_chrom.setdefault(r["chrom"], []).append(r)
    sv_order = sorted({r["sv_locus"] for r in weights})

    iids = None
    N = None
    dosage = {sv: None for sv in sv_order}
    qc = {sv: {"chrom": "", "n_total": 0, "n_used": 0, "n_swapped": 0} for sv in sv_order}
    for r in weights:
        qc[r["sv_locus"]]["n_total"] += 1
        qc[r["sv_locus"]]["chrom"] = r["chrom"]

    for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
        pref = f"{BED}/{chrom}/{chrom}.pericentromeric_snps.min_ac1000"
        if not os.path.exists(pref + ".bed"):
            print(f"  {chrom}: NO BED ({pref}.bed) -> skipped", flush=True)
            continue
        fam_iids = read_fam_iids(pref + ".fam")
        if iids is None:
            iids = fam_iids
            N = len(iids)
            dosage = {sv: np.zeros(N, dtype=np.float32) for sv in sv_order}
            print(f"samples N={N} (Josh fam order)", flush=True)
        else:
            assert fam_iids == iids, f"{chrom} fam sample order differs from first chrom"

        wr = by_chrom[chrom]
        want = {r["pos"] for r in wr}
        pos_map, nvar = index_bim(pref + ".bim", want)

        bytes_per_var = (N + 3) // 4
        with open(pref + ".bed", "rb") as bf:
            magic = bf.read(3)
            assert magic == b"\x6c\x1b\x01", f"{chrom}: bad .bed magic {magic!r} (need variant-major)"
            used = 0
            for row in wr:
                pos, ref, alt = row["pos"], row["ref"], row["alt"]
                beta, sv = float(row["beta_std"]), row["sv_locus"]
                hit = None
                for (vi, a1, a2) in pos_map.get(pos, []):
                    if a1 == alt and a2 == ref:
                        hit = (vi, False); break
                    if a1 == ref and a2 == alt:
                        hit = (vi, True); break
                if hit is None:
                    continue
                vi, sw = hit
                bf.seek(3 + vi * bytes_per_var)
                raw = np.frombuffer(bf.read(bytes_per_var), dtype=np.uint8)
                a1c = _LUT[raw].reshape(-1)[:N]
                g = (2.0 - a1c) if sw else a1c
                m = np.nanmean(g); s = np.nanstd(g)
                if not s > 0:
                    continue
                g = np.where(np.isnan(g), m, g)
                dosage[sv] += (beta * (g - m) / s).astype(np.float32)
                qc[sv]["n_used"] += 1; qc[sv]["n_swapped"] += int(sw); used += 1
        del pos_map
        print(f"  {chrom}: {used} SNPs scored (of {len(wr)} weights; bim had {nvar} variants)", flush=True)

    dpath = f"{OUT}/sv_dosages.tsv.gz"
    with gzip.open(dpath, "wt") as f:
        f.write("IID\t" + "\t".join(sv_order) + "\n")
        for i in range(N):
            f.write(iids[i] + "\t" + "\t".join(f"{dosage[sv][i]:.6g}" for sv in sv_order) + "\n")

    qpath = f"{OUT}/dosage_qc.tsv"
    with open(qpath, "w") as f:
        f.write("sv_locus\tchrom\tn_total\tn_used\tn_swapped\tdosage_mean\tdosage_sd\n")
        for sv in sv_order:
            d = dosage[sv]
            f.write(f"{sv}\t{qc[sv]['chrom']}\t{qc[sv]['n_total']}\t{qc[sv]['n_used']}\t{qc[sv]['n_swapped']}\t{d.mean():.4g}\t{d.std():.4g}\n")

    nz = sum(1 for sv in sv_order if qc[sv]["n_used"] > 0)
    print(f"\nDONE: {nz}/{len(sv_order)} SVs scored (n_used>0)")
    print(f"  dosages (individual-level, in-perimeter): {dpath}")
    print(f"  QC (aggregate, shareable): {qpath}")


if __name__ == "__main__":
    main()
