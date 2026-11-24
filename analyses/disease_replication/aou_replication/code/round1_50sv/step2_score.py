#!/usr/bin/env python3
"""Step 2 of the AoU imputed-SV replication: per-person imputed SV dosage."""
import pgenlib, numpy as np, csv, os, gzip

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
H = f"{BASE}/aou_handoff"
PGEN = f"{HOME}/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen"
OUT = f"{BASE}/aou_replication/results/step2"
os.makedirs(OUT, exist_ok=True)

weights = list(csv.DictReader(open(f"{H}/sv_weights_long.tsv"), delimiter="\t"))
by_chrom = {}
for r in weights:
    by_chrom.setdefault(r["chrom"], []).append(r)
sv_order = sorted({r["sv_locus"] for r in weights})

def read_iids(chrom):
    iids = []
    with open(f"{PGEN}/acaf_threshold.{chrom}.psam") as f:
        header = f.readline().rstrip("\n").split("\t")
        iid_col = header.index("#IID") if "#IID" in header else (header.index("IID") if "IID" in header else 1)
        for line in f:
            iids.append(line.rstrip("\n").split("\t")[iid_col])
    return iids

iids = read_iids("chr1")
N = len(iids)
print(f"samples N={N}", flush=True)

dosage = {sv: np.zeros(N, dtype=np.float32) for sv in sv_order}
qc = {sv: {"chrom": "", "n_total": 0, "n_used": 0, "n_swapped": 0} for sv in sv_order}
for r in weights:
    qc[r["sv_locus"]]["n_total"] += 1
    qc[r["sv_locus"]]["chrom"] = r["chrom"]

for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
    wr = by_chrom[chrom]
    want = {r["pos"] for r in wr}
    offsets = [0]; pos_map = {}; vidx = 0
    with open(f"{PGEN}/acaf_threshold.{chrom}.pvar") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t", 5)
            alts = c[4].split(",")
            offsets.append(offsets[-1] + 1 + len(alts))
            if c[1] in want:
                pos_map.setdefault(c[1], []).append((vidx, c[3], alts))
            vidx += 1
    offsets = np.array(offsets, dtype=np.uintp)
    assert read_iids(chrom) == iids, f"{chrom} psam sample ORDER differs from chr1"
    reader = pgenlib.PgenReader(
        f"{PGEN}/acaf_threshold.{chrom}.pgen".encode(),
        raw_sample_ct=N, variant_ct=vidx, allele_idx_offsets=offsets,
    )
    buf = np.empty(N, dtype=np.int32)
    used = 0
    for row in wr:
        pos, ref, alt, beta, sv = row["pos"], row["ref"], row["alt"], float(row["beta_std"]), row["sv_locus"]
        hit = None
        for (vi, pref, palts) in pos_map.get(pos, []):
            if pref == ref and alt in palts:
                hit = (vi, 1 + palts.index(alt), False); break
            if pref == alt and ref in palts:
                hit = (vi, 0, True); break
        if hit is None:
            continue
        vi, ai, sw = hit
        reader.read(vi, buf, allele_idx=ai)
        g = buf.astype(np.float64)
        g[g < 0] = np.nan
        m = np.nanmean(g); s = np.nanstd(g)
        if not s > 0:
            continue
        g = np.where(np.isnan(g), m, g)
        dosage[sv] += (beta * (g - m) / s).astype(np.float32)
        qc[sv]["n_used"] += 1; qc[sv]["n_swapped"] += int(sw); used += 1
    reader.close()
    del offsets, pos_map
    print(f"  {chrom}: {used} SNPs scored", flush=True)

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
