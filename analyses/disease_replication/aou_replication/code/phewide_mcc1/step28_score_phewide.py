#!/usr/bin/env python3
import numpy as np, csv, os, glob
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
PH = f"{BASE}/phewide_all_handoff"
BED = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
COV = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"
SCRATCH = "/tmp/phewide/phewide_pgen"
OUTQC = f"{BASE}/aou_replication/results_vds/step28_phewide"
os.makedirs(SCRATCH, exist_ok=True); os.makedirs(OUTQC, exist_ok=True)

_CODE2A1 = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)
_LUT = np.empty((256, 4), dtype=np.float64)
for _b in range(256):
    for _k in range(4):
        _LUT[_b, _k] = _CODE2A1[(_b >> (2 * _k)) & 3]

BASECOV = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
           "sex", *[f"gPC{i}" for i in range(1, 17)]]


def read_fam_iids(path):
    out = []
    with open(path) as f:
        for line in f:
            out.append(line.split()[1])
    return out


def index_bim(path, want_pos):
    pos_map = defaultdict(list); vidx = 0
    with open(path) as f:
        for line in f:
            c = line.split()
            if c[3] in want_pos:
                pos_map[c[3]].append((vidx, c[4], c[5]))
            vidx += 1
    return pos_map, vidx


def load_weights():
    rows = []
    for f in sorted(glob.glob(f"{PH}/sv_weights_long_chr*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            r["tier"] = "signal"; rows.append(r)
    for f in sorted(glob.glob(f"{PH}/discovery_null_ext/sv_weights_long_ext_chr*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            r["tier"] = "null_ext"; rows.append(r)
    return rows


def main():
    weights = load_weights()
    tier = {r["sv_locus"]: r["tier"] for r in weights}
    by_chrom = defaultdict(lambda: defaultdict(list))
    for r in weights:
        ch = r["chrom"] if str(r["chrom"]).startswith("chr") else f"chr{r['chrom']}"
        by_chrom[ch][r["sv_locus"]].append(r)
    print(f"{len(weights):,} weight-SNPs, {len(tier):,} SVs "
          f"({sum(v=='signal' for v in tier.values())} signal / {sum(v=='null_ext' for v in tier.values())} ext)", flush=True)

    fam_iids = None
    for ch in by_chrom:
        p = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000.fam"
        if os.path.exists(p):
            fam_iids = read_fam_iids(p); break
    covset = set()
    with open(COV) as f:
        h = f.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if all(p[idx[c]] not in ("", "NA") for c in BASECOV if c in idx):
                covset.add(p[idx["IID"]])
    Nfam = len(fam_iids)
    cols = np.array([i for i, iid in enumerate(fam_iids) if iid in covset], dtype=np.int64)
    samp_iids = [fam_iids[i] for i in cols]
    Nanl = len(cols)
    print(f"VDS fam N={Nfam}; analysis cohort (covar-complete ∩ fam) N={Nanl}", flush=True)

    if not os.path.exists(f"{SCRATCH}/sv.fam"):
        sex = {}
        with open(COV) as f:
            h = f.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
            for line in f:
                p = line.rstrip("\n").split("\t")
                sex[p[idx["IID"]]] = {"Male": "1", "Female": "2"}.get(p[idx["sex"]], "0")
        with open(f"{SCRATCH}/sv.fam", "w") as f:
            for iid in samp_iids:
                f.write(f"{iid}\t{iid}\t0\t0\t{sex.get(iid,'0')}\t-9\n")
    with open(f"{SCRATCH}/sv_tier.tsv", "w") as f:
        f.write("sv_locus\ttier\n")
        for sv in sorted(tier):
            f.write(f"{sv}\t{tier[sv]}\n")

    dosef = f"{SCRATCH}/sv.dose"
    done = set()
    if os.path.exists(dosef):
        with open(dosef) as f:
            for line in f:
                done.add(line.split("\t", 1)[0])
        print(f"resume: {len(done)} SVs already scored", flush=True)

    qc = {}
    fo = open(dosef, "a")
    for ch in sorted(by_chrom):
        pref = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000"
        if not os.path.exists(pref + ".bed"):
            print(f"  {ch}: NO BED -> skipped", flush=True); continue
        assert read_fam_iids(pref + ".fam") == fam_iids, f"{ch} fam order differs"
        svrows = by_chrom[ch]
        want = {r["pos"] for sv in svrows for r in svrows[sv]}
        pos_map, nvar = index_bim(pref + ".bim", want)
        bpv = (Nfam + 3) // 4
        nsv_done = 0
        with open(pref + ".bed", "rb") as bf:
            assert bf.read(3) == b"\x6c\x1b\x01", f"{ch}: bad .bed magic"
            for sv, rws in svrows.items():
                if sv in done:
                    continue
                acc = np.zeros(Nfam, dtype=np.float64); n_used = 0; n_sw = 0
                for row in rws:
                    hit = None
                    for (vi, a1, a2) in pos_map.get(row["pos"], []):
                        if a1 == row["alt"] and a2 == row["ref"]:
                            hit = (vi, False); break
                        if a1 == row["ref"] and a2 == row["alt"]:
                            hit = (vi, True); break
                    if hit is None:
                        continue
                    vi, sw = hit
                    bf.seek(3 + vi * bpv)
                    a1c = _LUT[np.frombuffer(bf.read(bpv), dtype=np.uint8)].reshape(-1)[:Nfam]
                    g = (2.0 - a1c) if sw else a1c
                    m = np.nanmean(g); s = np.nanstd(g)
                    if not s > 0:
                        continue
                    g = np.where(np.isnan(g), m, g)
                    acc += float(row["beta_std"]) * (g - m) / s
                    n_used += 1; n_sw += int(sw)
                d = acc[cols]
                lo, hi = float(d.min()), float(d.max())
                ds = 2.0 * (d - lo) / (hi - lo) if hi > lo else np.full(Nanl, 1.0)
                fo.write(f"{sv}\tA\tC\t" + "\t".join(f"{x:.5f}" for x in ds) + "\n"); fo.flush()
                qc[sv] = (ch, len(rws), n_used, n_sw, float(d.mean()), float(d.std()))
                nsv_done += 1
        print(f"  {ch}: scored {nsv_done} new SVs (of {len(svrows)}); bim {nvar} vars", flush=True)
        del pos_map
    fo.close()

    qcf = f"{OUTQC}/dosage_qc_phewide.tsv"
    prior = {}
    if os.path.exists(qcf):
        for r in csv.DictReader(open(qcf), delimiter="\t"):
            prior[r["sv_locus"]] = r
    with open(qcf, "w") as f:
        f.write("sv_locus\ttier\tchrom\tn_total\tn_used\tn_swapped\tdosage_mean\tdosage_sd\n")
        for sv in sorted(tier):
            if sv in qc:
                ch, nt, nu, nsw, mn, sd = qc[sv]
                f.write(f"{sv}\t{tier[sv]}\t{ch}\t{nt}\t{nu}\t{nsw}\t{mn:.5g}\t{sd:.5g}\n")
            elif sv in prior:
                p = prior[sv]
                f.write(f"{sv}\t{tier[sv]}\t{p['chrom']}\t{p['n_total']}\t{p['n_used']}\t{p['n_swapped']}\t{p['dosage_mean']}\t{p['dosage_sd']}\n")
    print(f"done. sv.dose rows now = {len(done)+len(qc)} ; QC -> {qcf}", flush=True)


if __name__ == "__main__":
    main()
