#!/usr/bin/env python3
"""step97b_dosage.py -- LANE D of the step97 discovery-request package."""

import csv
import gc
import os
import sys

import numpy as np
import pgenlib

R = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
     "disease_replication/aou_replication")
SRC = f"{R}/results_vds/step28_phewide_ext2"
PGDIR = f"{SRC}/pgen"
COVAR = f"{R}/results_vds/phewide_carryover/inputs/covar.txt"
HITS = f"{R}/results_vds/step30_phewide_ext2/disease_matrices_firth/disease_hits_final.tsv"
OUT = f"{R}/results_vds/step97_discovery_requests"

EXTRA_SV = "chr8_sv30_DEL_mon_best8_midpoint1000"
STRATA = ["EUR", "AFR", "AMR", "EAS", "MID", "SAS"]
MIN_CELL = 20
NEAR_CONST_K = 0.01
LOW_INFO_VR = 0.10


def die(msg):
    sys.exit(f"ABORT: {msg}")


def manifest_prefixes():
    rows = list(csv.DictReader(open(f"{SRC}/pgen_manifest.tsv"), delimiter="\t"))
    out = []
    for r in rows:
        p = f"{PGDIR}/{r['prefix']}"
        if not os.path.exists(p + ".pgen"):
            die(f"manifest lists a missing pgen: {p}.pgen")
        out.append((p, int(r["n_variants"])))
    return out


def build_variant_index(pre):
    loc = {}
    sig0 = None
    for pfx, nv in pre:
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        if len(ids) != nv:
            die(f"{pfx}: pvar has {len(ids)} variants, manifest says {nv}")
        for i, sv in enumerate(ids):
            if sv in loc:
                die(f"duplicate sv_locus across chunks: {sv}")
            loc[sv] = (pfx, i)
        iids = [l.split("\t")[1] for l in open(f"{pfx}.psam") if not l.startswith("#")]
        s = "\n".join(iids)
        if sig0 is None:
            sig0 = s
        elif s != sig0:
            die(f"psam sample order mismatch at {pfx}.psam -- cannot align to covar.txt")
    return loc, sig0.split("\n")


def load_strata(psam_iids):
    hdr = open(COVAR).readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    need = ["anc_amr", "anc_eas", "anc_eur", "anc_mid", "anc_sas"]
    for c in need:
        if c not in idx:
            die(f"covar.txt is missing column {c}")
    i_iid = idx["IID"]
    code = np.full(len(psam_iids), -1, np.int8)
    onehot_to_stratum = {
        "anc_eur": "EUR", "anc_amr": "AMR", "anc_eas": "EAS",
        "anc_mid": "MID", "anc_sas": "SAS",
    }
    sidx = {s: k for k, s in enumerate(STRATA)}
    n = 0
    with open(COVAR) as fh:
        fh.readline()
        for row, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if row >= len(psam_iids):
                die("covar.txt has more data rows than the psam")
            if f[i_iid] != psam_iids[row]:
                die(f"covar.txt row {row} IID does not match psam row {row}; "
                    "positional alignment is unsafe")
            hot = [c for c in need if float(f[idx[c]]) != 0.0]
            if len(hot) > 1:
                die(f"row {row}: more than one anc_* one-hot set ({hot})")
            code[row] = sidx["AFR"] if not hot else sidx[onehot_to_stratum[hot[0]]]
            n += 1
    if n != len(psam_iids):
        die(f"covar.txt has {n} data rows, psam has {len(psam_iids)}")
    if (code < 0).any():
        die("some participants were not assigned a stratum")
    return code


def load_targets():
    q05, q10 = set(), set()
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        q = float(r["q_bh_tier"])
        if q < 0.05:
            q05.add(r["sv_locus"])
        elif q < 0.10:
            q10.add(r["sv_locus"])
    q10 -= q05
    band = {s: "q_lt_0.05" for s in q05}
    band.update({s: "q_0.05_to_0.10" for s in q10})
    if EXTRA_SV not in band:
        band[EXTRA_SV] = "reference_not_a_hit"
    return band, len(q05), len(q10)


def load_sv_meta():
    cv = {}
    for r in csv.DictReader(open(f"{SRC}/sv_cv_r2_all.tsv"), delimiter="\t"):
        cv[r["sv_locus"]] = r["cv_r2"]
    tier, nsnp, unsc = {}, {}, {}
    for r in csv.DictReader(open(f"{SRC}/dosage_qc_ext2.tsv"), delimiter="\t"):
        tier[r["sv_locus"]] = r["tier"]
        nsnp[r["sv_locus"]] = (r.get("n_used", "NA"), r.get("n_total", "NA"))
        unsc[r["sv_locus"]] = r.get("unscoreable_mass_frac", "NA")
    return cv, tier, nsnp, unsc


def fmt(x, nd=6):
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def main():
    os.makedirs(OUT, exist_ok=True)
    pre = manifest_prefixes()
    loc, psam_iids = build_variant_index(pre)
    N = len(psam_iids)
    print(f"[i] {len(loc)} SVs across {len(pre)} pgen chunks; N={N}", flush=True)

    code = load_strata(psam_iids)
    del psam_iids
    gc.collect()

    sidx = [np.flatnonzero(code == k).astype(np.int32) for k in range(len(STRATA))]
    ns = [len(a) for a in sidx]
    print("[i] stratum N: " + ", ".join(f"{s}={n}" for s, n in zip(STRATA, ns)), flush=True)
    expected = {"EUR": 177230, "AFR": 64726, "AMR": 60984,
                "EAS": 7144, "SAS": 3944, "MID": 1176}
    for s, n in zip(STRATA, ns):
        if n != expected[s]:
            die(f"stratum {s}: derived N={n} but the arm's documented N={expected[s]}")
    del code
    gc.collect()

    band, n05, n10 = load_targets()
    cv, tier, nsnp, unsc = load_sv_meta()
    print(f"[i] targets: {n05} SVs q<0.05, {n10} SVs 0.05<=q<0.10, "
          f"+1 reference = {len(band)}", flush=True)

    missing = sorted(s for s in band if s not in loc)
    if missing:
        die(f"target SVs absent from the archived pgen: {missing}")

    by_chunk = {}
    for sv in band:
        pfx, vi = loc[sv]
        by_chunk.setdefault(pfx, []).append((vi, sv))

    buf = np.empty(N, np.float64)
    rows = []
    thr_rows = []
    THRESH = [0.25, 0.5, 1.0, 1.5]
    n_done = 0
    for pfx in sorted(by_chunk):
        rd = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for vi, sv in sorted(by_chunk[pfx]):
            rd.read_dosages(vi, buf)
            nonfinite = int((~np.isfinite(buf)).sum())
            if nonfinite:
                die(f"{sv}: {nonfinite} non-finite dosages -- the 100%-call-rate "
                    "assumption is violated, do not silently impute")
            a_mean = float(buf.mean())
            a_sd = float(buf.std(ddof=1))
            a_var = a_sd * a_sd
            a_lo, a_hi = float(buf.min()), float(buf.max())
            a_range = a_hi - a_lo
            wvar_acc = [0.0] * len(STRATA)

            for k, s in enumerate(STRATA):
                g = buf[sidx[k]]
                n = int(g.size)
                mean = float(g.mean())
                sd = float(g.std(ddof=1)) if n > 1 else float("nan")
                var = sd * sd if np.isfinite(sd) else float("nan")
                vr = (var / a_var) if a_var > 0 else float("nan")
                span = ((float(g.max()) - float(g.min())) / a_range * 100.0
                        if a_range > 0 else float("nan"))
                nc = bool(np.isfinite(sd) and a_sd > 0 and sd < NEAR_CONST_K * a_sd)
                li = bool(np.isfinite(vr) and vr < LOW_INFO_VR)
                relinf = (n * var) / (N * a_var) if a_var > 0 else float("nan")
                wvar_acc[k] += (n - 1) * var if n > 1 else 0.0
                rows.append([
                    sv, s,
                    (str(n) if n >= MIN_CELL else f"<{MIN_CELL}"),
                    fmt(mean), fmt(sd), fmt(var, 8), fmt(vr, 6), fmt(span, 4),
                    "TRUE" if nc else "FALSE", "TRUE" if li else "FALSE",
                    fmt(relinf, 6),
                    fmt(a_mean), fmt(a_sd), fmt(a_var, 8),
                    band[sv], tier.get(sv, "NA"), cv.get(sv, "NA"),
                    nsnp.get(sv, ("NA", "NA"))[0], nsnp.get(sv, ("NA", "NA"))[1],
                    unsc.get(sv, "NA"),
                ])
                del g

            if a_var > 0:
                a_var_pop = a_var * (N - 1) / N
                within = sum(wvar_acc) / N
                between = sum(
                    len(sidx[k]) * (float(buf[sidx[k]].mean()) - a_mean) ** 2
                    for k in range(len(STRATA))) / N
                if abs((within + between) / a_var_pop - 1.0) > 1e-8:
                    die(f"{sv}: variance decomposition failed "
                        f"(within+between)/total = {(within+between)/a_var_pop:.9f}")

            counts = []
            for t in THRESH:
                c = int((buf > t).sum())
                counts.append(str(c) if c >= MIN_CELL else f"<{MIN_CELL}")
            thr_rows.append([sv, band[sv], fmt(a_mean), fmt(a_sd)] + counts)

            n_done += 1
            if n_done % 20 == 0:
                print(f"  {n_done}/{len(band)} SVs", flush=True)
        rd.close()
        gc.collect()

    hdr = ["sv_locus", "stratum", "n_participants", "dosage_mean", "dosage_sd",
           "dosage_var", "var_ratio_vs_allcohort", "pct_dosage_range_spanned",
           "near_constant", "low_information",
           "frac_cohort_dosage_information",
           "allcohort_dosage_mean", "allcohort_dosage_sd",
           "allcohort_dosage_var", "q_band", "tier",
           "cv_r2_global_NOT_per_stratum", "n_tag_snps_used", "n_tag_snps_total",
           "unscoreable_weight_mass_frac"]
    rows.sort(key=lambda r: (r[0], STRATA.index(r[1])))
    p1 = f"{OUT}/perstratum_dosage_denominators.tsv"
    with open(p1, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    p2 = f"{OUT}/_threshold_sensitivity.tsv"
    with open(p2, "w") as f:
        f.write("sv_locus\tq_band\tallcohort_dosage_mean\tallcohort_dosage_sd\t"
                + "\t".join(f"n_dosage_gt_{t}" for t in THRESH) + "\n")
        for r in sorted(thr_rows):
            f.write("\t".join(map(str, r)) + "\n")

    for p in (p1, p2):
        mb = os.path.getsize(p) / 1e6
        print(f"[w] {p}  ({mb:.2f} MB, {'OK' if mb < 4.5 else 'TOO BIG'})", flush=True)
        if mb >= 4.5:
            die(f"{p} exceeds the 4.5 MB cap")

    nnc = sum(1 for r in rows if r[8] == "TRUE")
    nli = sum(1 for r in rows if r[9] == "TRUE")
    print(f"[i] {len(rows)} SV x stratum cells; "
          f"near_constant(sd<{NEAR_CONST_K}*allcohort_sd) = {nnc}; "
          f"low_information(var_ratio<{LOW_INFO_VR}) = {nli}", flush=True)


if __name__ == "__main__":
    main()
