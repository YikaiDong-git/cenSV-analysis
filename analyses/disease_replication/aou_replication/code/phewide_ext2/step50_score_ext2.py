#!/usr/bin/env python3
"""Step 50 (ext2 arm) — score ALL 11,129 ext2 predictors per-person from the VDS .bed into a NEW pgen tree."""
import os, sys, csv, glob, json, time, argparse, subprocess, shutil
import numpy as np
from collections import defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phewide_ext"))
import step42_gate2_headline_rescore as G

HOME = os.path.expanduser("~")
BASE = G.BASE
PH = G.PH
EXT2 = G.EXT2
BED = G.BED
RV = G.RV
OUT = f"{RV}/step28_phewide_ext2"
PGDIR = f"{OUT}/pgen"
SCRATCH = "/tmp/phewide_ext2"
PLINK2 = f"{HOME}/.local/bin/plink2"
CARRY = G.CARRY

KNOWN_EMPTY_NOTE = "VDS .bim has 0 variants (extraction gap, not a weights problem)"


def chrom_order(ch):
    b = ch.replace("chr", "")
    return 99 if b == "X" else (100 if b == "Y" else int(b))


def preflight():
    info = {}
    for ch in sorted({f"chr{c}" for c in list(range(1, 23))} | {"chrX"}, key=chrom_order):
        pref = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000"
        if not os.path.exists(pref + ".bed"):
            info[ch] = dict(ok=False, nvar=0, reason="no .bed"); continue
        nb = os.path.getsize(pref + ".bim")
        nvar = sum(1 for _ in open(pref + ".bim")) if nb else 0
        info[ch] = dict(ok=nvar > 0, nvar=nvar, bed_bytes=os.path.getsize(pref + ".bed"),
                        reason="" if nvar > 0 else KNOWN_EMPTY_NOTE)
    return info


def build_tier_map(svs):
    old = {}
    with open(f"{CARRY}/sv_tier.tsv") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            old[p[0]] = p[1]
    tier = {}
    for sv, ch in svs.items():
        if sv in old:
            tier[sv] = old[sv]
        elif ch == "chrX":
            tier[sv] = "ext2_chrX"
        else:
            tier[sv] = "ext2_new"
    assert len(tier) == len(svs), "tier map incomplete"
    return tier, old


def score_chrom(task):
    ch, blk, svlist, fam_iids, cols, chunk_size = task
    Nfam = len(fam_iids); Nanl = len(cols)
    pref = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000"
    wpath = f"{EXT2}/sv_weights_long_ext2_{ch}.tsv"
    rows_by_sv = G.load_weight_rows([wpath], set(svlist))

    want = {r["pos"] for sv in rows_by_sv for r in rows_by_sv[sv]}
    pos_map, nvar = G.index_bim(pref + ".bim", want)
    bpv = (Nfam + 3) // 4
    assert G.read_fam_iids(pref + ".fam") == fam_iids, f"{ch} fam order differs"

    cdir = f"{SCRATCH}/{ch}_b{blk:02d}"; os.makedirs(cdir, exist_ok=True)
    qc = []; prefixes = []
    fo = None; nin = 0; cidx = 0

    def close_and_import():
        nonlocal fo, nin, cidx
        if fo is None:
            return
        fo.close(); fo = None
        cf = f"{cdir}/chunk_{cidx:03d}"
        outp = f"{PGDIR}/svpgen_{ch}_b{blk:02d}_{cidx:03d}"
        subprocess.run([PLINK2, "--import-dosage", cf, "noheader", "format=1", "single-chr=1",
                        "--fam", f"{OUT}/sv.fam", "--make-pgen", "--out", outp, "--memory", "4000"],
                       check=True, capture_output=True)
        if os.path.exists(outp + "-temporary.pgen"):
            raise RuntimeError(f"{ch} b{blk} chunk {cidx} import INCOMPLETE (temp remains)")
        nv = sum(1 for l in open(outp + ".pvar") if not l.startswith("#"))
        assert nv == nin, f"{ch} b{blk} chunk {cidx}: pgen has {nv} vars, wrote {nin}"
        prefixes.append((outp, nv))
        os.remove(cf)
        cidx += 1; nin = 0

    with open(pref + ".bed", "rb") as bf:
        assert bf.read(3) == b"\x6c\x1b\x01", f"{ch}: bad .bed magic"
        for sv in svlist:
            rws = rows_by_sv.get(sv, [])
            acc = np.zeros(Nfam, dtype=np.float64)
            n_used = 0; n_sw = 0; miss_mass = 0.0; tot_mass = 0.0
            for row in rws:
                b = abs(float(row["beta_std"])); tot_mass += b
                hit = None
                for (vi, a1, a2) in pos_map.get(row["pos"], []):
                    if a1 == row["alt"] and a2 == row["ref"]:
                        hit = (vi, False); break
                    if a1 == row["ref"] and a2 == row["alt"]:
                        hit = (vi, True); break
                if hit is None:
                    miss_mass += b; continue
                vi, sw = hit
                bf.seek(3 + vi * bpv)
                a1c = G._LUT[np.frombuffer(bf.read(bpv), dtype=np.uint8)].reshape(-1)[:Nfam]
                g = (2.0 - a1c) if sw else a1c
                m = np.nanmean(g); s = np.nanstd(g)
                if not s > 0:
                    miss_mass += b; continue
                g = np.where(np.isnan(g), m, g)
                acc += float(row["beta_std"]) * (g - m) / s
                n_used += 1; n_sw += int(sw)

            d = acc[cols]
            lo, hi = float(d.min()), float(d.max())
            if n_used == 0 or not (hi > lo):
                qc.append(dict(sv_locus=sv, chrom=ch, n_total=len(rws), n_used=n_used, n_swapped=n_sw,
                               unscoreable_mass_frac=(miss_mass / tot_mass) if tot_mass else 1.0,
                               dosage_mean="", dosage_sd="", scored=0,
                               skip_reason="no usable SNP" if n_used == 0 else "degenerate (hi==lo)"))
                continue
            ds = 2.0 * (d - lo) / (hi - lo)
            if fo is None:
                fo = open(f"{cdir}/chunk_{cidx:03d}", "w")
            fo.write(f"{sv}\tA\tC\t" + "\t".join(f"{x:.5f}" for x in ds) + "\n")
            nin += 1
            qc.append(dict(sv_locus=sv, chrom=ch, n_total=len(rws), n_used=n_used, n_swapped=n_sw,
                           unscoreable_mass_frac=(miss_mass / tot_mass) if tot_mass else 0.0,
                           dosage_mean=f"{float(d.mean()):.5g}", dosage_sd=f"{float(d.std()):.5g}",
                           scored=1, skip_reason=""))
            if nin >= chunk_size:
                close_and_import()
    close_and_import()
    shutil.rmtree(cdir, ignore_errors=True)
    ns = sum(r["scored"] for r in qc)
    print(f"  [{ch} b{blk:02d}] {ns}/{len(svlist)} scored, {len(prefixes)} pgen chunks", flush=True)
    return ch, qc, prefixes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--chunk", type=int, default=300)
    ap.add_argument("--chroms", default="")
    ap.add_argument("--block", type=int, default=300)
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(PGDIR, exist_ok=True); os.makedirs(SCRATCH, exist_ok=True)
    G.guard(OUT); G.guard(PGDIR)

    fam_iids, cols, sv_fam, cov_ok = G.build_cohort()
    if not os.path.exists(f"{OUT}/sv.fam"):
        shutil.copy(f"{CARRY}/sv.fam", f"{OUT}/sv.fam")
    assert G.read_fam_iids(f"{OUT}/sv.fam") == sv_fam, "sv.fam mismatch"
    print(f"cohort N={len(cols):,}  (covariate re-derivation agrees: {cov_ok})", flush=True)

    info = preflight()
    print("chromosome preflight:", flush=True)
    for ch in sorted(info, key=chrom_order):
        i = info[ch]
        print(f"  {ch:6s} nvar={i['nvar']:>7,}  {'OK' if i['ok'] else 'EXCLUDED: ' + i['reason']}", flush=True)

    by_chrom = defaultdict(list); svchrom = {}
    meta_cv = {}
    for r in csv.DictReader(open(f"{EXT2}/sv_metadata_ext2.tsv", newline=""), delimiter="\t"):
        meta_cv[r["sv_locus"].strip()] = r.get("cv_r2", "")
    for p in sorted(glob.glob(f"{EXT2}/sv_weights_long_ext2_chr*.tsv")):
        ch = "chr" + os.path.basename(p).split("_ext2_chr")[1].replace(".tsv", "")
        seen = []
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sv = r["sv_locus"].strip()
                if sv not in svchrom:
                    svchrom[sv] = ch; seen.append(sv)
        by_chrom[ch] = seen
    total = sum(len(v) for v in by_chrom.values())
    print(f"\next2 predictors: {total:,} across {len(by_chrom)} chromosomes", flush=True)
    assert total == 11129, f"expected 11,129 ext2 predictors, got {total}"

    excluded = {ch: v for ch, v in by_chrom.items() if not info.get(ch, {}).get("ok")}
    for ch in excluded:
        print(f"  EXCLUDING {ch}: {len(by_chrom[ch])} predictors — {info[ch]['reason']}", flush=True)
    run_chroms = [ch for ch in sorted(by_chrom, key=chrom_order) if info.get(ch, {}).get("ok")]
    if a.chroms:
        want = set(a.chroms.split(",")); run_chroms = [c for c in run_chroms if c in want]
    n_run = sum(len(by_chrom[c]) for c in run_chroms)
    print(f"scoring {n_run:,} predictors on {len(run_chroms)} chromosomes "
          f"({total - n_run:,} excluded)\n", flush=True)

    tier, old_tier = build_tier_map(svchrom)
    from collections import Counter
    print("tier composition (all 11,129):", dict(Counter(tier.values())), flush=True)

    BLOCK = a.block
    tasks = []
    for ch in run_chroms:
        svl = by_chrom[ch]
        for b, i in enumerate(range(0, len(svl), BLOCK)):
            tasks.append((ch, b, svl[i:i + BLOCK], fam_iids, cols, a.chunk))
    tasks.sort(key=lambda t: -len(t[2]))
    print(f"{len(tasks)} blocks of <= {BLOCK} SVs across {len(run_chroms)} chromosomes "
          f"on {a.workers} workers", flush=True)
    allqc = []; allpre = []
    with Pool(a.workers) as pool:
        for ch, qc, pre in pool.imap_unordered(score_chrom, tasks):
            allqc.extend(qc); allpre.extend([(ch, p, n) for p, n in pre])

    scored = [r for r in allqc if r["scored"]]
    print(f"\nscored {len(scored):,} of {n_run:,} attempted "
          f"({n_run - len(scored):,} skipped as unscoreable/degenerate)", flush=True)

    ks = ["sv_locus", "tier", "chrom", "n_total", "n_used", "n_swapped", "unscoreable_mass_frac",
          "dosage_mean", "dosage_sd", "cv_r2", "scored", "skip_reason"]
    with open(G.guard(f"{OUT}/dosage_qc_ext2.tsv"), "w") as f:
        f.write("\t".join(ks) + "\n")
        for r in sorted(allqc, key=lambda r: r["sv_locus"]):
            r = dict(r, tier=tier[r["sv_locus"]], cv_r2=meta_cv.get(r["sv_locus"], ""))
            r["unscoreable_mass_frac"] = f"{r['unscoreable_mass_frac']:.5g}"
            f.write("\t".join(str(r[k]) for k in ks) + "\n")

    sset = {r["sv_locus"] for r in scored}
    with open(G.guard(f"{OUT}/sv_tier.tsv"), "w") as f:
        f.write("sv_locus\ttier\n")
        for sv in sorted(sset):
            f.write(f"{sv}\t{tier[sv]}\n")
    with open(G.guard(f"{OUT}/sv_cv_r2_all.tsv"), "w") as f:
        f.write("sv_locus\tcv_r2\n")
        for sv in sorted(sset):
            f.write(f"{sv}\t{meta_cv.get(sv,'')}\n")

    allpre.sort(key=lambda t: (chrom_order(t[0]), os.path.basename(t[1])))
    with open(G.guard(f"{OUT}/pgen_manifest.tsv"), "w") as f:
        f.write("order\tchrom\tprefix\tn_variants\n")
        for i, (ch, p, n) in enumerate(allpre):
            f.write(f"{i}\t{ch}\t{os.path.basename(p)}\t{n}\n")
    tot_pgen = sum(n for _, _, n in allpre)
    assert tot_pgen == len(scored), f"pgen variants {tot_pgen} != scored {len(scored)}"

    summ = dict(n_ext2=total, n_attempted=n_run, n_scored=len(scored),
                n_excluded_chrom=sum(len(by_chrom[c]) for c in excluded),
                excluded_chroms={c: info[c]["reason"] for c in excluded},
                n_pgen_chunks=len(allpre), tiers=dict(Counter(tier[s] for s in sset)),
                secs=time.time() - t0)
    with open(G.guard(f"{OUT}/step50_summary.json"), "w") as f:
        json.dump(summ, f, indent=1)
    print(json.dumps(summ, indent=1), flush=True)
    print(f"\nSTEP50_DONE in {(time.time()-t0)/60:.1f} min -> {PGDIR}", flush=True)


if __name__ == "__main__":
    main()
