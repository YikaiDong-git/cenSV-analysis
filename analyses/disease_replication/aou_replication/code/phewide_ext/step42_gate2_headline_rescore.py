#!/usr/bin/env python3
"""Step 42 — GATE 2: does the 11,129 ext2 SUPERSET weighting actually move our published headline?"""
import os, sys, csv, glob, json, time, hashlib
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8"); os.environ.setdefault("OMP_NUM_THREADS", "8")
import numpy as np
import pgenlib
from collections import defaultdict
from scipy.stats import chi2 as chi2dist

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
PH = f"{BASE}/phewide_all_handoff"
EXT2 = f"{PH}/full_11129_ext2"
BED = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
COVSRC = f"{HOME}/workspace/working/scripts/covariates/output/covariates_v3_all_ancestries.tsv"
RV = f"{BASE}/aou_replication/results_vds"
PGDIR = f"{RV}/step28_phewide/pgen"
CARRY = f"{RV}/phewide_carryover"
HITS = f"{RV}/step30_phewide/disease_matrices_firth/disease_hits_final.tsv"
CELLS = f"{RV}/step30_phewide/disease_matrices_firth/disease_firth_cells.tsv"
OUT = f"{RV}/step42_gate2"
GATES = f"{OUT}/gates"
RESULTS = f"{OUT}/results"
THRU = f"{OUT}/throughput"
ROWLVL = f"{OUT}/rowlevel"

PGEN_FP = "2da66b5a70a7789686a62c5ce76a7d60"
R_MIN = 0.999
FROZEN = [f"{RV}/step30_phewide", f"{RV}/step28_phewide", CARRY]

_CODE2A1 = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)
_LUT = np.empty((256, 4), dtype=np.float64)
for _b in range(256):
    for _k in range(4):
        _LUT[_b, _k] = _CODE2A1[(_b >> (2 * _k)) & 3]

BASECOV = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
           "sex", *[f"gPC{i}" for i in range(1, 17)]]
FIRTH_MAXIT = 25; FIRTH_PEN = 0.5
IDENT_TOL = 1.0; VREL = 1e-6


def guard(path):
    if os.path.islink(path):
        raise SystemExit(f"ABORT V5: output target is a symlink: {path}")
    d = os.path.dirname(path)
    if os.path.islink(d):
        raise SystemExit(f"ABORT V5: output directory is a symlink: {d}")
    return path


def tree_manifest(root):
    man = {}
    for dp, _dn, fn in os.walk(root):
        for f in fn:
            p = os.path.join(dp, f)
            try:
                st = os.stat(p)
            except OSError:
                continue
            man[os.path.relpath(p, root)] = [st.st_size, st.st_mtime_ns]
    return man


def manifest_sig(man):
    return hashlib.md5(json.dumps(man, sort_keys=True).encode()).hexdigest()


def _pgen_fingerprint():
    fs = sorted(glob.glob(f"{PGDIR}/svpgen300_*.pgen"),
                key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
    acc = hashlib.md5()
    for f in fs:
        h = hashlib.md5()
        with open(f, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                h.update(blk)
        acc.update(h.hexdigest().encode())
    return acc.hexdigest(), len(fs)


def read_fam_iids(path):
    out = []
    with open(path) as f:
        for line in f:
            out.append(line.split()[1])
    return out


def chunk_prefixes():
    fs = sorted(glob.glob(f"{PGDIR}/svpgen300_*.pgen"),
                key=lambda p: int(p.split("svpgen300_")[1].split(".pgen")[0]))
    return [f[:-len(".pgen")] for f in fs]


def psam_iids(pfx):
    out = []
    for line in open(f"{pfx}.psam"):
        if line.startswith("#"):
            continue
        out.append(line.split("\t")[1] if "\t" in line else line.split()[1])
    return out


def build_cohort():
    fam_iids = None
    for ch in sorted(os.listdir(BED)):
        p = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000.fam"
        if os.path.exists(p):
            fam_iids = read_fam_iids(p); break
    if fam_iids is None:
        raise SystemExit("ABORT: no VDS .fam found")

    sv_fam = read_fam_iids(f"{CARRY}/sv.fam")
    posn = {iid: i for i, iid in enumerate(fam_iids)}
    missing = [i for i in sv_fam if i not in posn]
    if missing:
        raise SystemExit(f"ABORT: {len(missing)} sv.fam IIDs absent from the VDS fam")
    cols = np.array([posn[i] for i in sv_fam], dtype=np.int64)
    if not np.all(np.diff(cols) > 0):
        raise SystemExit("ABORT: sv.fam is not an order-preserving subset of the VDS fam — "
                         "step28's `cols` filter cannot be reproduced")

    covset = set()
    with open(COVSRC) as f:
        h = f.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if all(p[idx[c]] not in ("", "NA") for c in BASECOV if c in idx):
                covset.add(p[idx["IID"]])
    cols_b = np.array([i for i, iid in enumerate(fam_iids) if iid in covset], dtype=np.int64)
    agree = (len(cols_b) == len(cols)) and bool(np.array_equal(cols_b, cols))
    print(f"  cohort: VDS fam N={len(fam_iids):,}  analysis N={len(cols):,}  "
          f"covariate-rederivation {'AGREES' if agree else '*** DISAGREES ***'}", flush=True)
    if not agree:
        print(f"  NOTE: re-derived N={len(cols_b):,} != frozen N={len(cols):,}; using the FROZEN sv.fam "
              f"(authoritative for what was published).", flush=True)
    return fam_iids, cols, sv_fam, agree


def load_weight_rows(paths, want_svs):
    out = defaultdict(list)
    for p in paths:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sv = (r.get("sv_locus") or "").strip()
                if sv in want_svs:
                    out[sv].append({k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()})
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


def score_svs(rows_by_sv, fam_iids, cols, label):
    Nfam = len(fam_iids); Nanl = len(cols)
    by_chrom = defaultdict(dict)
    for sv, rws in rows_by_sv.items():
        for r in rws:
            ch = r["chrom"] if str(r["chrom"]).startswith("chr") else f"chr{r['chrom']}"
            by_chrom[ch].setdefault(sv, []).append(r)
    res = {}
    for ch in sorted(by_chrom):
        pref = f"{BED}/{ch}/{ch}.pericentromeric_snps.min_ac1000"
        if not os.path.exists(pref + ".bed"):
            print(f"  [{label}] {ch}: NO BED -> skipped", flush=True); continue
        assert read_fam_iids(pref + ".fam") == fam_iids, f"{ch} fam order differs"
        svrows = by_chrom[ch]
        want = {r["pos"] for sv in svrows for r in svrows[sv]}
        pos_map, nvar = index_bim(pref + ".bim", want)
        bpv = (Nfam + 3) // 4
        with open(pref + ".bed", "rb") as bf:
            assert bf.read(3) == b"\x6c\x1b\x01", f"{ch}: bad .bed magic"
            for sv, rws in svrows.items():
                t0 = time.time()
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
                res[sv] = (ds, len(rws), n_used, n_sw, time.time() - t0)
        print(f"  [{label}] {ch}: scored {len(svrows)} SVs; bim {nvar} vars", flush=True)
    return res


def read_archived_dosages(loci):
    want = set(loci); out = {}
    for pfx in chunk_prefixes():
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        local = {sid: i for i, sid in enumerate(ids) if sid in want}
        if not local:
            continue
        n = len(psam_iids(pfx)); buf = np.empty(n, np.float64)
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for sid, v in local.items():
            r.read_dosages(v, buf); out[sid] = buf.copy()
        r.close()
    return out


def load_X(iidx):
    N = len(iidx)
    with open(f"{CARRY}/inputs/covar.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ci = [j for j, c in enumerate(hdr) if c not in ("FID", "IID")]
        names = ["intercept"] + [hdr[j] for j in ci]
        iidcol = hdr.index("IID")
        X = np.zeros((N, len(ci) + 1), np.float64); X[:, 0] = 1.0
        seen = 0
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[iidcol])
            if kk is None:
                continue
            for jj, j in enumerate(ci):
                X[kk, jj + 1] = float(p[j])
            seen += 1
    assert seen == N, f"covar.txt covered {seen}/{N} cohort samples"
    return X, names


def load_pheno_cols(phecodes, iidx):
    N = len(iidx)
    hdr = open(f"{CARRY}/inputs/disease.pheno").readline().rstrip("\n").split("\t")
    colof = {ph: hdr.index(ph) for ph in phecodes if ph in hdr}
    miss = [p for p in phecodes if p not in colof]
    if miss:
        raise SystemExit(f"ABORT: phecodes absent from disease.pheno: {miss}")
    Y = {ph: np.full(N, -1, np.int8) for ph in colof}
    with open(f"{CARRY}/inputs/disease.pheno") as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t"); kk = iidx.get(p[1])
            if kk is None:
                continue
            for ph, c in colof.items():
                if p[c] == "0":
                    Y[ph][kk] = 0
                elif p[c] == "1":
                    Y[ph][kk] = 1
    return Y


def fit_null(Xk, y, present):
    kk = Xk.shape[1]
    beta = np.zeros(kk)
    ybar = float(y[present].mean())
    ybar = min(max(ybar, 1e-6), 1 - 1e-6)
    beta[0] = np.log(ybar / (1 - ybar))
    converged = False
    for _ in range(25):
        eta = np.clip(Xk @ beta, -30, 30)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1 - p); w[~present] = 0.0
        r = y - p; r[~present] = 0.0
        WX = Xk * w[:, None]
        XtWX = Xk.T @ WX
        grad = Xk.T @ r
        try:
            delta = np.linalg.solve(XtWX, grad)
        except np.linalg.LinAlgError:
            delta = np.linalg.lstsq(XtWX, grad, rcond=None)[0]
        beta += delta
        if np.max(np.abs(delta)) < 1e-6:
            converged = True
            break
    eta = np.clip(Xk @ beta, -30, 30); p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1 - p); w[~present] = 0.0
    r = y - p; r[~present] = 0.0
    return p, w, r, converged


def null_cache(X, Y):
    cache = {}
    for ph, yb in Y.items():
        t0 = time.time()
        present = yb >= 0
        nca = int((yb == 1).sum()); npr = int(present.sum())
        y = np.where(present, yb, 0).astype(np.float64)
        _, w, r, conv = fit_null(X, y, present)
        degen = abs(float(r.sum())) > IDENT_TOL
        WX = X * w[:, None]
        XtWXi = np.linalg.pinv(X.T @ WX)
        cache[ph] = dict(w=w, r=r, WX=WX, XtWXi=XtWXi, present=present, ncase=nca,
                         npres=npr, conv=conv, degen=degen)
        print(f"    null {ph}: cases={nca:,} n={npr:,} conv={conv} degen={degen} "
              f"({time.time()-t0:.0f}s)", flush=True)
    return cache


def score_test_one(g, nul):
    if nul["degen"]:
        return None
    U = float(g @ nul["r"])
    A = g @ nul["WX"]
    Vraw = float((g * g) @ nul["w"])
    V = Vraw - float(A @ nul["XtWXi"] @ A)
    if not (np.isfinite(V) and V > VREL * Vraw and Vraw > 0):
        return None
    stat = U * U / V
    return float(U / V), float(chi2dist.sf(stat, 1)), float(stat)


_XCOV = None; _DOS = None; _YV = None


def firth(Xcov, g, y, present):
    from firthmodels import FirthLogisticRegression
    Xp = Xcov[present]; gp = g[present]; yp = y[present].astype(int)
    keep = np.array([np.ptp(Xp[:, j]) > 0 for j in range(Xp.shape[1])])
    Xf = np.column_stack([Xp[:, keep], gp]); vi = int(keep.sum())
    m = FirthLogisticRegression(fit_intercept=True, max_iter=FIRTH_MAXIT, penalty_weight=FIRTH_PEN)
    try:
        m.fit(Xf, yp)
    except Exception:
        return np.nan, np.nan, False
    beta = float(m.coef_[vi])
    try:
        m.lrt(vi); p = float(m.lrt_pvalues_[vi])
    except Exception:
        p = float(m.pvalues_[vi])
    return beta, p, bool(getattr(m, "converged_", True))


def _firth_task(t):
    idx, which, sv, ph = t
    g = _DOS[which][sv]; y = _YV[ph]
    b, p, c = firth(_XCOV, g, y, y >= 0)
    return idx, (b, p, c)


def pearson(a, b):
    a = np.asarray(a, np.float64); b = np.asarray(b, np.float64)
    a = a - a.mean(); b = b - b.mean()
    da = np.sqrt((a * a).sum()); db = np.sqrt((b * b).sum())
    if da == 0 or db == 0:
        return np.nan
    return float((a * b).sum() / (da * db))


def cosine_w(rows_a, rows_b):
    ka = {(r["pos"], r["ref"], r["alt"]): float(r["beta_std"]) for r in rows_a}
    kb = {(r["pos"], r["ref"], r["alt"]): float(r["beta_std"]) for r in rows_b}
    keys = sorted(set(ka) | set(kb))
    va = np.array([ka.get(k, 0.0) for k in keys]); vb = np.array([kb.get(k, 0.0) for k in keys])
    na = np.linalg.norm(va); nb = np.linalg.norm(vb)
    return float(va @ vb / (na * nb)) if na and nb else np.nan


def main():
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    assert stage in ("verify", "score", "refit", "all"), f"unknown stage {stage}"
    t00 = time.time()
    for d in (OUT, ROWLVL, GATES, RESULTS, THRU):
        os.makedirs(d, exist_ok=True)
    guard(OUT); guard(ROWLVL)

    print("== V4-pre: frozen-tree manifests ==", flush=True)
    pre = {}
    for root in FROZEN:
        m = tree_manifest(root); pre[root] = m
        print(f"  {os.path.basename(root):22s} {len(m):5d} files  sig={manifest_sig(m)}", flush=True)
    with open(guard(f"{GATES}/_frozen_manifest_pre.json"), "w") as f:
        json.dump({k: manifest_sig(v) for k, v in pre.items()}, f, indent=1)

    print("== V2: archived pgen provenance ==", flush=True)
    fp, nchunk = _pgen_fingerprint()
    print(f"  fingerprint {fp}  ({nchunk} chunks)", flush=True)
    if fp != PGEN_FP:
        raise SystemExit(f"ABORT V2: pgen fingerprint {fp} != expected {PGEN_FP}")
    pfxs = chunk_prefixes()
    p0 = psam_iids(pfxs[0])
    sv_fam = read_fam_iids(f"{CARRY}/sv.fam")
    if p0 != sv_fam:
        raise SystemExit(f"ABORT V2: psam order != sv.fam order ({len(p0)} vs {len(sv_fam)})")
    sig0 = hashlib.md5("\n".join(p0).encode()).hexdigest()
    for pfx in pfxs[1:]:
        if hashlib.md5("\n".join(psam_iids(pfx)).encode()).hexdigest() != sig0:
            raise SystemExit(f"ABORT V2: psam order mismatch at {pfx}")
    print(f"  psam order == sv.fam for all {nchunk} chunks (N={len(p0):,})  OK", flush=True)

    hits = [r for r in csv.DictReader(open(HITS), delimiter="\t")]
    sig = [r for r in hits if float(r["q_bh_tier"]) < 0.05]
    SVS = sorted({r["sv_locus"] for r in sig}); PHS = sorted({r["phecode"] for r in sig})
    print(f"== targets: {len(sig)} q<0.05 cells = {len(SVS)} SVs x {len(PHS)} phecodes ==", flush=True)

    fam_iids, cols, sv_fam2, cov_agree = build_cohort()
    assert sv_fam2 == sv_fam

    old_paths = sorted(glob.glob(f"{PH}/sv_weights_long_chr*.tsv")) + \
        sorted(glob.glob(f"{PH}/discovery_null_ext/sv_weights_long_ext_chr*.tsv"))
    new_paths = sorted(glob.glob(f"{EXT2}/sv_weights_long_ext2_chr*.tsv"))
    W_old = load_weight_rows(old_paths, set(SVS))
    W_new = load_weight_rows(new_paths, set(SVS))
    for sv in SVS:
        if sv not in W_old:
            raise SystemExit(f"ABORT: {sv} absent from OLD weights")
        if sv not in W_new:
            raise SystemExit(f"ABORT: {sv} absent from ext2 weights")
    print(f"  weights: old {sum(len(v) for v in W_old.values())} rows / "
          f"new {sum(len(v) for v in W_new.values())} rows over {len(SVS)} SVs", flush=True)

    print("== V1: scorer fidelity (re-score the 8 SVs from OLD weights) ==", flush=True)
    arch = read_archived_dosages(SVS)
    missing = [s for s in SVS if s not in arch]
    if missing:
        raise SystemExit(f"ABORT V1: SVs not in archived pgen: {missing}")
    myold = score_svs(W_old, fam_iids, cols, "old")
    v1 = {}
    for sv in SVS:
        r = pearson(arch[sv], myold[sv][0])
        v1[sv] = r
        flag = "OK " if (np.isfinite(r) and r > R_MIN) else "FAIL"
        print(f"  {flag} r={r:.8f}  {sv}", flush=True)
    bad = [s for s, r in v1.items() if not (np.isfinite(r) and r > R_MIN)]
    if bad:
        raise SystemExit(f"ABORT V1: {len(bad)} SVs below r>{R_MIN}: {bad}\n"
                         f"The reimplementation is not faithful to step28 — every new-weight number "
                         f"downstream would be meaningless. Fix before proceeding.")
    print(f"  V1 PASSED: all {len(SVS)} SVs r>{R_MIN} (min {min(v1.values()):.8f})", flush=True)
    with open(guard(f"{GATES}/gate2_v1_scorer_fidelity.tsv"), "w") as f:
        f.write("sv_locus\tn_snps_old\tn_used_old\tn_swapped_old\tr_archived_vs_rescored\n")
        for sv in SVS:
            _, nt, nu, nsw, _ = myold[sv]
            f.write(f"{sv}\t{nt}\t{nu}\t{nsw}\t{v1[sv]:.10f}\n")
    if stage == "verify":
        _v4post(pre, t00); return

    print("== new-weight scoring (ext2) ==", flush=True)
    t_new0 = time.time()
    mynew = score_svs(W_new, fam_iids, cols, "new")
    t_new = time.time() - t_new0
    for sv in SVS:
        if sv not in mynew:
            raise SystemExit(f"ABORT: {sv} could not be scored from ext2 weights")

    rows = []
    for sv in SVS:
        do, nt_o, nu_o, nsw_o, _ = myold[sv]
        dn, nt_n, nu_n, nsw_n, sec = mynew[sv]
        a = arch[sv]
        wa = np.array([abs(float(r["beta_std"])) for r in W_old[sv]]).sum()
        keys_old = {(r["pos"], r["ref"], r["alt"]) for r in W_old[sv]}
        added = sum(abs(float(r["beta_std"])) for r in W_new[sv]
                    if (r["pos"], r["ref"], r["alt"]) not in keys_old)
        tot_new = sum(abs(float(r["beta_std"])) for r in W_new[sv])
        rows.append(dict(
            sv_locus=sv, n_snps_old=nt_o, n_snps_new=nt_n, n_used_old=nu_o, n_used_new=nu_n,
            n_swapped_old=nsw_o, n_swapped_new=nsw_n,
            added_weight_mass_frac=added / tot_new if tot_new else np.nan,
            cosine_weights=cosine_w(W_old[sv], W_new[sv]),
            r_archived_vs_rescored_old=v1[sv],
            r_archived_vs_new=pearson(a, dn),
            r_rescored_old_vs_new=pearson(do, dn),
            sd_archived=float(np.std(a)), sd_new=float(np.std(dn)),
            secs_new_scoring=sec))
    cols_out = list(rows[0].keys())
    with open(guard(f"{RESULTS}/gate2_dosage_correlation.tsv"), "w") as f:
        f.write("\t".join(cols_out) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols_out) + "\n")
    print("  --- empirical dosage correlation (archived published vs ext2) ---", flush=True)
    for r in sorted(rows, key=lambda r: r["r_archived_vs_new"]):
        print(f"   r={r['r_archived_vs_new']:+.6f}  cos={r['cosine_weights']:.4f}  "
              f"snps {r['n_snps_old']}->{r['n_snps_new']}  {r['sv_locus']}", flush=True)
    np.savez_compressed(guard(f"{ROWLVL}/dosages.npz"),
                        **{f"arch::{s}": arch[s].astype(np.float32) for s in SVS},
                        **{f"myold::{s}": myold[s][0].astype(np.float32) for s in SVS},
                        **{f"new::{s}": mynew[s][0].astype(np.float32) for s in SVS})
    with open(guard(f"{ROWLVL}/README.md"), "w") as f:
        f.write("# ROW-LEVEL — controlled tier, bucket only\n\n"
                "`dosages.npz` holds per-person SV dosage vectors (N=315,204) for the 8 headline SVs "
                "under three weightings. This is individual-level derived genotype data: it must NEVER "
                "be committed to git, emailed, or copied outside the controlled workspace. Everything "
                "else written by step42 is aggregate.\n")
    print(f"  new-weight scoring wall-clock: {t_new:.1f}s for {len(SVS)} SVs "
          f"({t_new/len(SVS):.2f}s/SV) -> 11,129 SVs ~= {t_new/len(SVS)*11129/3600:.2f} h single-process",
          flush=True)
    if stage == "score":
        _v4post(pre, t00); return

    print("== loading covariates + phenotypes ==", flush=True)
    iidx = {i: k for k, i in enumerate(sv_fam)}
    global _XCOV, _DOS, _YV
    X, xnames = load_X(iidx)
    print(f"  X {X.shape}", flush=True)
    t0 = time.time(); Y = load_pheno_cols(PHS, iidx)
    print(f"  Y {len(Y)} phecodes ({time.time()-t0:.0f}s)", flush=True)

    print("== covariate-null fits (one per phecode) ==", flush=True)
    nulls = null_cache(X, Y)

    published = {(r["sv_locus"], r["phecode"]): r for r in csv.DictReader(open(CELLS), delimiter="\t")}
    final = {(r["sv_locus"], r["phecode"]): r for r in sig}

    print("== score tests: archived (V3) vs new dosages ==", flush=True)
    st = {}
    for sv in SVS:
        for ph in PHS:
            if (sv, ph) not in final:
                continue
            st[(sv, ph, "old")] = score_test_one(arch[sv], nulls[ph])
            st[(sv, ph, "new")] = score_test_one(mynew[sv][0], nulls[ph])

    v3a = []
    for (sv, ph), pr in published.items():
        if (sv, ph) not in final:
            continue
        got = st.get((sv, ph, "old"))
        try:
            exp = float(pr["p_score"])
        except (TypeError, ValueError):
            continue
        if got is None:
            v3a.append((sv, ph, exp, np.nan, np.inf)); continue
        _, p, _ = got
        rel = abs(p - exp) / exp if exp > 0 else (0.0 if p == 0 else np.inf)
        v3a.append((sv, ph, exp, p, rel))
    worst = max((x[4] for x in v3a), default=0.0)
    print(f"  V3a score-test fidelity: {len(v3a)} cells, worst relative deviation {worst:.3e}", flush=True)
    for sv, ph, exp, p, rel in sorted(v3a, key=lambda x: -x[4])[:5]:
        print(f"    rel={rel:.2e}  published={exp:.6g} recomputed={p:.6g}  {sv} {ph}", flush=True)

    print("== exact Firth refits (28 cells x 2 dosage sets) ==", flush=True)
    from multiprocessing import Pool
    _XCOV = X[:, 1:]
    _DOS = {"old": {s: arch[s] for s in SVS}, "new": {s: mynew[s][0] for s in SVS}}
    _YV = Y
    tasks = []
    for (sv, ph) in sorted(final):
        for which in ("old", "new"):
            tasks.append((len(tasks), which, sv, ph))
    fres = [None] * len(tasks)
    t0 = time.time()
    with Pool(6) as pool:
        done = 0
        for idx, r in pool.imap_unordered(_firth_task, tasks, chunksize=1):
            fres[idx] = r; done += 1
            if done % 8 == 0:
                print(f"    firth {done}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)
    fm = {(t[2], t[3], t[1]): fres[t[0]] for t in tasks}
    print(f"  firth done ({time.time()-t0:.0f}s)", flush=True)

    v3b = []
    for (sv, ph) in sorted(final):
        fr = final[(sv, ph)]
        if fr["p_source"] != "exact":
            continue
        exp = float(fr["p_final"]); expb = float(fr["beta_final_logOR"])
        b, p, c = fm[(sv, ph, "old")]
        rel = abs(p - exp) / exp if exp > 0 else (0.0 if p == 0 else np.inf)
        relb = abs(b - expb) / abs(expb) if expb else np.nan
        v3b.append((sv, ph, exp, p, rel, expb, b, relb))
    worstb = max((x[4] for x in v3b), default=0.0)
    worstbb = max((x[7] for x in v3b), default=0.0)
    n_approx = sum(1 for k in final if final[k]["p_source"] != "exact")
    print(f"  V3b exact-Firth fidelity vs p_final: {len(v3b)} exact-source cells, "
          f"worst rel dev p={worstb:.3e} beta={worstbb:.3e}   ({n_approx} approx-source cells excluded "
          f"— their p_final is step30h2's APPROXIMATE Firth, a different estimator)", flush=True)
    for sv, ph, exp, p, rel, eb, b, rb in sorted(v3b, key=lambda x: -x[4])[:5]:
        print(f"    rel={rel:.2e}  published={exp:.6g} recomputed={p:.6g}  {sv} {ph}", flush=True)

    sig_p = [float(r["p_final"]) for r in hits if r["tier"] == "signal" and float(r["q_bh_tier"]) < 0.05]
    nonsig_p = [float(r["p_final"]) for r in hits if r["tier"] == "signal" and float(r["q_bh_tier"]) >= 0.05]
    bh_hi = max(sig_p); bh_lo = min(nonsig_p) if nonsig_p else np.inf
    print(f"  published signal-tier BH q<0.05 band: largest passing p={bh_hi:.4g}, "
          f"smallest failing p={bh_lo:.4g}", flush=True)

    outrows = []
    for (sv, ph) in sorted(final, key=lambda k: float(final[k]["p_final"])):
        fr = final[(sv, ph)]
        so = st[(sv, ph, "old")]; sn = st[(sv, ph, "new")]
        bo, po, co = fm[(sv, ph, "old")]; bn, pn, cn = fm[(sv, ph, "new")]
        src = fr["p_source"]
        p_pub = float(fr["p_final"])
        p_old_eq = po
        p_new_eq = pn
        outrows.append(dict(
            sv_locus=sv, phecode=ph, label=fr["label"], cases=fr["cases"], n=fr["n"],
            p_published=p_pub, q_published=float(fr["q_bh_tier"]), p_source=src,
            beta_published=float(fr["beta_final_logOR"]),
            p_score_old=so[1] if so else np.nan, p_score_new=sn[1] if sn else np.nan,
            beta_score_old=so[0] if so else np.nan, beta_score_new=sn[0] if sn else np.nan,
            p_firth_old=po, p_firth_new=pn, beta_firth_old=bo, beta_firth_new=bn,
            firth_conv_old=int(bool(co)), firth_conv_new=int(bool(cn)),
            p_old_likeforlike=p_old_eq, p_new_likeforlike=p_new_eq,
            log10_p_shift=(np.log10(p_new_eq) - np.log10(p_old_eq))
            if (np.isfinite(p_new_eq) and np.isfinite(p_old_eq) and p_new_eq > 0 and p_old_eq > 0) else np.nan,
            beta_sign_same=int(np.isfinite(bn) and np.isfinite(bo) and np.sign(bn) == np.sign(bo)),
            beta_ratio_new_over_old=(bn / bo) if (np.isfinite(bn) and np.isfinite(bo) and bo != 0) else np.nan,
            cutoff_band=("KEEP" if (np.isfinite(p_new_eq) and p_new_eq <= bh_hi * (1 + 1e-3))
                         else "BORDER" if (np.isfinite(p_new_eq) and p_new_eq <= bh_lo) else "LOST"),
            r_dosage=pearson(arch[sv], mynew[sv][0]),
        ))
    ck = list(outrows[0].keys())
    with open(guard(f"{RESULTS}/gate2_cell_refit.tsv"), "w") as f:
        f.write("\t".join(ck) + "\n")
        for r in outrows:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                              for c in ck) + "\n")
    with open(guard(f"{GATES}/gate2_v3_refit_fidelity.tsv"), "w") as f:
        f.write("check\tsv_locus\tphecode\tpublished\trecomputed\trel_dev\n")
        for sv, ph, e, p, rl in v3a:
            f.write(f"V3a_score\t{sv}\t{ph}\t{e:.10g}\t{p:.10g}\t{rl:.3e}\n")
        for sv, ph, e, p, rl, eb, b, rb in v3b:
            f.write(f"V3b_firth_p\t{sv}\t{ph}\t{e:.10g}\t{p:.10g}\t{rl:.3e}\n")
            f.write(f"V3b_firth_beta\t{sv}\t{ph}\t{eb:.10g}\t{b:.10g}\t{rb:.3e}\n")

    from collections import Counter
    bands = Counter(r["cutoff_band"] for r in outrows)
    n_clear = bands["KEEP"]
    print("\n  --- 28 cells: exact Firth on ARCHIVED -> ext2 dosages ---", flush=True)
    for r in outrows:
        print(f"   {r['p_old_likeforlike']:>10.3g} -> {r['p_new_likeforlike']:>10.3g}  "
              f"dlog10={r['log10_p_shift']:+6.2f}  r={r['r_dosage']:.4f}  "
              f"{r['cutoff_band']:6s} {r['sv_locus'][:34]:34s} {r['phecode']}", flush=True)
    print(f"\n  vs the published BH cutoff (KEEP p<={bh_hi:.3g} / BORDER <={bh_lo:.3g} / LOST above): "
          f"KEEP={bands['KEEP']}  BORDER={bands['BORDER']}  LOST={bands['LOST']}  of {len(outrows)}",
          flush=True)
    per_sv = defaultdict(Counter)
    for r in outrows:
        per_sv[r["sv_locus"]][r["cutoff_band"]] += 1
    print("  per SV:", flush=True)
    for sv in sorted(per_sv, key=lambda s: -per_sv[s]["KEEP"]):
        rr = next(r["r_dosage"] for r in outrows if r["sv_locus"] == sv)
        c = per_sv[sv]
        print(f"    r={rr:.4f}  KEEP={c['KEEP']} BORDER={c['BORDER']} LOST={c['LOST']}   {sv}", flush=True)

    with open(guard(f"{RESULTS}/gate2_summary.json"), "w") as f:
        json.dump(dict(
            n_cells=len(outrows), n_svs=len(SVS), n_phecodes=len(PHS),
            v1_min_r=min(v1.values()), v3a_worst_rel=worst,
            v3b_worst_rel_p=worstb, v3b_worst_rel_beta=worstbb, n_approx_source=n_approx,
            cov_rederivation_agrees=cov_agree,
            bh_cutoff_hi=bh_hi, bh_cutoff_lo=bh_lo, n_cells_clearing_cutoff=n_clear, bands=dict(bands),
            secs_new_scoring_total=t_new, secs_per_sv=t_new / len(SVS),
            pgen_fingerprint=fp), f, indent=1)
    _v4post(pre, t00)


def _v4post(pre, t00):
    print("== V4-post: frozen-tree integrity ==", flush=True)
    ok = True
    report = {}
    for root in FROZEN:
        post = tree_manifest(root); prem = pre[root]
        added = sorted(set(post) - set(prem)); removed = sorted(set(prem) - set(post))
        changed = sorted(k for k in set(post) & set(prem) if post[k] != prem[k])
        report[os.path.basename(root)] = dict(n=len(post), added=added[:20],
                                              removed=removed[:20], changed=changed[:20])
        status = "OK" if not (added or removed or changed) else "*** MODIFIED ***"
        print(f"  {os.path.basename(root):22s} {len(post):5d} files  "
              f"+{len(added)} -{len(removed)} ~{len(changed)}  {status}", flush=True)
        if added or removed or changed:
            ok = False
            for k in (added + removed + changed)[:10]:
                print(f"      {k}", flush=True)
    with open(guard(f"{GATES}/_frozen_manifest_post.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(f"\n{'GATE2 V4 PASSED — 0 files changed in the frozen trees' if ok else 'GATE2 V4 FAILED'}"
          f"   total {time.time()-t00:.0f}s", flush=True)
    if not ok:
        raise SystemExit("ABORT V4: a frozen tree was modified")


if __name__ == "__main__":
    main()
