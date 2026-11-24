#!/usr/bin/env python3
import os, sys, gc, gzip, csv, numpy as np
import pgenlib
from scipy.stats import rankdata, norm, t as tdist

BASE = os.path.expanduser("~/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication")
H = f"{BASE}/phewide_222_handoff"
R6 = f"{BASE}/aou_replication/results/step6_222"
PGEN = os.path.expanduser("~/workspace/vwb-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen")
COVDIR = os.path.expanduser("~/workspace/working/scripts/covariates/output")
OUT = f"{BASE}/aou_replication/results/step10_cond"; os.makedirs(OUT, exist_ok=True)

WINDOWS_MB = [1, 2, 5]
PCS = [10, 20, 40, 100, 200]
MMAX = {"eur": 500, "all": 250}
PCKEEP = {"eur": 200, "all": 40}
CENTRAL_W, CENTRAL_K, FLANK_KB = 2, 20, 10

SVS = [
    {"sv": "chr3_sv20_DEL_non_best6_alt_edge1000", "chrom": "chr3", "cv_r2": 0.201,
     "traits": [{"loinc": "29463-7", "name": "weight", "disc": "-"}]},
    {"sv": "chr4_sv349_DEL_non_best8_fksize200k", "chrom": "chr4", "cv_r2": 0.116,
     "traits": [{"loinc": "8302-2", "name": "height", "disc": "-"}]},
    {"sv": "chr11_sv263_DEL_non_best8_edge500", "chrom": "chr11", "cv_r2": 0.425,
     "traits": [{"loinc": "39156-5", "name": "BMI", "disc": "-"}]},
    {"sv": "chr4_sv398_INS_non_best8_fksize100k", "chrom": "chr4", "cv_r2": 0.120,
     "traits": [{"loinc": "32623-1", "name": "MPV_30100", "disc": "-"},
                {"loinc": "776-5", "name": "PDW_30110", "disc": "-"}]},
]
REF_EUR = {
    ("chr3_sv20_DEL_non_best6_alt_edge1000", "weight"): -0.087168,
    ("chr4_sv349_DEL_non_best8_fksize200k", "height"): -0.047159,
    ("chr11_sv263_DEL_non_best8_edge500", "BMI"): -0.059771,
    ("chr4_sv398_INS_non_best8_fksize100k", "MPV_30100"): -0.18691,
    ("chr4_sv398_INS_non_best8_fksize100k", "PDW_30110"): -0.26122,
}
HIT_SVS = [s["sv"] for s in SVS]
ALL_LOINCS = sorted({t["loinc"] for s in SVS for t in s["traits"]})


def irnt(x):
    m = np.isfinite(x); out = np.full(x.shape, np.nan); r = rankdata(x[m])
    out[m] = norm.ppf((r - 0.5) / m.sum()); return out


def ols_dose(X, y):
    XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ (X.T @ y); resid = y - X @ beta
    n, p = X.shape; s2 = float(resid @ resid) / (n - p); se = np.sqrt(np.diag(XtXi) * s2)
    return beta[1], se[1], 2 * tdist.sf(abs(beta[1] / se[1]), n - p)


def fit(dose, y, C, PC, k):
    mask = np.isfinite(y) & np.isfinite(dose)
    n = int(mask.sum())
    cols = [np.ones(n), dose[mask], C[mask]]
    if k > 0 and PC is not None:
        cols.append(PC[mask, :k])
    X = np.column_stack(cols).astype(np.float64)
    keep = np.ones(X.shape[1], bool); keep[2:] = X[:, 2:].std(0) > 1e-9
    X = X[:, keep]
    yy = y[mask].astype(np.float64)
    b, se, p = ols_dose(X, yy)
    del X, yy
    return n, b, se, p


_chrom_cache = {}


def load_chrom(chrom, center):
    key = chrom
    if key in _chrom_cache:
        offs, cp, cv, vct, c0 = _chrom_cache[key]
        if abs(center - c0) <= 1_000_000:
            return offs, cp, cv, vct
    lo, hi = center - 6_000_000, center + 6_000_000
    offsets = [0]; cpos = []; cvidx = []; vidx = 0
    with open(f"{PGEN}/acaf_threshold.{chrom}.pvar") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t", 5)
            alts = c[4].split(",")
            offsets.append(offsets[-1] + 1 + len(alts))
            pos = int(c[1])
            if len(alts) == 1 and lo <= pos <= hi:
                cpos.append(pos); cvidx.append(vidx)
            vidx += 1
    offs = np.array(offsets, dtype=np.uintp)
    cp = np.array(cpos, dtype=np.int64); cv = np.array(cvidx, dtype=np.int64)
    _chrom_cache[key] = (offs, cp, cv, vidx, center)
    print(f"    [{chrom}] scanned {vidx} variants; {len(cp)} biallelic in +/-6Mb of {center}", flush=True)
    return offs, cp, cv, vidx


def build_pcs(chrom, center, exclude_pos, W_mb, flank_kb, sel, Npgen, Mmax, pckeep):
    offs, cpos, cvidx, vct = load_chrom(chrom, center)
    W = W_mb * 1_000_000
    inwin = np.abs(cpos - center) <= W
    pos_w, vidx_w = cpos[inwin], cvidx[inwin]
    if flank_kb > 0:
        fk = flank_kb * 1000
        keepm = np.ones(len(pos_w), bool)
        for ep in exclude_pos:
            keepm &= np.abs(pos_w - ep) > fk
    else:
        exset = set(exclude_pos)
        keepm = np.array([p not in exset for p in pos_w], bool)
    pos_w, vidx_w = pos_w[keepm], vidx_w[keepm]
    if len(vidx_w) > Mmax:
        idx = np.linspace(0, len(vidx_w) - 1, Mmax).round().astype(int)
        vidx_w = vidx_w[idx]
    N = len(sel)
    reader = pgenlib.PgenReader(f"{PGEN}/acaf_threshold.{chrom}.pgen".encode(),
                                raw_sample_ct=Npgen, variant_ct=vct, allele_idx_offsets=offs)
    buf = np.empty(Npgen, dtype=np.int32)
    cols = []
    for vi in vidx_w:
        reader.read(int(vi), buf, allele_idx=1)
        g = buf[sel].astype(np.float32)
        g[g < 0] = np.nan
        m = np.nanmean(g)
        if not np.isfinite(m):
            continue
        maf = min(m / 2, 1 - m / 2)
        s = np.nanstd(g)
        if maf < 0.01 or not s > 0:
            continue
        g = np.where(np.isnan(g), m, g)
        cols.append(((g - m) / s).astype(np.float32))
    reader.close(); del buf
    if len(cols) < 2:
        return None, 0
    G = np.column_stack(cols); del cols; gc.collect()
    M = G.shape[1]
    S = G.T @ G
    evals, evecs = np.linalg.eigh(S)
    K = min(pckeep, M)
    Vtop = evecs[:, ::-1][:, :K]
    scores = (G @ Vtop).astype(np.float32); del G, S, evecs; gc.collect()
    return scores, M


def read_pgen_iids():
    with open(f"{PGEN}/acaf_threshold.chr1.psam") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ic = hdr.index("#IID") if "#IID" in hdr else (hdr.index("IID") if "IID" in hdr else 1)
        return [l.rstrip("\n").split("\t")[ic] for l in f]


def load_mode(mode):
    import pandas as pd
    covfile = {"eur": f"{COVDIR}/covariates_v3_eur.tsv",
               "all": f"{COVDIR}/covariates_v3_all_ancestries.tsv"}[mode]
    base = ["age_at_recruitment_spline_1", "age_at_recruitment_spline_2", "age_at_recruitment_spline_3",
            "sex", *[f"gPC{i}" for i in range(1, 17)]]
    use = ["IID", *base] + (["ancestry_pred"] if mode == "all" else [])
    cov = pd.read_csv(covfile, sep="\t", dtype={"IID": str}, usecols=lambda c: c in use)
    cov["sex"] = cov["sex"].map({"Male": 1.0, "Female": 0.0})
    for c in base:
        cov[c] = pd.to_numeric(cov[c], errors="coerce")
    cov = cov.dropna(subset=base).reset_index(drop=True)
    iids = cov["IID"].tolist(); idx = {i: k for k, i in enumerate(iids)}; N = len(iids)
    C = cov[base].to_numpy(np.float64)
    if mode == "all":
        C = np.hstack([C, pd.get_dummies(cov["ancestry_pred"], drop_first=True).astype(np.float64).to_numpy()])
    del cov; gc.collect()
    dose = {}
    with gzip.open(f"{R6}/sv_dosages_222.tsv.gz", "rt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        col = {s: hdr.index(s) for s in HIT_SVS}
        for s in HIT_SVS:
            dose[s] = np.full(N, np.nan, np.float32)
        for line in f:
            p = line.rstrip("\n").split("\t"); i = idx.get(p[0])
            if i is not None:
                for s in HIT_SVS:
                    dose[s][i] = float(p[col[s]])
    lv = {lc: np.full(N, np.nan, np.float32) for lc in ALL_LOINCS}
    for r in csv.DictReader(open(f"{R6}/measurement_medians_222.tsv"), delimiter="\t"):
        i = idx.get(r["IID"])
        if i is not None and r["loinc"] in lv:
            lv[r["loinc"]][i] = float(r["median_value"])
    lirnt = {lc: irnt(v.astype(float)) for lc, v in lv.items()}; del lv; gc.collect()
    return iids, idx, C, N, dose, lirnt


def predictor_positions():
    d = {}
    for r in csv.DictReader(open(f"{H}/sv_weights_long_222.tsv"), delimiter="\t"):
        if r["sv_locus"] in HIT_SVS:
            d.setdefault(r["sv_locus"], []).append(int(r["pos"]))
    return d


HDR = ("sv\ttrait\tloinc\tcv_r2\tdisc_sign\tmode\twindow_mb\tflank_kb\tn_pca_snps\tk_pcs\tn\t"
       "beta_m0\tse_m0\tp_m0\tbeta_m1\tse_m1\tp_m1\tpct_atten\taou_sign_m1\tref_beta_step8\n")


def main(mode):
    print(f"=== SNP-conditioning mode={mode} ===", flush=True)
    pgen_iids = read_pgen_iids(); Npgen = len(pgen_iids)
    ppos = predictor_positions()
    centers = {sv: (min(ppos[sv]) + max(ppos[sv])) // 2 for sv in HIT_SVS}
    iids, idx, C, N, dose, lirnt = load_mode(mode)
    pgen_pos = {i: k for k, i in enumerate(pgen_iids)}
    missing = [i for i in iids if i not in pgen_pos]
    assert not missing, f"{len(missing)} mode IIDs absent from pgen psam"
    sel = np.array([pgen_pos[i] for i in iids], dtype=np.int64)
    del pgen_pos, pgen_iids; gc.collect()

    if mode == "eur":
        runs = [(w, 0) for w in WINDOWS_MB] + [(CENTRAL_W, FLANK_KB)]
    else:
        runs = [(CENTRAL_W, 0)]
    outpath = f"{OUT}/conditioning_{mode}.tsv"
    with open(outpath, "w") as fo:
        fo.write(HDR)
        for scfg in SVS:
            sv, chrom, cvr2 = scfg["sv"], scfg["chrom"], scfg["cv_r2"]
            center = centers[sv]; exclude = ppos[sv]
            dsv = dose[sv].astype(np.float64)
            m0 = {}
            for t in scfg["traits"]:
                y = lirnt[t["loinc"]]
                n0, b0, se0, p0 = fit(dsv, y, C, None, 0)
                m0[t["name"]] = (n0, b0, se0, p0)
                ref = REF_EUR.get((sv, t["name"]))
                tag = ""
                if mode == "eur" and ref is not None and abs(b0 - ref) > 5e-3:
                    tag = f"  !! M0 mismatch ref={ref}"
                print(f"  {sv[:26]} {t['name']:10} M0 n={n0} b={b0:+.5f} p={p0:.2e}{tag}", flush=True)
            for (W, flank) in runs:
                PC, nsnp = build_pcs(chrom, center, exclude, W, flank, sel, Npgen, MMAX[mode], PCKEEP[mode])
                if PC is None:
                    print(f"  {sv[:26]} W={W} flank={flank}: no PCs", flush=True); continue
                ks = [CENTRAL_K] if (mode != "eur" or flank > 0) else PCS
                ks = sorted({min(k, PC.shape[1]) for k in ks})
                for t in scfg["traits"]:
                    y = lirnt[t["loinc"]]; n0, b0, se0, p0 = m0[t["name"]]
                    for kk in ks:
                        n1, b1, se1, p1 = fit(dsv, y, C, PC, kk)
                        atten = (abs(b0) - abs(b1)) / abs(b0) * 100 if b0 != 0 else float("nan")
                        asign = "+" if b1 > 0 else "-"
                        ref = REF_EUR.get((sv, t["name"]), "")
                        fo.write(f"{sv}\t{t['name']}\t{t['loinc']}\t{cvr2}\t{t['disc']}\t{mode}\t{W}\t{flank}\t"
                                 f"{nsnp}\t{kk}\t{n1}\t{b0:.5g}\t{se0:.5g}\t{p0:.3g}\t{b1:.5g}\t{se1:.5g}\t"
                                 f"{p1:.3g}\t{atten:.1f}\t{asign}\t{ref}\n")
                        print(f"    {sv[:22]} {t['name']:10} W={W} fl={flank} k={kk:2d} "
                              f"nsnp={nsnp:3d}  b:{b0:+.4f}->{b1:+.4f}  atten={atten:5.1f}%  p:{p0:.1e}->{p1:.1e}",
                              flush=True)
                del PC
                gc.collect()
            fo.flush()
    print(f"wrote {outpath}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "eur")
