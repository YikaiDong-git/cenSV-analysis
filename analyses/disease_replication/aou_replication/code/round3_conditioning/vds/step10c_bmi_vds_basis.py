#!/usr/bin/env python3
import os, sys, gc, numpy as np
from step10_snp_conditioning_vds import (load_mode, read_pgen_iids, predictor_positions,
                                         irnt, ols_dose, fit, OUT)

HOME = os.path.expanduser("~")
BED = f"{HOME}/workspace/working/jwein22/centromere_phewas/results/vds_pericentromeric_snps_min_ac1000"
SV = "chr11_sv263_DEL_non_best8_edge500"
CHROM = "chr11"
LOINC = "39156-5"
CVR2 = 0.425

WINDOWS_MB = [1, 2, 5]
PCS = [10, 20, 40, 80]
MMAX = {"eur": 300, "all": 200}
PCKEEP = {"eur": 90, "all": 90}

_CODE2A1 = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)
_LUT = np.empty((256, 4), dtype=np.float64)
for _b in range(256):
    for _k in range(4):
        _LUT[_b, _k] = _CODE2A1[(_b >> (2 * _k)) & 3]


def read_fam_iids(path):
    with open(path) as f:
        return [ln.split()[1] for ln in f]


def scan_bim_window(path, center, W, exclude, flank_kb):
    lo, hi = center - W, center + W
    exset = set(exclude)
    fk = flank_kb * 1000
    keep = []
    vidx = 0
    with open(path) as f:
        for ln in f:
            c = ln.split()
            p = int(c[3])
            if lo <= p <= hi:
                drop = (p in exset)
                if not drop and fk > 0:
                    for ep in exclude:
                        if abs(p - ep) <= fk:
                            drop = True; break
                if not drop:
                    keep.append((vidx, p))
            vidx += 1
    return keep, vidx


def build_pcs_vds(center, exclude, W_mb, flank_kb, sel_vds, Nfam, Mmax, pckeep, bimct_cache):
    pref = f"{BED}/{CHROM}/{CHROM}.pericentromeric_snps.min_ac1000"
    W = W_mb * 1_000_000
    keep, nvar = scan_bim_window(pref + ".bim", center, W, exclude, flank_kb)
    if len(keep) > Mmax:
        idx = np.linspace(0, len(keep) - 1, Mmax).round().astype(int)
        keep = [keep[i] for i in idx]
    N = len(sel_vds)
    have = sel_vds >= 0
    sel_have = sel_vds[have]
    bytes_per_var = (Nfam + 3) // 4
    G = np.empty((N, len(keep)), dtype=np.float32)
    cnt = 0
    with open(pref + ".bed", "rb") as bf:
        assert bf.read(3) == b"\x6c\x1b\x01", "bad .bed magic"
        for (vi, _p) in keep:
            bf.seek(3 + vi * bytes_per_var)
            raw = np.frombuffer(bf.read(bytes_per_var), dtype=np.uint8)
            a1c = _LUT[raw].reshape(-1)[:Nfam]
            g = np.full(N, np.nan, dtype=np.float64)
            g[have] = a1c[sel_have]
            m = np.nanmean(g)
            if not np.isfinite(m):
                continue
            maf = min(m / 2, 1 - m / 2)
            s = np.nanstd(g)
            if maf < 0.01 or not s > 0:
                continue
            g = np.where(np.isnan(g), m, g)
            G[:, cnt] = ((g - m) / s).astype(np.float32)
            cnt += 1
    if cnt < 2:
        return None, 0
    G = G[:, :cnt]
    M = cnt
    S = G.T @ G
    evals, evecs = np.linalg.eigh(S)
    K = min(pckeep, M)
    Vtop = evecs[:, ::-1][:, :K]
    scores = (G @ Vtop).astype(np.float32); del G, S, evecs; gc.collect()
    return scores, M


HDR = ("sv\ttrait\tloinc\tcv_r2\tmode\tbasis\twindow_mb\tflank_kb\tmmax\tn_pca_snps\tk_pcs\tn\t"
       "beta_m0\tse_m0\tp_m0\tbeta_m1\tse_m1\tp_m1\tpct_atten\taou_sign_m1\n")


def main(mode):
    print(f"=== BMI VDS-basis conditioning mode={mode} ===", flush=True)
    ppos = predictor_positions()
    exclude = ppos[SV]
    center = (min(exclude) + max(exclude)) // 2
    iids, idx, C, N, dose, lirnt = load_mode(mode)
    dsv = dose[SV].astype(np.float64)
    y = lirnt[LOINC]

    pref = f"{BED}/{CHROM}/{CHROM}.pericentromeric_snps.min_ac1000"
    fam_iids = read_fam_iids(pref + ".fam"); Nfam = len(fam_iids)
    fam_pos = {i: k for k, i in enumerate(fam_iids)}
    sel_vds = np.array([fam_pos.get(i, -1) for i in iids], dtype=np.int64)
    n_absent = int((sel_vds < 0).sum())
    print(f"  N(mode)={N}  Nfam={Nfam}  mode-iids absent from VDS fam={n_absent}", flush=True)

    n0, b0, se0, p0 = fit(dsv, y, C, None, 0)
    print(f"  M0: n={n0} b={b0:+.5f} p={p0:.3e}", flush=True)

    cells = [(w, 0, MMAX[mode]) for w in WINDOWS_MB]
    cells += [(2, 10, MMAX[mode])]
    outp = f"{OUT}/bmi_vds_basis_{mode}.tsv"
    with open(outp, "w") as fo:
        fo.write(HDR); fo.flush()
        for (W, flank, mmax) in cells:
            PC, nsnp = build_pcs_vds(center, exclude, W, flank, sel_vds, Nfam, mmax, PCKEEP[mode], {})
            if PC is None:
                print(f"  W={W} flank={flank} mmax={mmax}: no PCs", flush=True); continue
            ks = sorted({min(k, PC.shape[1]) for k in PCS})
            for kk in ks:
                n1, b1, se1, p1 = fit(dsv, y, C, PC, kk)
                atten = (abs(b0) - abs(b1)) / abs(b0) * 100 if b0 != 0 else float("nan")
                asign = "+" if b1 > 0 else "-"
                fo.write(f"{SV}\tBMI\t{LOINC}\t{CVR2}\t{mode}\tVDS\t{W}\t{flank}\t{mmax}\t{nsnp}\t{kk}\t{n1}\t"
                         f"{b0:.5g}\t{se0:.5g}\t{p0:.3g}\t{b1:.5g}\t{se1:.5g}\t{p1:.3g}\t{atten:.1f}\t{asign}\n")
                fo.flush()
                print(f"    W={W} fl={flank} mmax={mmax} nsnp={nsnp:4d} k={kk:3d}  "
                      f"b:{b0:+.4f}->{b1:+.4f}  atten={atten:6.1f}%  p:{p0:.1e}->{p1:.1e}", flush=True)
            del PC; gc.collect()
    print(f"wrote {outp}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "eur")
