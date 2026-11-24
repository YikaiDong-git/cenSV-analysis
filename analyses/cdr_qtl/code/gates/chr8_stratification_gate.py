#!/usr/bin/env python3
"""chr8 DEL-342 -> cdr_pos_offset survivor: ancestry-stratification gate."""
import struct
import numpy as np
import pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
SUB = f"{ROOT}/subproject_methylation"
DISC = f"{ROOT}/analyses/cdr_qtl"
PHENO = f"{DISC}/output/per_hap_chrom_phenotypes.tsv"
GRM = f"{SUB}/results/greml/pericen_1m_all"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
CHROM = "chr8"
TARGET = dict(pos=44456532, svtype="DEL", svlen=342)
PHENO_NAME = "cdr_pos_offset"


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    ids = [l.split()[1] for l in open(f"{prefix}.grm.id")]
    return K, ids


def rint(a):
    r = stats.rankdata(a)
    return stats.norm.ppf((r - 0.5) / len(r))


def ols_slope(y, x, extra=None):
    cols = [np.ones(len(y)), x]
    if extra is not None:
        cols.append(extra)
    X = np.column_stack(cols)
    m = ~np.isnan(y) & ~np.isnan(x)
    if extra is not None:
        m &= ~np.isnan(extra).any(axis=1) if extra.ndim > 1 else ~np.isnan(extra)
    Xm, ym, n = X[m], y[m], int(m.sum())
    k = Xm.shape[1]
    if n <= k + 1 or Xm[:, 1].std() == 0:
        return None
    beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
    resid = ym - Xm @ beta
    sigma2 = (resid ** 2).sum() / (n - k)
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(Xm.T @ Xm)))
    t = beta[1] / se[1]
    p = 2 * (1 - stats.t.cdf(abs(t), n - k))
    return dict(beta=float(beta[1]), se=float(se[1]), t=float(t), p=float(p), n=n)


def main():
    svd = (f"{ROOT}/sv_extraction/decompose/"
           "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200/"
           f"{CHROM}")
    npz = np.load(f"{svd}/{CHROM}_sv_genotypes.npz")
    key = "arr_0" if "arr_0" in npz else list(npz.keys())[0]
    G = npz[key]
    sv_samples = open(f"{svd}/{CHROM}_sv_samples.txt").read().split()
    sv_var = pd.read_csv(f"{svd}/{CHROM}_sv_variants.tsv", sep="\t")
    if G.shape[0] != len(sv_samples) and G.shape[1] == len(sv_samples):
        G = G.T
    cand = sv_var[(sv_var.pos == TARGET["pos"]) & (sv_var.svtype == TARGET["svtype"]) &
                  (sv_var.svlen.between(TARGET["svlen"] - 2, TARGET["svlen"] + 2))]
    print(f"[match] {len(cand)} candidate SVs at {CHROM}:{TARGET['pos']} {TARGET['svtype']} "
          f"~{TARGET['svlen']}bp")
    if cand.empty:
        print("!! no match; nearby DELs:")
        print(sv_var[(sv_var.svtype == "DEL") & sv_var.pos.between(TARGET["pos"]-50, TARGET["pos"]+50)]
              [["pos", "svtype", "svlen", "af"]].to_string(index=False))
        return
    afcol = "af" if "af" in cand.columns else None
    print(cand[["pos", "svtype", "svlen"] + ([afcol] if afcol else [])].to_string(index=False))
    if afcol:
        cand = cand.iloc[(cand[afcol] - 0.245).abs().argsort()]
    sv_row_idx = cand.index[0]
    dos_full = G[:, sv_row_idx].astype(float)
    dos = pd.Series(dos_full, index=sv_samples)
    print(f"[pick] row {sv_row_idx} af={cand.iloc[0][afcol] if afcol else 'NA'}")

    ph = pd.read_csv(PHENO, sep="\t")
    ph8 = ph[ph.chm13_chrom == CHROM]
    y_ind = ph8.groupby("sample")[PHENO_NAME].mean()
    bc = ph8.groupby("sample")["basecaller"].first()

    ped = pd.read_csv(PED, sep=r"\s+")
    col_sample = [c for c in ped.columns if c.lower() in ("sample", "individual", "sampleid", "iid")]
    col_super = [c for c in ped.columns if "super" in c.lower() or c.lower() in ("superpop",)]
    print(f"[ped] cols={list(ped.columns)}")
    sc = col_sample[0] if col_sample else ped.columns[1]
    pc_ = col_super[0] if col_super else ped.columns[-1]
    superpop = ped.set_index(sc)[pc_].to_dict()

    K, gids = read_grm(GRM)
    ev, evec = np.linalg.eigh(K)
    order = np.argsort(ev)[::-1]
    PCs = pd.DataFrame(evec[:, order[:10]], columns=[f"PC{i+1}" for i in range(10)])
    PCs["sample"] = gids
    PCs = PCs.set_index("sample")

    keep = sorted(set(y_ind.index) & set(dos.index) & set(PCs.index))
    print(f"[align] n={len(keep)}")
    y = rint(y_ind.loc[keep].values)
    d = dos.loc[keep].values
    bcv = np.array([1.0 if bc.get(s) == "R1041" else 0.0 for s in keep])
    PCmat = PCs.loc[keep].values
    sp = np.array([superpop.get(s, "NA") for s in keep])

    print(f"\n[dosage] af={d.mean()/2:.3f} std={d.std():.3f}")
    print(f"[strat proxy] r(dosage, PC1) = {np.corrcoef(d, PCmat[:,0])[0,1]:+.3f}")

    print("\n=== POOLED (basecaller covar only) ===")
    r = ols_slope(y, d, bcv)
    print(f"  beta={r['beta']:+.3f} se={r['se']:.3f} z~t={r['t']:+.2f} p={r['p']:.2e} n={r['n']}")
    pooled_simple = r

    print("=== POOLED (10 ancestry PCs + basecaller) ===")
    r = ols_slope(y, d, np.column_stack([PCmat, bcv]))
    print(f"  beta={r['beta']:+.3f} se={r['se']:.3f} z~t={r['t']:+.2f} p={r['p']:.2e} n={r['n']}")
    pooled_pc = r

    print("\n=== WITHIN-SUPERPOP (cdr_pos_offset_RINT ~ dosage) ===")
    rows = []
    for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        m = sp == s
        if m.sum() < 12:
            print(f"  {s}: n={int(m.sum())} (skip, <12)")
            continue
        ys = rint(y_ind.loc[[k for i, k in enumerate(keep) if m[i]]].values)
        ds = d[m]
        rr = ols_slope(ys, ds)
        if rr is None:
            print(f"  {s}: n={int(m.sum())} degenerate")
            continue
        rr["superpop"] = s
        rows.append(rr)
        print(f"  {s}: n={rr['n']:3d} beta={rr['beta']:+.3f} se={rr['se']:.3f} "
              f"t={rr['t']:+.2f} p={rr['p']:.2e}")

    if len(rows) >= 2:
        betas = np.array([r["beta"] for r in rows])
        ses = np.array([r["se"] for r in rows])
        w = 1 / ses ** 2
        beta_fe = (w * betas).sum() / w.sum()
        se_fe = np.sqrt(1 / w.sum())
        z_fe = beta_fe / se_fe
        Q = (w * (betas - beta_fe) ** 2).sum()
        dfQ = len(rows) - 1
        pQ = 1 - stats.chi2.cdf(Q, dfQ)
        I2 = max(0, (Q - dfQ) / Q) if Q > 0 else 0
        tau2 = max(0, (Q - dfQ) / (w.sum() - (w**2).sum()/w.sum()))
        w_re = 1 / (ses ** 2 + tau2)
        beta_re = (w_re * betas).sum() / w_re.sum()
        se_re = np.sqrt(1 / w_re.sum())
        z_re = beta_re / se_re
        print(f"\n=== META across {len(rows)} superpops ===")
        print(f"  fixed-effect:  beta={beta_fe:+.3f} z={z_fe:+.2f} "
              f"p={2*(1-stats.norm.cdf(abs(z_fe))):.2e}")
        print(f"  Cochran Q={Q:.2f} df={dfQ} pQ={pQ:.2e}  I^2={I2:.2f}  tau2={tau2:.4f}")
        print(f"  random-effect: beta={beta_re:+.3f} z={z_re:+.2f} "
              f"p={2*(1-stats.norm.cdf(abs(z_re))):.2e}")
        print(f"\n=== CONFOUND SIGNATURE CHECK ===")
        print(f"  pooled(PC+bc) z={pooled_pc['t']:+.2f}  vs  within-superpop RE-meta z={z_re:+.2f}")
        signs = np.sign(betas)
        print(f"  per-superpop beta signs: {dict(zip([r['superpop'] for r in rows], signs.astype(int)))}")
        flip = (signs != signs[np.argmax(w)]).any()
        verdict = ("CONFOUND-SUSPECT" if (abs(pooled_pc['t']) > 1.5*abs(z_re) or flip or I2 > 0.5)
                   else "SURVIVES-GATE")
        print(f"  sign-flip across superpops: {flip}   I^2>0.5: {I2>0.5}   "
              f"pooled>>meta: {abs(pooled_pc['t'])>1.5*abs(z_re)}")
        print(f"\n  >>> VERDICT: {verdict}")

    import json
    summ = dict(target=f"{CHROM}:{TARGET['pos']}_{TARGET['svtype']}{TARGET['svlen']}",
                pooled_simple=pooled_simple, pooled_pc=pooled_pc,
                within=rows, r_dos_pc1=float(np.corrcoef(d, PCmat[:,0])[0,1]))
    with open(f"{DISC}/output/gates/chr8_stratification_gate.json", "w") as f:
        json.dump(summ, f, indent=2, default=str)
    print(f"\n[written] {DISC}/output/gates/chr8_stratification_gate.json")


if __name__ == "__main__":
    main()
