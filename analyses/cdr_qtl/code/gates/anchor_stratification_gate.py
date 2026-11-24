#!/usr/bin/env python3
import struct, json
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
SVDROOT = (f"{ROOT}/sv_extraction/decompose/"
           "mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200")
PHENO_NAME = "cdr_pos_offset"

TARGETS = [
    ("chr16", 37195737, "INS", 340, "chr16 INS-340 (2x monomer, af0.43)"),
    ("chr5", 49436587, "DEL", 19361, "chr5 DEL-19361 (large, af0.45)"),
    ("chr5", 49456085, "DUP", 6791, "chr5 DUP:TANDEM-6791 (af0.49)"),
    ("chr5", 49436650, "INS", 679, "chr5 INS-679 (4x monomer, af0.36)"),
    ("chr19", 29761459, "DEL", 680, "chr19 DEL-680 (low-freq af0.06)"),
    ("chr19", 29099871, "DEL", 339, "chr19 DEL-339 (low-freq af0.06)"),
]


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
        m &= (~np.isnan(extra).any(axis=1)) if extra.ndim > 1 else (~np.isnan(extra))
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


ph = pd.read_csv(PHENO, sep="\t")
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample", "individual", "iid")][0]
pc_ = [c for c in ped.columns if "super" in c.lower()][0]
superpop_map = ped.set_index(sc)[pc_].to_dict()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
order = np.argsort(ev)[::-1]
PCs = pd.DataFrame(evec[:, order[:10]], columns=[f"PC{i+1}" for i in range(10)])
PCs["sample"] = gids
PCs = PCs.set_index("sample")

results = []
for (chrom, pos, svt, svlen, label) in TARGETS:
    print(f"\n{'='*70}\n=== {label}  [{chrom}:{pos} {svt} {svlen}] ===")
    svd = f"{SVDROOT}/{chrom}"
    npz = np.load(f"{svd}/{chrom}_sv_genotypes.npz")
    key = "arr_0" if "arr_0" in npz else list(npz.keys())[0]
    G = npz[key]
    sv_samples = open(f"{svd}/{chrom}_sv_samples.txt").read().split()
    sv_var = pd.read_csv(f"{svd}/{chrom}_sv_variants.tsv", sep="\t")
    if G.shape[0] != len(sv_samples) and G.shape[1] == len(sv_samples):
        G = G.T
    cand = sv_var[(sv_var.pos == pos) & (sv_var.svtype.astype(str).str.startswith(svt)) &
                  (sv_var.svlen.between(svlen - 3, svlen + 3))]
    if cand.empty:
        print("  !! no SV match; skip")
        results.append(dict(label=label, verdict="NO_MATCH"))
        continue
    afcol = "af" if "af" in cand.columns else None
    ridx = cand.index[0]
    dos = pd.Series(G[:, ridx].astype(float), index=sv_samples)

    ph_c = ph[ph.chm13_chrom == chrom]
    y_ind = ph_c.groupby("sample")[PHENO_NAME].mean()
    keep = sorted(set(y_ind.index) & set(dos.index) & set(PCs.index))
    if len(keep) < 30:
        print(f"  n={len(keep)} <30 skip"); results.append(dict(label=label, verdict="N_LOW")); continue
    y = rint(y_ind.loc[keep].values)
    d = dos.loc[keep].values
    PCmat = PCs.loc[keep].values
    sp = np.array([superpop_map.get(s, "NA") for s in keep])
    rpc1 = float(np.corrcoef(d, PCmat[:, 0])[0, 1])
    pooled = ols_slope(y, d, PCmat)
    print(f"  af={d.mean()/2:.3f}  r(dos,PC1)={rpc1:+.3f}  pooled(10PC) z={pooled['t']:+.2f} p={pooled['p']:.2e} n={pooled['n']}")

    rows = []
    for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        m = sp == s
        if m.sum() < 12:
            continue
        ys = rint(y_ind.loc[[k for i, k in enumerate(keep) if m[i]]].values)
        rr = ols_slope(ys, d[m])
        if rr:
            rr["superpop"] = s; rows.append(rr)
            print(f"    {s}: n={rr['n']:3d} beta={rr['beta']:+.3f} t={rr['t']:+.2f} p={rr['p']:.2e}")
    verdict = "INDETERMINATE"
    if len(rows) >= 2:
        betas = np.array([r["beta"] for r in rows]); ses = np.array([r["se"] for r in rows])
        w = 1 / ses ** 2
        beta_fe = (w * betas).sum() / w.sum(); Q = (w * (betas - beta_fe) ** 2).sum()
        dfQ = len(rows) - 1; I2 = max(0, (Q - dfQ) / Q) if Q > 0 else 0
        tau2 = max(0, (Q - dfQ) / (w.sum() - (w ** 2).sum() / w.sum()))
        w_re = 1 / (ses ** 2 + tau2); beta_re = (w_re * betas).sum() / w_re.sum()
        z_re = beta_re / np.sqrt(1 / w_re.sum())
        signs = np.sign(betas); flip = (signs != signs[np.argmax(w)]).any()
        print(f"  META: Q={Q:.2f} pQ={1-stats.chi2.cdf(Q,dfQ):.2e} I2={I2:.2f}  RE z={z_re:+.2f}")
        print(f"  pooled z={pooled['t']:+.2f} vs within-meta z={z_re:+.2f}  sign-flip={flip}")
        confound = (abs(pooled['t']) > 1.5 * abs(z_re)) or flip or I2 > 0.5
        verdict = "CONFOUND-SUSPECT" if confound else "SURVIVES-GATE"
        print(f"  >>> VERDICT: {verdict}")
        results.append(dict(label=label, chrom=chrom, pos=pos, af=float(d.mean()/2),
                            r_dos_pc1=rpc1, pooled_z=pooled['t'], pooled_p=pooled['p'],
                            within_meta_z=float(z_re), Q=float(Q), I2=float(I2),
                            sign_flip=bool(flip), n_superpops=len(rows), verdict=verdict))
    else:
        print("  too few superpops for meta"); results.append(dict(label=label, verdict="FEW_POPS"))

print("\n\n=== SUMMARY ===")
for r in results:
    print(f"  {r['label']:42s} -> {r.get('verdict')}"
          + (f"  (pooled z={r.get('pooled_z',0):+.1f}, meta z={r.get('within_meta_z',0):+.1f}, "
             f"I2={r.get('I2',0):.2f}, flip={r.get('sign_flip')})" if 'pooled_z' in r else ""))
with open(f"{DISC}/output/gates/anchor_stratification_gate.json", "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"\n[written] {DISC}/output/gates/anchor_stratification_gate.json")
