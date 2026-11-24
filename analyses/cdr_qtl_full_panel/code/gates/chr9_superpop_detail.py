#!/usr/bin/env python
"""Why did only 2 superpopulations enter the ancestry meta for the chr9 cdr_fraction leads?"""
import numpy as np, pandas as pd, struct
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
CAN = f"{ROOT}/analyses/cdr_qtl"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
PHENO = "cdr_fraction"
_g = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="	")
_g = _g[(_g.overall == "SURVIVES") & (_g.chrom == "chr9")]
LEADS = [(r.chrom, int(r.pos), r.svtype, int(r.svlen), r.sv_locus, float(r.af), int(r.geno_group))
         for r in _g.itertuples()]


def read_grm(p):
    n = sum(1 for _ in open(f"{p}.grm.id")); K = np.zeros((n, n))
    with open(f"{p}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{p}.grm.id")]


def rint(a): return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def slope(y, x, extra=None):
    cols = [np.ones(len(y)), np.asarray(x, float)]
    if extra is not None:
        cols.append(np.asarray(extra, float))
    X = np.column_stack(cols)
    if len(y) <= X.shape[1] + 1 or X[:, 1].std() == 0:
        return None
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (r ** 2).sum() / (len(y) - X.shape[1])))
    t = b[1] / se[1]
    return dict(beta=b[1], se=se[1], t=t, p=2 * (1 - stats.t.cdf(abs(t), len(y) - X.shape[1])))


ph = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
p9 = ph[ph.chm13_chrom == "chr9"]
pheno = p9.groupby("sample")[PHENO].mean()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
spmap = ped.set_index(sc)[spc].to_dict()
G = np.load(f"{WU}/geno/perchrom/chr9/chr9_sv_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/perchrom/chr9/chr9_sv_variants.tsv", sep="\t")
svs = open(f"{WU}/geno/perchrom/chr9/chr9_sv_samples.txt").read().split()
svi = {s: i for i, s in enumerate(svs)}

rows = []
for chrom, pos, svt, svlen, locus, af_gate, ggroup in LEADS:
    cand = vt[vt.sv_locus == locus]
    if not len(cand):
        print(f"{locus} not found"); continue
    idx = cand.index[0]
    print("=" * 96)
    print(f"=== {chrom} {svt}-{svlen} @{pos:,}  geno_group={ggroup}  (af in gate table {af_gate:.4f})"
          f"   phenotype = {PHENO} ===")
    dos = pd.Series(G[idx, :], index=svs).replace(-1, np.nan)
    keep = sorted(set(pheno.dropna().index) & set(dos.dropna().index) & set(PCs.index))
    y = rint(pheno.loc[keep].values); d = dos.loc[keep].values
    sp = np.array([spmap.get(s, "NA") for s in keep])
    ov = slope(y, d, PCs.loc[keep].values)
    print(f"  pooled: n={len(keep)}  af={d.mean()/2:.4f}  t={ov['t']:+.2f}  p={ov['p']:.2e}")
    print(f"  {'pop':5s} {'n':>4s} {'carriers':>9s} {'MAF':>7s} {'beta':>8s} {'t':>7s} {'p':>9s}  admitted?")
    for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        m = sp == s
        n = int(m.sum())
        if n == 0:
            continue
        maf = d[m].mean() / 2
        carr = int((d[m] > 0).sum())
        r = slope(rint(pheno.loc[[k for k, f in zip(keep, m) if f]].values), d[m]) if d[m].std() > 0 else None
        why = []
        if n < 12:
            why.append("n<12")
        if not (0.05 < maf < 0.95):
            why.append(f"MAF={maf:.3f}<0.05")
        adm = "ADMITTED" if not why else "dropped: " + ", ".join(why)
        bs = f"{r['beta']:+8.3f} {r['t']:+7.2f} {r['p']:9.2e}" if r else f"{'-':>8s} {'-':>7s} {'-':>9s}"
        print(f"  {s:5s} {n:4d} {carr:9d} {maf:7.4f} {bs}  {adm}")
        rows.append(dict(sv=f"{svt}-{svlen}@{pos}#g{ggroup}", pop=s, n=n, carriers=carr, maf=float(maf),
                         beta=(r["beta"] if r else np.nan), t=(r["t"] if r else np.nan),
                         p=(r["p"] if r else np.nan), admitted=(not why),
                         drop_reason=("" if not why else ", ".join(why))))
    print("  NOTE: strata dropped for MAF are dropped on FREQUENCY, not on effect. Their exclusion says")
    print("        nothing about whether the association holds there - it says the gate cannot test it.")

out = f"{WU}/output/chr9/chr9_superpop_detail.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {out}")
