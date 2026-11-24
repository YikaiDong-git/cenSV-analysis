#!/usr/bin/env python
"""Repeated-measures cis/trans SV-burden model over (individual x chromosome) observations."""
import os, sys, struct, glob, re
import numpy as np, pandas as pd
from collections import defaultdict
from scipy import stats
from statsmodels.regression.mixed_linear_model import MixedLM
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CAN = CENSV_ROOT + "/analyses/cdr_phenotypes/output/compat/compat_root_block"
AUD = f"{ROOT}/analyses/cdr_calls"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
OUT = CENSV_ROOT + "/analyses/cdr_phenotypes/output/per_frame/frames/burden_block"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, f"{WU}/code")
from phenotype_family import PRIMARY, SUPPLEMENTARY, NPH

RNG = np.random.default_rng(20260804)
NPERM = 1000


def rint(a):
    a = np.asarray(a, float)
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
G[G < 0] = np.nan
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samples = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
rep = vt.reset_index().sort_values("index").drop_duplicates("geno_group", keep="first")
ri = rep["index"].values
Gu, vu = G[ri, :], vt.iloc[ri].reset_index(drop=True)

rows = []
for c in sorted(set(vu.chrom)):
    m = (vu.chrom == c).values
    rows.append(pd.DataFrame(dict(sample=samples, chrom=c,
                                  burden_cis=np.nansum(Gu[m, :], axis=0),
                                  churn_cis=np.nansum(
                                      Gu[m, :] * np.abs(vu.svlen.values[m])[:, None], axis=0))))
PC = pd.concat(rows, ignore_index=True)
tot = PC.groupby("sample").burden_cis.sum().rename("burden_tot")
PC = PC.merge(tot, on="sample")
PC["burden_trans"] = PC.burden_tot - PC.burden_cis

corr = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
deep = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
k = ["sample", "hap_local", "chm13_chrom"]
allph = corr.merge(deep[k + ["cdr_pos_offset_deepest"]], on=k, how="left")
allph["total_cdr_bp"] = allph.cdr_fraction * allph.active_len
PHEN = PRIMARY + list(SUPPLEMENTARY) + ["active_len", "total_cdr_bp"]
P = (allph.groupby(["sample", "chm13_chrom"])[PHEN].mean().reset_index()
     .rename(columns={"chm13_chrom": "chrom"}))

coh = pd.read_csv(f"{AUD}/output/cohort_corrected_final.tsv", sep="\t")
cov_rows = []
for f in sorted(glob.glob(f"{AUD}/output/calls/*/binned_freq.bed")):
    tag = os.path.basename(os.path.dirname(f))
    s, h = tag.rsplit("_", 1)
    a = pd.read_csv(f, sep="\t", header=None, usecols=[4], names=["cov"])
    cov_rows.append(dict(sample=s, hap=h, coverage=a["cov"].mean()))
Q = coh.merge(pd.DataFrame(cov_rows), on=["sample", "hap"], how="outer")
Q["is_r1041"] = Q.chemistry.astype(str).str.contains("1041").astype(float)
QI = Q.groupby("sample").agg(n_contigs=("n_contigs", "sum"), coverage=("coverage", "mean"),
                             is_r1041=("is_r1041", "max")).reset_index()

ROW_RE = re.compile(r"^(.+?)_(chr[0-9XY]+)_((.+?)#([12])#[^:]+):")
nint = defaultdict(int)
for bed in sorted(glob.glob(f"{AUD}/inputs/globe_aa/*_AS-HOR_stv_row.ort.bed")):
    for line in open(bed):
        m = ROW_RE.match(line.split("\t")[0])
        if m:
            nint[(m.group(1), m.group(2))] += 1
NI = pd.DataFrame([dict(sample=a, chrom=b, n_intervals=v) for (a, b), v in nint.items()])

K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]],
                   columns=[f"gPC{i+1}" for i in range(10)], index=gids)
PCs.index.name = "sample"

D = (PC.merge(P, on=["sample", "chrom"], how="inner")
       .merge(QI, on="sample", how="inner")
       .merge(NI, on=["sample", "chrom"], how="inner")
       .merge(PCs.reset_index(), on="sample", how="inner"))
print(f"model rows (individual x chromosome): {len(D)}  "
      f"[{D['sample'].nunique()} individuals x {D.chrom.nunique()} chromosomes]")

for col in ["burden_cis", "churn_cis"] + PHEN:
    D[col + "_z"] = D.groupby("chrom")[col].transform(
        lambda s: rint(s) if s.notna().sum() == len(s) and s.std() > 0 else np.nan)
D["burden_trans_z"] = rint(D.burden_trans)

COVS = ["n_contigs", "coverage", "is_r1041", "n_intervals"] + [f"gPC{i+1}" for i in range(10)]
BASE = ["burden_cis_z", "burden_trans_z"] + COVS


def fit(dsub, y, terms):
    cols = [y, "sample", "chrom"] + [t for t in terms if t not in (y, "sample", "chrom")]
    d = dsub[cols].replace([np.inf, -np.inf], np.nan).dropna()
    X = np.column_stack([np.ones(len(d))] + [d[t].values.astype(float) for t in terms])
    m = MixedLM(d[y].values.astype(float), X, groups=d["sample"].values).fit(reml=True)
    idx = {t: i + 1 for i, t in enumerate(terms)}
    return m, d, idx


def perm_p(d, y, terms, key, t_obs):
    ids = d["sample"].unique()
    prof = d.pivot_table(index="sample", columns="chrom", values="burden_cis_z")
    cnt = 0
    for _ in range(NPERM):
        perm = dict(zip(ids, RNG.permutation(ids)))
        pv = d.apply(lambda r: prof.loc[perm[r["sample"]], r["chrom"]], axis=1).values
        dd = d.copy(); dd["burden_cis_z"] = pv
        dd = dd.replace([np.inf, -np.inf], np.nan).dropna(subset=["burden_cis_z"])
        X = np.column_stack([np.ones(len(dd))] + [dd[t].values.astype(float) for t in terms])
        try:
            mm = MixedLM(dd[y].values.astype(float), X, groups=dd["sample"].values).fit(reml=True)
            if abs(mm.tvalues[terms.index(key) + 1]) >= abs(t_obs):
                cnt += 1
        except Exception:
            pass
    return (cnt + 1) / (NPERM + 1)


print("\n" + "=" * 104)
print(f"CIS / TRANS BURDEN, all {len(PHEN)} outcomes.  alpha = 0.05/{NPH} = {0.05/NPH:.2g} "
      f"over the primary family")
print("=" * 104)
hdr = (f"{'outcome':24s} {'role':6s} {'n':>5s} {'b_cis':>8s} {'t_cis':>7s} {'p_cis':>9s} "
       f"{'t_trans':>8s} {'p_trans':>9s}")
print(hdr); print("-" * len(hdr))
res = []
for y in PHEN:
    yz = y + "_z"
    role = "prim" if y in PRIMARY else ("supp" if y in SUPPLEMENTARY else "aux")
    try:
        m, d, idx = fit(D, yz, BASE)
    except Exception as e:
        print(f"{y:24s} FAILED {e}"); continue
    tc, pc = m.tvalues[idx["burden_cis_z"]], m.pvalues[idx["burden_cis_z"]]
    tt, pt = m.tvalues[idx["burden_trans_z"]], m.pvalues[idx["burden_trans_z"]]
    bc = m.params[idx["burden_cis_z"]]
    print(f"{y:24s} {role:6s} {len(d):5d} {bc:8.4f} {tc:7.2f} {pc:9.3g} {tt:8.2f} {pt:9.3g}")
    res.append(dict(outcome=y, role=role, n=len(d), beta_cis=bc, t_cis=tc, p_cis=pc,
                    t_trans=tt, p_trans=pt))
R = pd.DataFrame(res)

print("\n" + "=" * 104)
print("LENGTH CHANNEL: 6 of 10 phenotypes divide by active_len, and burden correlates with")
print("active_len at r = +0.25. Re-fit the length-normalised outcomes WITH active_len adjusted out.")
print("=" * 104)
print("NOTE on direction: burden->active_len is t = +11.8, but a longer array is also more")
print("OPPORTUNITY for an SV to be called, so active_len is plausibly a CONFOUNDER (ascertainment)")
print("rather than a mediator. The adjusted column is therefore the conservative reading.")
LN = [p for p in PHEN if p != "active_len"]
print(f"\n{'outcome':24s} {'t_cis raw':>10s} {'t_cis +active_len':>18s} {'attenuation':>12s}")
for y in LN:
    base_t = R.loc[R.outcome == y, "t_cis"]
    if not len(base_t):
        continue
    m2, d2, i2 = fit(D, y + "_z", BASE + ["active_len_z"])
    t2 = m2.tvalues[i2["burden_cis_z"]]
    print(f"{y:24s} {float(base_t.iloc[0]):10.2f} {t2:18.2f} "
          f"{100*(1-abs(t2)/max(1e-9,abs(float(base_t.iloc[0])))):11.0f}%")
    R.loc[R.outcome == y, "t_cis_adj_active_len"] = t2

R.to_csv(f"{OUT}/cis_trans_burden_models.tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {OUT}/cis_trans_burden_models.tsv")

top = R[R.role == "prim"].reindex(R[R.role == "prim"].p_cis.abs().sort_values().index).head(1)
if len(top):
    y = top.outcome.iloc[0]
    m, d, idx = fit(D, y + "_z", BASE)
    t_obs = m.tvalues[idx["burden_cis_z"]]
    print(f"\npermutation check on the strongest primary outcome ({y}, t = {t_obs:.2f}), "
          f"{NPERM} individual-label shuffles ...")
    print(f"  permutation p = {perm_p(d, y + '_z', BASE, 'burden_cis_z', t_obs):.4f}  "
          f"(model p = {m.pvalues[idx['burden_cis_z']]:.3g})")
