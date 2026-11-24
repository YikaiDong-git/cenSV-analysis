#!/usr/bin/env python
"""Is the SV-vs-tagSNP pair separable ANYWHERE?"""
import numpy as np, pandas as pd, struct, glob, os
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CDR = f"{ROOT}/analyses/cdr_qtl"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
CHROM, PHENO = "chr8", "cdr_pos_offset"


def read_grm(p):
    n = sum(1 for _ in open(f"{p}.grm.id")); K = np.zeros((n, n))
    with open(f"{p}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]; K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{p}.grm.id")]


def rint(a): return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def fit(y, cols):
    X = np.column_stack([np.ones(len(y))] + [c if c.ndim == 2 else c[:, None] for c in cols])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; df = max(1, len(y) - np.linalg.matrix_rank(X))
    cov = np.linalg.pinv(X.T @ X) * (r ** 2).sum() / df
    se = np.sqrt(np.maximum(0, np.diag(cov))); t = b / np.where(se == 0, np.nan, se)
    return b, se, t, 2 * (1 - stats.t.cdf(np.abs(t), df))


ph = pd.read_csv(f"{CDR}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
p8 = ph[ph.chm13_chrom == CHROM]
pheno = p8.groupby("sample")[PHENO].mean()
bc = p8.groupby("sample")["basecaller"].agg(lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
K, gids = read_grm(GRM); ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
spmap = ped.set_index(sc)[spc].to_dict()

G = np.load(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_variants.tsv", sep="\t")
svs = open(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_samples.txt").read().split()
PANEL = sorted(glob.glob(f"{ROOT}/snp_extraction/midpoint1000*"))[0]
Gs = np.load(f"{PANEL}/{CHROM}_snp_genotypes.npz", allow_pickle=True)
Gs = Gs[list(Gs.keys())[0]]
f = f"{PANEL}/{CHROM}_snp_samples.txt"
sns_ = open(f).read().split() if os.path.exists(f) else svs
if Gs.shape[0] != len(sns_) and Gs.shape[1] == len(sns_):
    Gs = Gs.T
keep = sorted(set(pheno.dropna().index) & set(bc.index) & set(PCs.index) & set(svs) & set(sns_))
svi = {s: i for i, s in enumerate(svs)}; sni = {s: i for i, s in enumerate(sns_)}
y = rint(pheno.loc[keep].values)
COV = np.column_stack([PCs.loc[keep].values, bc.loc[keep].values[:, None]])
row = vt[(vt.svtype == "DEL") & (vt.svlen == 113) & (vt.pos == 44300221)]
sv = None
for idx in row.index:
    v = np.array([G[idx, svi[s]] for s in keep], float)
    _, _, t, _ = fit(y, [v, COV])
    if abs(t[1] + 6.17) < 0.05:
        sv = v
assert sv is not None
snpG = np.array([[Gs[sni[s], j] for j in range(Gs.shape[1])] for s in keep], float)
cols = np.where(np.nanstd(snpG, axis=0) > 0)[0]
rr = np.array([abs(np.corrcoef(sv, snpG[:, j])[0, 1]) for j in cols])
tag = snpG[:, cols[int(np.nanargmax(rr))]]
sp = np.array([spmap.get(s, "NA") for s in keep])

print(f"n={len(keep)}  pooled |r|(SV,tag)={abs(np.corrcoef(sv,tag)[0,1]):.4f}  "
      f"r^2={np.corrcoef(sv,tag)[0,1]**2:.4f}")
print(f"variance inflation factor 1/(1-r^2) = {1/(1-np.corrcoef(sv,tag)[0,1]**2):.1f}\n")
print("=== how many individuals actually DISCORDANT (the entire information to separate them)? ===")
svb = (sv > 0).astype(int); tgb = (tag > 0).astype(int)
disc = int((svb != tgb).sum())
print(f"  carrier-status discordant: {disc}/{len(keep)}  ({100*disc/len(keep):.1f}%)")
print(f"  exact dosage discordant  : {int((sv != tag).sum())}/{len(keep)}")
ct = pd.crosstab(pd.Series(sv, name="SV_dosage"), pd.Series(tag, name="tag_dosage"))
print(ct.to_string())

print("\n=== does the LD break down in any superpopulation? ===")
print(f"{'pop':5s} {'n':>4s} {'|r|':>6s} {'r2':>6s} {'VIF':>6s} | {'SV alone p':>11s} {'joint p(SV)':>11s} {'joint p(tag)':>12s}")
rows = []
for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
    m = sp == s
    if m.sum() < 15:
        continue
    ys = rint(pheno.loc[[k for k, f_ in zip(keep, m) if f_]].values)
    a, b = sv[m], tag[m]
    if a.std() == 0 or b.std() == 0:
        continue
    r = abs(np.corrcoef(a, b)[0, 1]); vif = 1 / max(1 - r ** 2, 1e-9)
    _, _, t1, p1 = fit(ys, [a])
    _, _, t2, p2 = fit(ys, [a, b])
    print(f"{s:5s} {int(m.sum()):4d} {r:6.3f} {r**2:6.3f} {vif:6.1f} | {p1[1]:11.3e} {p2[1]:11.3f} {p2[2]:12.3f}")
    rows.append(dict(pop=s, n=int(m.sum()), r=float(r), r2=float(r ** 2), vif=float(vif),
                     p_sv_alone=float(p1[1]), p_joint_sv=float(p2[1]), p_joint_tag=float(p2[2])))
out = f"{WU}/output/chr8_multisignal/ld_breakdown_by_superpop.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {out}")
print("\nNOTE: a stratum can only adjudicate if its r^2 drops substantially AND it retains power;")
print("with n<=46 per superpopulation neither is likely, which is itself the answer.")
