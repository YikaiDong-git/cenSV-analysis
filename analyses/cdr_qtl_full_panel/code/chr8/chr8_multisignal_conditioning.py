#!/usr/bin/env python
"""Does ONE haplotype explain the whole chr8 CDR-position locus?"""
import numpy as np, pandas as pd, os, glob, struct, json
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
CDR = f"{ROOT}/analyses/cdr_qtl"
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
CHROM, PHENO = "chr8", "cdr_pos_offset"
OUT = f"{WU}/output/chr8_multisignal"
os.makedirs(OUT, exist_ok=True)


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n))
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


def rint(a):
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def fit(y, cols):
    X = np.column_stack([np.ones(len(y))] + [c if c.ndim == 2 else c[:, None] for c in cols])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = max(1, len(y) - np.linalg.matrix_rank(X))
    cov = np.linalg.pinv(X.T @ X) * (resid ** 2).sum() / df
    se = np.sqrt(np.maximum(0, np.diag(cov)))
    t = beta / np.where(se == 0, np.nan, se)
    p = 2 * (1 - stats.t.cdf(np.abs(t), df))
    return beta, se, t, p


ph = pd.read_csv(f"{CDR}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
p8 = ph[ph.chm13_chrom == CHROM]
pheno_ind = p8.groupby("sample")[PHENO].mean()
bc_map = p8.groupby("sample")["basecaller"].agg(
    lambda v: 1.0 if v.astype(str).str.contains("1041").any() else 0.0)
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], index=gids)

G = np.load(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_variants.tsv", sep="\t")
sv_samples = open(f"{WU}/geno/perchrom/{CHROM}/{CHROM}_sv_samples.txt").read().split()

cands = sorted(glob.glob(f"{ROOT}/snp_extraction/midpoint1000*"))
PANEL = cands[0]
Gs = np.load(f"{PANEL}/{CHROM}_snp_genotypes.npz", allow_pickle=True)
Gs = Gs[list(Gs.keys())[0]]
snp_samples_f = f"{PANEL}/{CHROM}_snp_samples.txt"
snp_samples = open(snp_samples_f).read().split() if os.path.exists(snp_samples_f) else sv_samples
snpmeta = pd.read_csv(f"{PANEL}/{CHROM}_snp_variants.tsv", sep="\t")
if Gs.shape[0] != len(snp_samples) and Gs.shape[1] == len(snp_samples):
    Gs = Gs.T
print(f"panel={os.path.basename(PANEL)}  SNP matrix {Gs.shape}  meta rows {len(snpmeta)}")

keep = sorted(set(pheno_ind.dropna().index) & set(bc_map.index) & set(PCs.index)
              & set(sv_samples) & set(snp_samples))
svi = {s: i for i, s in enumerate(sv_samples)}
sni = {s: i for i, s in enumerate(snp_samples)}
y = rint(pheno_ind.loc[keep].values)
PCm = PCs.loc[keep].values
bc = bc_map.loc[keep].values
COV = np.column_stack([PCm, bc[:, None]]) if bc.std() > 1e-10 else PCm
n = len(keep)
print(f"n = {n}")

LEADS = [("DEL", 113, 44300221), ("INS", 59, 44302147), ("DEL", 347, 44460941), ("INS", 1868, 44247324)]


def sv_vec(svtype, svlen, pos, which=0):
    rows = vt[(vt.svtype == svtype) & (vt.svlen == svlen) & (vt.pos == pos)]
    if not len(rows):
        return None, None
    out = []
    for idx in rows.index:
        v = np.array([G[idx, svi[s]] for s in keep], float)
        out.append((v, vt.loc[idx]))
    uniq, seen = [], set()
    for v, meta in out:
        k = v.tobytes()
        if k not in seen:
            seen.add(k); uniq.append((v, meta))
    return uniq, rows


u, _ = sv_vec(*LEADS[0])
gate = None
for v, meta in u:
    b, se, t, p = fit(y, [v, COV])
    if abs(t[1] + 6.17) < 0.05:
        gate = (v, meta, t[1])
print(f"\n=== A. GATE: DEL-113 published t=-6.17 -> reproduced t={gate[2]:.4f} "
      f"({'PASS' if gate else 'FAIL'})")
assert gate is not None, "gate failed - design does not match the published scan"
del113 = gate[0]

snpG = np.array([[Gs[sni[s], j] for j in range(Gs.shape[1])] for s in keep], float)
poly = np.nanstd(snpG, axis=0) > 0
cols = np.where(poly)[0]
r_tag = np.array([abs(np.corrcoef(del113, snpG[:, j])[0, 1]) for j in cols])
best = cols[int(np.nanargmax(r_tag))]
tag = snpG[:, best]
tag_name = None
if len(snpmeta) >= Gs.shape[1]:
    row = snpmeta.iloc[best]
    tag_name = ":".join(str(row[c]) for c in snpmeta.columns[:4] if c in snpmeta)
print(f"    tag SNP col {best}  |r| to DEL-113 = {np.nanmax(r_tag):.4f}   ({tag_name})")

res = []
print("\n=== B/C/D. per surviving genotype-group: tag conditioning, own-best tag, and SYMMETRY ===")
print(f"{'SV':22s} {'af':>6s} {'|r|tag':>7s} | {'alone t':>8s} {'alone p':>9s} | "
      f"{'joint p(SV)':>11s} {'joint p(TAG)':>12s} {'SEinfl':>7s} | {'ownTag|r|':>9s} {'own joint p':>11s}")
for svtype, svlen, pos in LEADS:
    u, _ = sv_vec(svtype, svlen, pos)
    if not u:
        print(f"  {svtype}-{svlen}@{pos}: not found"); continue
    for gi, (v, meta) in enumerate(u):
        rt = abs(np.corrcoef(v, tag)[0, 1])
        _, _, t1, p1 = fit(y, [v, COV])
        b2, se2, t2, p2 = fit(y, [v, tag, COV])
        _, se1, _, _ = fit(y, [v, COV])
        se_infl = se2[1] / se1[1]
        r_own = np.array([abs(np.corrcoef(v, snpG[:, j])[0, 1]) for j in cols])
        bo = cols[int(np.nanargmax(r_own))]
        _, _, t3, p3 = fit(y, [v, snpG[:, bo], COV])
        lab = f"{svtype}-{svlen}@{pos}" + (f"#{gi+1}" if len(u) > 1 else "")
        print(f"  {lab:20s} {v.mean()/2:6.3f} {rt:7.3f} | {t1[1]:8.2f} {p1[1]:9.2e} | "
              f"{p2[1]:11.3f} {p2[2]:12.3f} {se_infl:7.2f} | {np.nanmax(r_own):9.3f} {p3[1]:11.3f}")
        res.append(dict(sv=lab, af=float(v.mean()/2), r_tag=float(rt), t_alone=float(t1[1]),
                        p_alone=float(p1[1]), p_joint_sv=float(p2[1]), p_joint_tag=float(p2[2]),
                        se_inflation=float(se_infl), r_own_tag=float(np.nanmax(r_own)),
                        p_joint_own=float(p3[1])))

print("\n=== E. JOINT SV MODEL: DEL-113 + INS-1868 (low mutual r^2) ===")
u1868, _ = sv_vec("INS", 1868, 44247324)
u347, _ = sv_vec("DEL", 347, 44460941)
ins1868 = u1868[0][0]; del347 = u347[0][0]
print(f"  r^2(DEL-113, INS-1868) = {np.corrcoef(del113, ins1868)[0,1]**2:.3f}")
print(f"  r^2(DEL-113, DEL-347)  = {np.corrcoef(del113, del347)[0,1]**2:.3f}")
for lab, pair in [("DEL-113 + INS-1868", [del113, ins1868]), ("DEL-113 + DEL-347", [del113, del347])]:
    b, se, t, p = fit(y, pair + [COV])
    print(f"  {lab:22s} -> p(first)={p[1]:.3e}  p(second)={p[2]:.3e}")
b, se, t, p = fit(y, [del113, ins1868, tag, COV])
print(f"  DEL-113 + INS-1868 + TAG -> p(DEL-113)={p[1]:.3f}  p(INS-1868)={p[2]:.3e}  p(TAG)={p[3]:.3f}")

print("\n=== F. RESIDUAL SCAN: all chr8 SVs after removing covariates + the tag SNP ===")
Xr = np.column_stack([np.ones(n), tag, COV])
yr = y - Xr @ np.linalg.lstsq(Xr, y, rcond=None)[0]
hits = []
for idx in range(len(vt)):
    v = np.array([G[idx, svi[s]] for s in keep], float)
    if v.std() == 0:
        continue
    vr = v - Xr @ np.linalg.lstsq(Xr, v, rcond=None)[0]
    if vr.std() == 0:
        continue
    r = np.corrcoef(vr, yr)[0, 1]
    df = n - Xr.shape[1] - 1
    t = r * np.sqrt(df) / np.sqrt(max(1 - r ** 2, 1e-30))
    hits.append((float(abs(t)), float(2 * (1 - stats.t.cdf(abs(t), df))), idx))
hits.sort(reverse=True)
print(f"  chr8 SVs tested: {len(hits)}   Bonferroni within chr8 (0.05/{len(hits)}) = {0.05/len(hits):.2e}")
for at, pv, idx in hits[:8]:
    r_ = vt.loc[idx]
    print(f"    {r_.svtype}-{r_.svlen}@{r_.pos} af={r_.af:.3f}  |t|={at:.2f} p={pv:.2e}"
          + ("  <-- survives chr8-wide Bonferroni" if pv < 0.05 / len(hits) else ""))

json.dump(res, open(f"{OUT}/chr8_multisignal_conditioning.json", "w"), indent=2)
pd.DataFrame(res).to_csv(f"{OUT}/chr8_multisignal_conditioning.tsv", sep="\t", index=False,
                         float_format="%.6g")
print(f"\nwrote {OUT}/chr8_multisignal_conditioning.tsv")
