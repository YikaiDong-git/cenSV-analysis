#!/usr/bin/env python
import pandas as pd, numpy as np, os, sys, struct
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
CAN = f"{ROOT}/analyses/cdr_qtl"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
SEED = 42


def read_grm(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    K = np.zeros((n, n), np.float32)
    with open(f"{prefix}.grm.bin", "rb") as f:
        for i in range(n):
            for j in range(i + 1):
                v = struct.unpack("f", f.read(4))[0]
                K[i, j] = v; K[j, i] = v
    return K, [l.split()[1] for l in open(f"{prefix}.grm.id")]


def rint(a):
    r = stats.rankdata(a)
    return stats.norm.ppf((r - 0.5) / len(r))


def resid(Y, X):
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X @ b


G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}
phc = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
phd = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample", "individual", "iid")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
sxc = [c for c in ped.columns if c.lower() == "sex"][0]
superpop = ped.set_index(sc)[spc].to_dict()
sexmap = ped.set_index(sc)[sxc].to_dict()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
order = np.argsort(ev)[::-1]
PCs = pd.DataFrame(evec[:, order[:10]], columns=[f"PC{i+1}" for i in range(10)])
PCs["sample"] = gids
PCs = PCs.set_index("sample")

gate = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
xrow = vt[(vt.chrom == "chrX") & (vt.pos == 57964988)].sort_values("af", ascending=False).head(1)
extra = []
if len(xrow):
    r = xrow.iloc[0]
    extra.append(dict(chrom="chrX", pos=int(r.pos), svtype=r.svtype, svlen=int(r.svlen), af=float(r.af),
                      sv_locus=r.sv_locus, phenos="cdr_pos_offset_deepest", p_best=np.nan,
                      boundary_verdict="(near-miss, not a hit)", overall="NEAR_MISS"))
leads = pd.concat([gate, pd.DataFrame(extra)], ignore_index=True) if extra else gate


def build_cell(chrom, pheno):
    src = phd if pheno == "cdr_pos_offset_deepest" else phc
    s = src[src.chm13_chrom == chrom]
    agg = s.groupby("sample").agg(y=(pheno, "mean"), aa_s=("active_start", "mean"),
                                  aa_e=("active_end", "mean"))
    bc = phc[phc.chm13_chrom == chrom].groupby("sample")["basecaller"].agg(
        lambda v: 1.0 if (v.astype(str).str.contains("1041").any()) else 0.0)
    agg = agg.join(bc.rename("bc"))
    return agg


print(f"reviewer follow-ups: NPERM={NPERM}  leads={len(leads)}")
results = []
for (chrom, pheno_str), grp in leads.groupby(["chrom", "phenos"]):
    phenos = [p for p in str(pheno_str).split(",") if p]
    for pheno in phenos:
        agg = build_cell(chrom, pheno)
        for r in grp.itertuples():
            i = loc2row[r.sv_locus]
            dos = pd.Series(G[i, :], index=samp).replace(-1, np.nan)
            keep = sorted(set(agg.dropna(subset=["y"]).index) & set(dos.dropna().index) & set(PCs.index))
            if len(keep) < 30:
                continue
            y = rint(agg.loc[keep, "y"].values)
            d = dos.loc[keep].values
            n = len(keep)
            one = np.ones((n, 1))
            PCm = PCs.loc[keep].values
            bc = agg.loc[keep, "bc"].values[:, None]
            if np.nanstd(bc) == 0:
                bc = np.zeros((n, 0))
            bnd = agg.loc[keep, ["aa_s", "aa_e"]].values
            bnd = (bnd - bnd.mean(0)) / np.where(bnd.std(0) == 0, 1, bnd.std(0))
            sx = np.array([sexmap.get(s, np.nan) for s in keep], float)[:, None]
            sets = {"C0_scan_model": np.hstack([one, PCm, bc]),
                    "C1_plus_boundary": np.hstack([one, PCm, bc, bnd])}
            if chrom == "chrX" and np.isfinite(sx).all() and np.nanstd(sx) > 0:
                sets["C2_plus_boundary_sex"] = np.hstack([one, PCm, bc, bnd, sx])
            row = dict(chrom=chrom, pos=int(r.pos), svtype=r.svtype, svlen=int(r.svlen),
                       af=float(r.af), phenotype=pheno, n=n,
                       gate_overall=getattr(r, "overall", ""), p_scan=getattr(r, "p_best", np.nan))
            for nm, X in sets.items():
                yr = resid(y, X); dr = resid(d, X)
                if dr.std() == 0:
                    continue
                yrs = (yr - yr.mean()) / yr.std(); drs = (dr - dr.mean()) / dr.std()
                dfr = n - X.shape[1] - 1
                r_obs = float(drs @ yrs / n)
                t_obs = r_obs * np.sqrt(dfr) / np.sqrt(max(1 - r_obs ** 2, 1e-30))
                rng = np.random.default_rng(SEED)
                cnt = 0; done = 0
                while done < NPERM:
                    b = min(2000, NPERM - done)
                    Yp = np.empty((b, n))
                    for k in range(b):
                        Yp[k] = yrs[rng.permutation(n)]
                    Rp = (Yp @ drs) / n
                    Tp = Rp * np.sqrt(dfr) / np.sqrt(np.maximum(1 - Rp ** 2, 1e-30))
                    cnt += int((np.abs(Tp) >= abs(t_obs)).sum()); done += b
                row[f"t_{nm}"] = round(t_obs, 4)
                row[f"p_{nm}"] = (cnt + 1) / (NPERM + 1)
            results.append(row)
            print(f"  {chrom:6s} {r.svtype}-{r.svlen}@{r.pos} {pheno:24s} n={n} "
                  + "  ".join(f"{k.replace('p_',''):>18s} t={row.get('t_'+k[2:],float('nan')):+.2f} p={v:.1e}"
                              for k, v in row.items() if k.startswith("p_") and k != "p_scan"))

res = pd.DataFrame(results)
os.makedirs(f"{WU}/output/reviewer_followups", exist_ok=True)
p = f"{WU}/output/reviewer_followups/boundary_covariate_sensitivity.tsv"
tmp = f"{p}.tmp.{os.getpid()}"
res.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
os.replace(tmp, p)
print(f"\nwrote {p}")

print("\n=== C. chr8 leverage checks (leave-one-superpop-out, basecaller strata) ===")
lev = []
c8 = leads[(leads.chrom == "chr8")]
agg8 = build_cell("chr8", "cdr_pos_offset")
for r in c8.itertuples():
    i = loc2row[r.sv_locus]
    dos = pd.Series(G[i, :], index=samp).replace(-1, np.nan)
    keep = sorted(set(agg8.dropna(subset=["y"]).index) & set(dos.dropna().index) & set(PCs.index))
    y = rint(agg8.loc[keep, "y"].values); d = dos.loc[keep].values
    sp = np.array([superpop.get(s, "NA") for s in keep])
    bc = agg8.loc[keep, "bc"].values
    X = np.hstack([np.ones((len(keep), 1)), PCs.loc[keep].values, bc[:, None]])
    yr = resid(y, X); dr = resid(d, X)
    base_t = float(np.corrcoef(dr, yr)[0, 1]) * np.sqrt(len(keep) - X.shape[1] - 1)
    line = f"  {r.svtype}-{r.svlen}@{r.pos} af={r.af:.3f} full t~{base_t:+.2f} | LOO:"
    worst = None
    for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        m = sp != s
        if m.sum() < 40:
            continue
        yy = rint(agg8.loc[[k for k, f in zip(keep, m) if f], "y"].values)
        dd = d[m]; XX = X[m]
        yrr = resid(yy, XX); drr = resid(dd, XX)
        tt = float(np.corrcoef(drr, yrr)[0, 1]) * np.sqrt(int(m.sum()) - XX.shape[1] - 1)
        line += f" -{s}:{tt:+.1f}"
        worst = tt if worst is None or abs(tt) < abs(worst) else worst
    for bval, bname in [(0.0, "R941"), (1.0, "R1041")]:
        m = bc == bval
        if m.sum() < 30:
            continue
        yy = rint(agg8.loc[[k for k, f in zip(keep, m) if f], "y"].values)
        XX = np.hstack([np.ones((int(m.sum()), 1)), PCs.loc[[k for k, f in zip(keep, m) if f]].values])
        yrr = resid(yy, XX); drr = resid(d[m], XX)
        tt = float(np.corrcoef(drr, yrr)[0, 1]) * np.sqrt(int(m.sum()) - XX.shape[1] - 1)
        line += f" | {bname}(n={int(m.sum())}):{tt:+.1f}"
    print(line)
    lev.append(dict(sv=f"{r.svtype}-{r.svlen}@{r.pos}", af=float(r.af), t_full=base_t, t_worst_loo=worst))
pd.DataFrame(lev).to_csv(f"{WU}/output/reviewer_followups/chr8_leverage.tsv", sep="\t", index=False,
                         float_format="%.6g")
print(f"wrote {WU}/output/reviewer_followups/chr8_leverage.tsv")
