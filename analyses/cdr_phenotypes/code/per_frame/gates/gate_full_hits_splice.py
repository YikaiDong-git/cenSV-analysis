#!/usr/bin/env python
"""SECOND-PASS gates on the full-substrate raw hits."""
import pandas as pd, numpy as np, os, sys, struct
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

_D = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_D, os.path.dirname(_D)]
import sys as _sys
_sys.path[:0] = [CENSV_ROOT + "/analyses/cdr_qtl_full_panel/code"]
from phenotype_family import PRIMARY, denominators, alpha, headline

ROOT = CENSV_ROOT + ""
WU = f"{ROOT}/analyses/cdr_qtl_full_panel"
CAN = CENSV_ROOT + "/analyses/cdr_phenotypes/output/compat/compat_root_splice"
PED = f"{ROOT}/resources/1000G/20130606_g1k_3202_samples_ped_population.txt"
GRM = f"{ROOT}/subproject_methylation/results/greml/pericen_1m_all"
HITS = sys.argv[1] if len(sys.argv) > 1 else f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/gate_rederive/agg_splice/raw_hits_1e7.tsv"
TAG = sys.argv[2] if len(sys.argv) > 2 else "1e7"

LENGTH_NORMALISED = {"cdr_pos_offset", "cdr_pos_offset_deepest", "cdr_fraction", "cdrs_per_mb",
                     "cdr_pos_var", "cdr_pos_iqr"}


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
    return stats.norm.ppf((stats.rankdata(a) - 0.5) / len(a))


def ols_slope(y, x, extra=None):
    cols = [np.ones(len(y)), np.asarray(x, float)]
    if extra is not None:
        e = np.asarray(extra, float)
        cols.append(e if e.ndim > 1 else e[:, None])
    X = np.column_stack(cols)
    m = np.isfinite(X).all(axis=1) & np.isfinite(y)
    Xm, ym, n = X[m], np.asarray(y, float)[m], int(m.sum())
    k = Xm.shape[1]
    if n <= k + 1 or Xm[:, 1].std() == 0:
        return None
    beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
    resid = ym - Xm @ beta
    se = np.sqrt(np.diag(np.linalg.pinv(Xm.T @ Xm) * (resid ** 2).sum() / (n - k)))
    t = beta[1] / se[1]
    return dict(beta=float(beta[1]), se=float(se[1]), t=float(t),
                p=float(2 * (1 - stats.t.cdf(abs(t), n - k))), n=n)


G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]].astype(float)
vt = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
samp = np.array(open(f"{WU}/geno/full_samples.txt").read().split())
loc2row = {l: i for i, l in enumerate(vt.sv_locus)}
PH_CORR = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes.tsv", sep="\t")
PH_DEEP = pd.read_csv(f"{CAN}/output/per_hap_chrom_phenotypes_deepest.tsv", sep="\t")
ped = pd.read_csv(PED, sep=r"\s+")
sc = [c for c in ped.columns if c.lower() in ("sampleid", "sample")][0]
spc = [c for c in ped.columns if "super" in c.lower()][0]
sxc = [c for c in ped.columns if c.lower() == "sex"][0]
superpop_map = ped.set_index(sc)[spc].to_dict()
sex_map = ped.set_index(sc)[sxc].to_dict()
K, gids = read_grm(GRM)
ev, evec = np.linalg.eigh(K)
PCs = pd.DataFrame(evec[:, np.argsort(ev)[::-1][:10]], columns=[f"PC{i+1}" for i in range(10)], index=gids)


def cell(chrom, pheno):
    src = PH_DEEP if pheno == "cdr_pos_offset_deepest" else PH_CORR
    s = src[src.chm13_chrom == chrom]
    if pheno not in s.columns:
        return None
    agg = s.groupby("sample").agg(y=(pheno, "mean"))
    b = PH_CORR[PH_CORR.chm13_chrom == chrom].groupby("sample").agg(
        aa_s=("active_start", "mean"), aa_e=("active_end", "mean"))
    return agg.join(b)


hits = pd.read_csv(HITS, sep="\t")
if not len(hits):
    print("no raw hits to gate"); sys.exit(0)
leads = (hits.sort_values("p_emp_perm")
         .groupby(["chrom", "geno_group", "phenotype"], as_index=False)
         .agg(sv_locus=("sv_locus", "first"), pos=("pos", "first"), svtype=("svtype", "first"),
              svlen=("svlen", "first"), af=("af", "first"), status=("status", "first"),
              p_best=("p_emp_perm", "min"), t_best=("t_obs", "first"), n_loci=("sv_locus", "nunique"))
         ).sort_values("p_best")
vt_ = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t")
FRAMINGS = denominators(vt_)
HL_LAB, HL_N, HL_A = headline(vt_)
for nm, n in FRAMINGS:
    k = ("row" if nm.startswith("row") else "uniq") + ("_ADMIT" if "ADMIT" in nm else "_FULL")
    leads[f"pass_{k}"] = leads.p_best <= alpha(n)
leads["pass_headline"] = leads.p_best <= HL_A
PASSCOLS = [c for c in leads.columns if c.startswith("pass_")]
print(f"raw hits: {len(hits)} tests -> {len(leads)} (chrom, geno_group, phenotype) leads "
      f"over {leads.phenotype.nunique()} phenotypes")
print(f"gate-coverage alpha = {max(alpha(n) for _, n in FRAMINGS):.3e}; "
      f"HEADLINE alpha = {HL_A:.3e} ({HL_LAB}, n={HL_N}) -> "
      f"{int(leads.pass_headline.sum())} of these leads clear the headline\n")

results = []
for r in leads.itertuples():
    lennorm = r.phenotype in LENGTH_NORMALISED
    print("=" * 92)
    print(f"=== {r.chrom} {r.svtype}-{r.svlen}@{r.pos} af={r.af:.3f} | {r.phenotype}"
          f"{' [length-normalised]' if lennorm else ' [not length-normalised]'} p={r.p_best:.2e} "
          f"status={r.status}")
    ag = cell(r.chrom, r.phenotype)
    if ag is None:
        print("  phenotype column absent - skipped"); continue
    i = loc2row[r.sv_locus]
    dos = pd.Series(G[i, :], index=samp).replace(-1, np.nan)
    keep = sorted(set(ag.dropna(subset=["y"]).index) & set(dos.dropna().index) & set(PCs.index))
    if len(keep) < 30:
        print(f"  n={len(keep)} <30 - skipped"); continue
    d = dos.loc[keep].values
    rec = dict(chrom=r.chrom, pos=int(r.pos), svtype=r.svtype, svlen=int(r.svlen), af=float(r.af),
               phenotype=r.phenotype, length_normalised=lennorm, sv_locus=r.sv_locus,
               geno_group=int(r.geno_group), status=r.status, p_best=float(r.p_best), n=len(keep),
               primary=r.phenotype in PRIMARY)
    rec.update({c: bool(getattr(r, c)) for c in PASSCOLS})

    rs, ps = stats.spearmanr(d, ag.loc[keep, "aa_s"].values, nan_policy="omit")
    re_, pe = stats.spearmanr(d, ag.loc[keep, "aa_e"].values, nan_policy="omit")
    rl, pl = stats.spearmanr(d, (ag.loc[keep, "aa_e"] - ag.loc[keep, "aa_s"]).values, nan_policy="omit")
    stable = (ps > 0.05 and pe > 0.05 and pl > 0.05)
    bverdict = ("GENUINE_WITHIN_ARRAY" if stable else "BOUNDARY_CONFOUNDED") if lennorm else \
               ("boundary_stable" if stable else "boundary_moves_NOT_DECISIVE")
    print(f"  boundary: start r={rs:+.2f} p={ps:.1e} [{'STABLE' if ps>0.05 else 'MOVES'}]  "
          f"end r={re_:+.2f} p={pe:.1e} [{'STABLE' if pe>0.05 else 'MOVES'}]  "
          f"LEN r={rl:+.2f} p={pl:.1e} [{'STABLE' if pl>0.05 else 'MOVES'}] -> {bverdict}")
    rec.update(bnd_start_r=float(rs), bnd_start_p=float(ps), bnd_end_r=float(re_), bnd_end_p=float(pe),
               bnd_len_r=float(rl), bnd_len_p=float(pl), boundary_verdict=bverdict)

    yv = ag.loc[keep, "y"].values
    ok = np.isfinite(yv) & np.isfinite(d)
    keep2 = [k for k, f in zip(keep, ok) if f]
    y = rint(yv[ok]); dd = d[ok]; PCm = PCs.loc[keep2].values
    sp = np.array([superpop_map.get(s, "NA") for s in keep2])
    pooled = ols_slope(y, dd, PCm)
    print(f"  pooled(10PC): z={pooled['t']:+.2f} p={pooled['p']:.2e} n={pooled['n']}")
    rows = []
    for s in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        m = sp == s
        if m.sum() < 12 or not (0.05 < dd[m].mean() / 2 < 0.95):
            continue
        rr = ols_slope(rint(yv[ok][m]), dd[m])
        if rr:
            rr["superpop"] = s; rows.append(rr)
    if len(rows) >= 2:
        betas = np.array([x["beta"] for x in rows]); ses = np.array([x["se"] for x in rows])
        w = 1 / ses ** 2
        Q = (w * (betas - (w * betas).sum() / w.sum()) ** 2).sum(); dfQ = len(rows) - 1
        I2 = max(0.0, (Q - dfQ) / Q) if Q > 0 else 0.0
        tau2 = max(0.0, (Q - dfQ) / (w.sum() - (w ** 2).sum() / w.sum()))
        w_re = 1 / (ses ** 2 + tau2)
        z_re = ((w_re * betas).sum() / w_re.sum()) / np.sqrt(1 / w_re.sum())
        flip = bool((np.sign(betas) != np.sign(betas[np.argmax(w)])).any())
        conf = (abs(pooled["t"]) > 1.5 * abs(z_re)) or flip or I2 > 0.5
        anc = "CONFOUND_SUSPECT" if conf else "SURVIVES_ANCESTRY"
        print(f"  META: Q={Q:.2f} I2={I2:.2f} RE z={z_re:+.2f} flip={flip} npops={len(rows)} -> {anc}")
        rec.update(I2=float(I2), meta_z=float(z_re), sign_flip=flip, n_superpops=len(rows),
                   pooled_z=float(pooled["t"]), ancestry_verdict=anc)
    else:
        print(f"  META: only {len(rows)} gateable superpop(s) -> UNGATEABLE_RARE")
        rec.update(n_superpops=len(rows), pooled_z=float(pooled["t"]),
                   ancestry_verdict="UNGATEABLE_RARE")

    if r.chrom == "chrX":
        sx = np.array([sex_map.get(s, np.nan) for s in keep2], float)
        ws = ols_slope(y, dd, np.column_stack([PCm, sx]))
        if ws:
            print(f"  chrX +Sex: z={ws['t']:+.2f} p={ws['p']:.2e} (without Sex z={pooled['t']:+.2f})")
            rec.update(sex_adj_z=float(ws["t"]), sex_adj_p=float(ws["p"]))

    anc_ok = rec.get("ancestry_verdict") == "SURVIVES_ANCESTRY"
    if lennorm:
        rec["overall"] = "SURVIVES" if (stable and anc_ok) else "DEMOTED"
    else:
        rec["overall"] = "SURVIVES" if anc_ok else "DEMOTED"
    print(f"  >>> OVERALL: {rec['overall']}")
    results.append(rec)

out = pd.DataFrame(results)
os.makedirs(f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/gate_rederive/gates_full_splice", exist_ok=True)
p = f"{CENSV_ROOT}/analyses/cdr_phenotypes/output/per_frame/gate_rederive/gates_full_splice/full_hits_gate_{TAG}.tsv"
tmp = f"{p}.tmp.{os.getpid()}"
out.to_csv(tmp, sep="\t", index=False, float_format="%.6g")
os.replace(tmp, p)
print("\n\n=== SUMMARY ===")
for r in results:
    print(f"  {r['chrom']:6s} {r['svtype']:10s}-{r['svlen']:<6d}@{r['pos']:<10d} af={r['af']:.3f} "
          f"{r['phenotype']:24s} p={r['p_best']:.1e}  {r['boundary_verdict']:28s} "
          f"{r.get('ancestry_verdict','-'):18s} -> {r['overall']}")
if len(out):
    print("\nby phenotype x outcome (gate-coverage set):")
    print(pd.crosstab(out.phenotype, out.overall).to_string())
    hl = out[out.pass_headline]
    surv = hl[hl.overall == "SURVIVES"]
    print(f"\n=== AT THE HEADLINE ALPHA ({HL_A:.3e}, {HL_LAB}) ===")
    print(f"  leads clearing the headline: {len(hl)}   surviving every gate: {len(surv)}")
    if len(surv):
        print(f"  surviving chromosomes: {sorted(set(surv.chrom))}")
        print(f"  surviving phenotypes:  {sorted(set(surv.phenotype))}")
        print(f"  distinct genotype-groups among survivors: {surv.geno_group.nunique()}")
        print("\n  " + pd.crosstab(surv.phenotype, surv.chrom).to_string().replace("\n", "\n  "))
    dem = hl[hl.overall == "DEMOTED"]
    if len(dem):
        print(f"\n  demoted at the headline: {len(dem)} on {sorted(set(dem.chrom))}")
print(f"\nwrote {p}")
