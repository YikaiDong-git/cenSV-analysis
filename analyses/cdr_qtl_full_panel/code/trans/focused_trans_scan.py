#!/usr/bin/env python
"""FOCUSED TRANS SCAN -- the head-to-head analogue of jwein22's focused trans-meQTL."""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

P1 = CENSV_ROOT + ""
SUB = f"{P1}/subproject_methylation"
WU = f"{P1}/analyses/cdr_qtl_full_panel"
CANON = f"{P1}/analyses/cdr_qtl"
OUT = f"{WU}/output/trans"
sys.path.insert(0, f"{WU}/code")
from phenotype_family import PRIMARY

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 20260805
os.makedirs(OUT, exist_ok=True)

PH_CORR = f"{CANON}/output/per_hap_chrom_phenotypes.tsv"
PH_DEEP = f"{CANON}/output/per_hap_chrom_phenotypes_deepest.tsv"
GRM = f"{SUB}/results/greml/pericen_1m_all"


def rint(a):
    r = stats.rankdata(a)
    return stats.norm.ppf((r - 0.5) / len(r))


def read_gcta_grm_bin(prefix):
    n = sum(1 for _ in open(f"{prefix}.grm.id"))
    vals = np.fromfile(f"{prefix}.grm.bin", dtype=np.float32)
    G = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            G[i, j] = G[j, i] = vals[k]; k += 1
    ids = [line.strip().split() for line in open(f"{prefix}.grm.id")]
    return G, [x[1] for x in ids]


GATE = pd.read_csv(f"{WU}/output/gates/full_hits_gate_1e7.tsv", sep="\t")
CAND = GATE[(GATE.pass_headline) & (GATE.overall == "SURVIVES")].copy()
print(f"headline-surviving rows: {len(CAND)}  distinct sv_locus: {CAND.sv_locus.nunique()}  "
      f"distinct geno_group: {CAND.geno_group.nunique()}")

V = pd.read_csv(f"{WU}/geno/full_variants.tsv", sep="\t", engine="python")
G = np.load(f"{WU}/geno/full_genotypes.npz")
G = G[list(G.keys())[0]]
GS = open(f"{WU}/geno/full_samples.txt").read().split()
if G.shape[0] != len(GS):
    G = G.T
print(f"genotype matrix {G.shape}  samples {len(GS)}  variants {len(V)}")

vi = {s: i for i, s in enumerate(V.sv_locus)}
seen, cands = {}, []
for r in CAND.drop_duplicates("sv_locus").itertuples():
    if r.sv_locus not in vi:
        continue
    d = G[:, vi[r.sv_locus]]
    key = np.nan_to_num(d.astype(np.float64), nan=-9.0).tobytes()
    if key in seen:
        seen[key].append(r.sv_locus); continue
    seen[key] = [r.sv_locus]
    cands.append(dict(sv_locus=r.sv_locus, sv_chrom=r.chrom, pos=r.pos, svtype=r.svtype,
                      svlen=r.svlen, af=r.af, geno_group=r.geno_group, col=vi[r.sv_locus]))
CD = pd.DataFrame(cands)
CD["aliases"] = [";".join(seen[np.nan_to_num(G[:, c].astype(np.float64), nan=-9.0).tobytes()])
                 for c in CD.col]
print(f"\ncandidate SVs after dosage dedup: {len(CD)}")
print(CD[["sv_locus", "sv_chrom", "pos", "svtype", "svlen", "af", "geno_group"]].to_string(index=False))

PC = pd.read_csv(PH_CORR, sep="\t")
PD_ = pd.read_csv(PH_DEEP, sep="\t")
BC = PC.groupby("sample")["basecaller"].first().to_dict()
Gm, gids = read_gcta_grm_bin(GRM)
w, v = np.linalg.eigh(Gm)
PCs = v[:, np.argsort(w)[::-1][:10]]
PCMAP = {s: PCs[i] for i, s in enumerate(gids)}

TARGETS = sorted(set(PC.chm13_chrom.dropna()) | set(PD_.chm13_chrom.dropna()))
print(f"\ntarget chromosomes: {len(TARGETS)}   phenotypes: {len(PRIMARY)}")

FULL = GS
fpos = {s: i for i, s in enumerate(FULL)}
N = len(FULL)

cells = []
for ph in PRIMARY:
    src = PD_ if ph == "cdr_pos_offset_deepest" else PC
    if ph not in src.columns:
        print(f"  [skip] {ph} not in its phenotype file")
        continue
    for tc in TARGETS:
        sub = src[src.chm13_chrom == tc]
        if not len(sub):
            continue
        ind = sub.groupby("sample")[ph].mean().dropna()
        samples = [s for s in ind.index if s in fpos and s in PCMAP]
        if len(samples) < 40:
            continue
        y = rint(ind.loc[samples].values)
        bc = np.array([1.0 if str(BC.get(s, "")) == "R1041" else 0.0 for s in samples])
        cov = [np.ones(len(y)), np.array([PCMAP[s] for s in samples])]
        if bc.std() > 1e-10:
            cov.append(bc[:, None])
        X = np.column_stack(cov)
        df = len(y) - X.shape[1] - 1
        yr = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        if yr.std() < 1e-12:
            continue
        cells.append(dict(pheno=ph, target=tc, samples=samples,
                          idx=np.array([fpos[s] for s in samples]), X=X, df=df,
                          y=(yr - yr.mean()) / yr.std()))
print(f"cells built: {len(cells)}")

tests, rows = [], []
for ci, c in enumerate(cells):
    for r in CD.itertuples():
        d = G[c["idx"], r.col].astype(float)
        m = ~np.isnan(d)
        if m.sum() < 40:
            continue
        x, Xc, yv = d[m], c["X"][m], c["y"][m]
        if x.std() < 1e-10:
            continue
        xr = x - Xc @ np.linalg.lstsq(Xc, x, rcond=None)[0]
        if xr.std() < 1e-10:
            continue
        xr = (xr - xr.mean()) / xr.std()
        yv = (yv - yv.mean()) / yv.std()
        n = len(xr)
        dfc = n - Xc.shape[1] - 1
        rr = float(xr @ yv / n)
        t = rr * np.sqrt(dfc) / np.sqrt(max(1 - rr ** 2, 1e-30))
        tests.append(dict(idx_m=c["idx"][m], xr=xr, y=yv, df=dfc, n=n, sv=r.sv_locus))
        rows.append(dict(sv_locus=r.sv_locus, sv_chrom=r.sv_chrom, pos=r.pos,
                         svtype=r.svtype, af=r.af, geno_group=r.geno_group,
                         phenotype=c["pheno"], target_chrom=c["target"],
                         relation="cis" if c["target"] == r.sv_chrom else "trans",
                         n=n, r=rr, t=t, abs_t=abs(t),
                         p_param=2 * stats.t.sf(abs(t), dfc)))
R = pd.DataFrame(rows)
R["test_i"] = np.arange(len(R))
print(f"\ntests: {len(R)}   cis {int((R.relation=='cis').sum())}   trans {int((R.relation=='trans').sum())}")

groups = {}
for i, t in enumerate(tests):
    groups.setdefault(t["idx_m"].tobytes(), []).append(i)
print(f"distinct individual-sets among tests: {len(groups)}")

rng = np.random.default_rng(SEED)
cnt = np.zeros(len(R), dtype=np.int64)
maxT = np.zeros(NPERM)
maxT_sv = {s: np.zeros(NPERM) for s in CD.sv_locus}
abs_t = R.abs_t.values
BATCH, done = 500, 0
while done < NPERM:
    B = min(BATCH, NPERM - done)
    rank = np.argsort(np.argsort(rng.random((B, N)), axis=1), axis=1)
    bmax = np.zeros(B)
    bmax_sv = {s: np.zeros(B) for s in CD.sv_locus}
    for gkey, members in groups.items():
        idx_m = tests[members[0]]["idx_m"]
        perm = np.argsort(rank[:, idx_m], axis=1)
        for i in members:
            t = tests[i]
            yp = t["y"][perm]
            yp = (yp - yp.mean(axis=1, keepdims=True)) / yp.std(axis=1, keepdims=True)
            rp = yp @ t["xr"] / t["n"]
            tp = np.abs(rp * np.sqrt(t["df"]) / np.sqrt(np.maximum(1 - rp ** 2, 1e-30)))
            cnt[i] += int((tp >= abs_t[i]).sum())
            bmax = np.maximum(bmax, tp)
            bmax_sv[t["sv"]] = np.maximum(bmax_sv[t["sv"]], tp)
    maxT[done:done + B] = bmax
    for s in CD.sv_locus:
        maxT_sv[s][done:done + B] = bmax_sv[s]
    done += B
    print(f"  perms {done}/{NPERM}")
cnt_cell = cnt

R["p_perm"] = (cnt_cell + 1) / (NPERM + 1)
R["fwer_sv"] = [float((maxT_sv[s] >= t).sum() + 1) / (NPERM + 1)
                for s, t in zip(R.sv_locus, R.abs_t)]
R["fwer_global"] = [float((maxT >= t).sum() + 1) / (NPERM + 1) for t in R.abs_t]
tr = R.relation == "trans"
q = np.full(len(R), np.nan)
p = R.p_perm.values[tr.values]
o = np.argsort(p)
qq = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
q[np.where(tr.values)[0][o]] = np.minimum(qq, 1.0)
R["bh_q_trans"] = q
R = R.sort_values(["relation", "p_perm"])
R.to_csv(f"{OUT}/focused_trans_results.tsv", sep="\t", index=False, float_format="%.6g")
print(f"\nwrote {OUT}/focused_trans_results.tsv")

print("\n=== CIS positive control (must reproduce the known chr8 result) ===")
print(R[R.relation == "cis"].nsmallest(10, "p_perm")[
    ["sv_locus", "pos", "svtype", "phenotype", "target_chrom", "n", "t", "p_perm", "fwer_sv"]
].to_string(index=False))
print("\n=== TRANS, strongest 15 ===")
print(R[tr].nsmallest(15, "p_perm")[
    ["sv_locus", "pos", "svtype", "phenotype", "target_chrom", "n", "t", "p_perm",
     "fwer_sv", "fwer_global", "bh_q_trans"]].to_string(index=False))
n_sig = int((R.loc[tr, "bh_q_trans"] < 0.05).sum())
print(f"\ntrans tests: {int(tr.sum())}   BH q<0.05: {n_sig}   "
      f"global-FWER<0.05: {int((R.loc[tr,'fwer_global'] < 0.05).sum())}")
print(f"permutation floor = {1/(NPERM+1):.2e}")

ESC_P, ESC_N = 0.005, 1_000_000
sel = R.index[(R.relation == "trans") & (R.p_perm <= ESC_P)].tolist()
print(f"\n=== stage 2: escalating {len(sel)} trans cells to {ESC_N:,} permutations ===")
esc = {}
rng2 = np.random.default_rng(SEED + 1)
for k in sel:
    i = int(R.test_i.loc[k])
    t = tests[i]
    obs, cnt2, done2 = float(R.abs_t.loc[k]), 0, 0
    while done2 < ESC_N:
        B = min(20_000, ESC_N - done2)
        yp = t["y"][np.argsort(rng2.random((B, t["n"])), axis=1)]
        yp = (yp - yp.mean(axis=1, keepdims=True)) / yp.std(axis=1, keepdims=True)
        rp = yp @ t["xr"] / t["n"]
        tp = np.abs(rp * np.sqrt(t["df"]) / np.sqrt(np.maximum(1 - rp ** 2, 1e-30)))
        cnt2 += int((tp >= obs).sum()); done2 += B
    esc[k] = (cnt2 + 1) / (ESC_N + 1)
R["p_perm_1e6"] = R.index.map(esc).astype(float)

pv = R.p_perm.copy()
pv[R.p_perm_1e6.notna()] = R.p_perm_1e6[R.p_perm_1e6.notna()]
R["p_best"] = pv
p = R.p_best.values[tr.values]
o = np.argsort(p)
qq = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
q2 = np.full(len(R), np.nan)
q2[np.where(tr.values)[0][o]] = np.minimum(qq, 1.0)
R["bh_q_trans_escalated"] = q2

loo = {}
for k in sel:
    t = tests[int(R.test_i.loc[k])]
    x, y, n = t["xr"], t["y"], t["n"]
    ts = []
    for j in range(n):
        m = np.ones(n, bool); m[j] = False
        xa, ya = x[m], y[m]
        xa = (xa - xa.mean()) / xa.std(); ya = (ya - ya.mean()) / ya.std()
        rr = float(xa @ ya / (n - 1))
        ts.append(rr * np.sqrt(t["df"] - 1) / np.sqrt(max(1 - rr ** 2, 1e-30)))
    ts = np.array(ts)
    loo[k] = (ts.min(), ts.max(), float(np.abs(ts - R.t.loc[k]).max()),
              int((np.sign(ts) != np.sign(R.t.loc[k])).sum()))
R["loo_t_min"] = R.index.map({k: v[0] for k, v in loo.items()})
R["loo_t_max"] = R.index.map({k: v[1] for k, v in loo.items()})
R["loo_max_shift"] = R.index.map({k: v[2] for k, v in loo.items()})
R["loo_sign_flips"] = R.index.map({k: v[3] for k, v in loo.items()})

R.to_csv(f"{OUT}/focused_trans_results.tsv", sep="\t", index=False, float_format="%.6g")
E = R.loc[sel].sort_values("p_best")
print(E[["sv_locus", "phenotype", "target_chrom", "n", "t", "p_perm", "p_perm_1e6",
         "bh_q_trans_escalated", "fwer_sv", "fwer_global",
         "loo_t_min", "loo_t_max", "loo_max_shift"]].to_string(index=False))
print(f"\nAFTER ESCALATION -- trans BH q<0.05: "
      f"{int((R.loc[tr,'bh_q_trans_escalated'] < 0.05).sum())}   "
      f"global-FWER<0.05: {int((R.loc[tr,'fwer_global'] < 0.05).sum())}")
print("\ndone")
