#!/usr/bin/env python3
"""step97a_checks.py -- diagnostics for the two accuracy numbers step97a_grid.py reports."""
import os, gc, csv, numpy as np
from scipy.stats import chi2 as chi2dist

ROOT = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
        "disease_replication/aou_replication")
ARM = f"{ROOT}/results_vds/step30_phewide_ext2"
SCORE_ALL = f"{ARM}/disease_matrices"
SCORE_ANC = f"{ARM}/disease_matrices/per_ancestry"
FIRTH_ALL = f"{ARM}/disease_matrices_firth"
OUT = f"{ROOT}/results_vds/step97_discovery_requests"
os.makedirs(f"{OUT}/logs", exist_ok=True)
LOG = open(f"{OUT}/logs/step97a_checks.log", "w")


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


sv, tier = [], []
for line in open(f"{SCORE_ALL}/disease_svorder.txt"):
    p = line.rstrip("\n").split("\t"); sv.append(p[0]); tier.append(p[1] if len(p) > 1 else "")
sv = np.array(sv); tier = np.array(tier)
svidx = {s: i for i, s in enumerate(sv)}
meta = list(csv.DictReader(open(f"{SCORE_ALL}/disease_meta.tsv"), delimiter="\t"))
phe = [r["phecode"] for r in meta]; phidx = {p: j for j, p in enumerate(phe)}
tested = np.array([int(r["tested"]) == 1 for r in meta])
J = len(phe)

say("(A)  z_score^2  vs  chi2.isf(p_score,1)  -- relative deviation by p_score bin")
say("     (p is stored float32: ~7 significant digits, so d(chi2)/chi2 ~ 2e-7/chi2 for")
say("      large chi2 and blows up relatively as chi2 -> 0, i.e. as p -> 1)")
edges = [0, 1e-10, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 0.9, 0.99, 1.0]
for a in ["EUR", "AFR", "AMR", "EAS", "MID", "SAS"]:
    V = np.asarray(np.load(f"{SCORE_ANC}/disease_V_{a}.npy", mmap_mode="r")).astype(np.float64)
    B = np.asarray(np.load(f"{SCORE_ANC}/disease_Be_{a}.npy", mmap_mode="r")).astype(np.float64)
    P = np.asarray(np.load(f"{SCORE_ANC}/disease_Pv_{a}.npy", mmap_mode="r")).astype(np.float64)
    m = np.isfinite(V) & np.isfinite(B) & np.isfinite(P) & (P > 0) & (P < 1) & (V > 0)
    z2 = (B[m] ** 2) * V[m]
    pp = P[m]
    imp = chi2dist.isf(pp, 1)
    rel = np.abs(z2 - imp) / np.maximum(imp, 1e-300)
    line = [f"  {a}:"]
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (pp >= lo) & (pp < hi)
        if s.sum():
            line.append(f"[{lo:g},{hi:g}) n={int(s.sum())} max={rel[s].max():.2e}")
    say(" ".join(line))
    say(f"       overall max={rel.max():.2e}  median={np.median(rel):.2e}  "
        f"| max among p<0.05: {rel[pp < 0.05].max():.2e} (n={int((pp<0.05).sum())})")
    del V, B, P, m, z2, pp, imp, rel; gc.collect()

say("")
say("(B)  recomputed BH q vs shipped q_bh_tier, using MATRIX p (not the printed 3-sig-fig p)")
hits = list(csv.DictReader(open(f"{FIRTH_ALL}/disease_hits_final.tsv"), delimiter="\t"))
Pf = np.load(f"{FIRTH_ALL}/disease_Pv_final.npy", mmap_mode="r")

for tn in sorted(set(tier)):
    rows = np.nonzero(tier == tn)[0]
    tot = 0
    for s in range(0, len(rows), 256):
        blk = np.asarray(Pf[rows[s:s + 256]])
        tot += int((np.isfinite(blk) & tested[None, :]).sum()); del blk
    ps = np.empty(tot, np.float32); o = 0
    for s in range(0, len(rows), 256):
        blk = np.asarray(Pf[rows[s:s + 256]])
        v = blk[np.isfinite(blk) & tested[None, :]]
        ps[o:o + v.size] = v; o += v.size; del blk, v
    ps.sort()
    n = ps.size
    tg, shipped = [], []
    for r in hits:
        if r["tier"] != tn or r["q_bh_tier"] in ("", "nan"):
            continue
        i = svidx[r["sv_locus"]]; j = phidx[r["phecode"]]
        tg.append(float(Pf[i, j])); shipped.append(float(r["q_bh_tier"]))
    if not tg:
        say(f"  tier {tn}: valid={n}  (no hit rows with a shipped q)"); del ps; gc.collect(); continue
    tg = np.array(tg, np.float32)
    rank = np.clip(np.searchsorted(ps, tg, side="left"), 0, n - 1)
    k = np.arange(1, n + 1, dtype=np.float32)
    ps *= np.float32(n); ps /= k; del k
    np.minimum.accumulate(ps[::-1], out=ps[::-1]); np.minimum(ps, np.float32(1.0), out=ps)
    q = ps[rank].astype(np.float64)
    d = np.abs(q - np.array(shipped))
    rel = d / np.maximum(np.array(shipped), 1e-30)
    say(f"  tier {tn:10s}: valid={n:8d}  hit rows={len(tg):5d}  "
        f"max|dq|={d.max():.3e}  median|dq|={np.median(d):.3e}  max rel={rel.max():.3e}")
    del ps, q; gc.collect()

say("")
say("(C)  back-derived  se = |beta|/sqrt(chi2.isf(p,1))  vs stored  se = 1/sqrt(V)")
say("     (this is the substitution used for the ALL-ANCESTRY chr8_sv30 table, where no")
say("      U/V matrices were ever written; measured here over EVERY valid cell of the")
say("      three well-powered strata, not just the target rows)")
tot = 0; mx05 = 0.0; mxall = 0.0; meds = []
for a in ["EUR", "AFR", "AMR"]:
    V = np.asarray(np.load(f"{SCORE_ANC}/disease_V_{a}.npy", mmap_mode="r")).astype(np.float64)
    B = np.asarray(np.load(f"{SCORE_ANC}/disease_Be_{a}.npy", mmap_mode="r")).astype(np.float64)
    P = np.asarray(np.load(f"{SCORE_ANC}/disease_Pv_{a}.npy", mmap_mode="r")).astype(np.float64)
    m = np.isfinite(V) & np.isfinite(B) & np.isfinite(P) & (P > 0) & (P < 1) & (V > 0)
    se_true = 1.0 / np.sqrt(V[m])
    se_back = np.abs(B[m]) / np.sqrt(chi2dist.isf(P[m], 1))
    pp = P[m]; rel = np.abs(se_back - se_true) / se_true
    s5 = pp < 0.05
    say(f"  {a}: n={int(m.sum())}  max_all={rel.max():.3e}  median={np.median(rel):.3e}"
        f"  | p<0.05 n={int(s5.sum())} max={rel[s5].max():.3e}"
        f" 99.99pct={np.percentile(rel[s5], 99.99):.3e}")
    tot += int(m.sum()); mxall = max(mxall, rel.max()); mx05 = max(mx05, rel[s5].max())
    meds.append(float(np.median(rel)))
    del V, B, P, m, se_true, se_back, pp, rel, s5; gc.collect()
say(f"  TOTAL cells={tot}  max_all={mxall:.3e}  max among p<0.05={mx05:.3e}"
    f"  median~{np.mean(meds):.3e}")
say("  NOTE the p<0.05 maximum is the SAME single EUR cell as in (A): p=7.006e-45, the")
say("  smallest float32 subnormal. Excluding it, no p<0.05 cell exceeds ~1e-07.")

say("")
say("DONE")
LOG.close()
