#!/usr/bin/env python3
"""Step 42c — EXACT tiered-BH recount + per-SD effect sizes for the Gate-2 refits."""
import os, sys, csv, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step42_gate2_headline_rescore as G

RV = G.RV
MAT = f"{RV}/step30_phewide/disease_matrices"
MATF = f"{RV}/step30_phewide/disease_matrices_firth"
OUT = G.OUT
ALPHA = 0.05


def bh_reject(p, alpha=ALPHA):
    ps = np.sort(p[np.isfinite(p)])
    m = ps.size
    if m == 0:
        return 0, np.nan, 0
    k = np.arange(1, m + 1)
    rej = np.where(ps <= k / m * alpha)[0]
    if not rej.size:
        return 0, np.nan, m
    kmax = int(rej[-1]) + 1
    return kmax, float(ps[kmax - 1]), m


def main():
    pre = {r: G.tree_manifest(r) for r in G.FROZEN}

    Pv = np.load(f"{MATF}/disease_Pv_final.npy")
    rows = [l.rstrip("\n").split("\t") for l in open(f"{MAT}/disease_svorder.txt") if l.strip()]
    svo = [r[0] for r in rows]
    tier = [r[1] if len(r) > 1 else "" for r in rows]
    meta = list(csv.DictReader(open(f"{MAT}/disease_meta.tsv"), delimiter="\t"))
    phe = [r["phecode"] for r in meta]
    sig = np.array([t == "signal" for t in tier])
    print(f"Pv {Pv.shape}; signal tier {sig.sum()} SVs x {len(phe)} phecodes", flush=True)

    P = Pv[sig, :].astype(np.float64)
    n0, cut0, m0 = bh_reject(P)
    print(f"  self-check (unmodified): m={m0:,}  rejects={n0}  cutoff={cut0:.6g}", flush=True)
    if n0 != 28:
        raise SystemExit(f"ABORT: BH self-check rejects {n0}, expected the published 28 — "
                         f"the recount cannot be trusted.")

    cells = list(csv.DictReader(open(f"{G.RESULTS}/gate2_cell_refit.tsv"), delimiter="\t"))
    dos = {r["sv_locus"]: r for r in csv.DictReader(open(f"{G.RESULTS}/gate2_dosage_correlation.tsv"),
                                                    delimiter="\t")}
    svi = {s: i for i, s in enumerate(np.array(svo)[sig])}
    phi = {p: j for j, p in enumerate(phe)}

    P2 = P.copy(); nsub = 0
    for r in cells:
        i, j = svi.get(r["sv_locus"]), phi.get(r["phecode"])
        if i is None or j is None:
            raise SystemExit(f"ABORT: cell not on the frozen axes: {r['sv_locus']} {r['phecode']}")
        P2[i, j] = float(r["p_new_likeforlike"]); nsub += 1
    n1, cut1, m1 = bh_reject(P2)
    print(f"  substituted {nsub}/28 cells", flush=True)
    print(f"  RECOUNTED:               m={m1:,}  rejects={n1}  cutoff={cut1:.6g}", flush=True)

    out = []
    for r in cells:
        sv = r["sv_locus"]; d = dos[sv]
        sd_o = float(d["sd_archived"]); sd_n = float(d["sd_new"])
        bo = float(r["beta_firth_old"]); bn = float(r["beta_firth_new"])
        pn = float(r["p_new_likeforlike"]); po = float(r["p_old_likeforlike"])
        out.append(dict(
            sv_locus=sv, phecode=r["phecode"], label=r["label"], cases=r["cases"],
            p_published=float(r["p_published"]), p_old=po, p_new=pn,
            survives_exact_bh=int(pn <= cut1),
            bracket_band=r["cutoff_band"],
            beta_old=bo, beta_new=bn,
            sd_dosage_old=sd_o, sd_dosage_new=sd_n,
            beta_perSD_old=bo * sd_o, beta_perSD_new=bn * sd_n,
            perSD_ratio=(bn * sd_n) / (bo * sd_o) if bo * sd_o != 0 else np.nan,
            r_dosage=float(r["r_dosage"])))
    ks = list(out[0].keys())
    with open(G.guard(f"{G.RESULTS}/gate2_bh_recount_persd.tsv"), "w") as f:
        f.write("\t".join(ks) + "\n")
        for r in out:
            f.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in ks) + "\n")

    nsurv = sum(r["survives_exact_bh"] for r in out)
    print(f"\n  EXACT BH: {nsurv}/28 of the published cells survive (bracket said "
          f"{sum(1 for r in out if r['bracket_band']=='KEEP')})", flush=True)
    print("  --- per-SD effect change (scale-free; raw beta ratios are meaningless after the [0,2] "
          "rescale) ---", flush=True)
    for r in sorted(out, key=lambda r: r["p_new"]):
        print(f"   {'KEEP' if r['survives_exact_bh'] else 'LOST'}  p {r['p_old']:.3g}->{r['p_new']:.3g}  "
              f"beta/SD {r['beta_perSD_old']:+.4f}->{r['beta_perSD_new']:+.4f} "
              f"({r['perSD_ratio']:+.2f}x)  {r['sv_locus'][:34]:34s} {r['phecode']}", flush=True)

    per_sv = {}
    for r in out:
        a = per_sv.setdefault(r["sv_locus"], [0, 0, r["r_dosage"]])
        a[0] += r["survives_exact_bh"]; a[1] += 1
    print("\n  per SV (exact BH):", flush=True)
    for sv, (k, n, rr) in sorted(per_sv.items(), key=lambda kv: -kv[1][0] / kv[1][1]):
        print(f"    {k}/{n}  r={rr:.4f}  {sv}", flush=True)

    with open(G.guard(f"{G.RESULTS}/gate2_bh_recount_summary.json"), "w") as f:
        json.dump(dict(m_signal_tier=m1, published_rejects=n0, published_cutoff=cut0,
                       recounted_rejects=n1, recounted_cutoff=cut1,
                       of_28_surviving_exact=nsurv,
                       of_28_surviving_bracket=sum(1 for r in out if r["bracket_band"] == "KEEP"),
                       caveat="all other ~8.26M cells held at their OLD p; a real rescore perturbs all 6,513"),
                  f, indent=1)

    ok = all(G.tree_manifest(r) == pre[r] for r in G.FROZEN)
    print(f"\n  frozen trees unchanged: {ok}", flush=True)
    if not ok:
        raise SystemExit("ABORT: a frozen tree was modified")


if __name__ == "__main__":
    main()
