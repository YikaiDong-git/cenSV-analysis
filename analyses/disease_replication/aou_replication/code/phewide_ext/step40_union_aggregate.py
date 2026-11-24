#!/usr/bin/env python3
import os, sys, csv, argparse
import numpy as np
from collections import Counter


def bh(pvals):
    ps = np.asarray(pvals, float); ok = np.isfinite(ps); n = int(ok.sum())
    q = np.full(len(ps), np.nan)
    if n == 0:
        return q
    ok_idx = np.nonzero(ok)[0]
    order = ok_idx[np.argsort(ps[ok_idx])]
    prev = 1.0
    for k, i in enumerate(order[::-1]):
        prev = min(prev, ps[i] * n / (n - k)); q[i] = prev
    return q


def load_arm(name, root, tier_override=None):
    mat = f"{root}/disease_matrices"; fir = f"{root}/disease_matrices_firth"
    pv_p = f"{fir}/disease_Pv_final.npy"; be_p = f"{fir}/disease_Be_final.npy"
    if not os.path.exists(pv_p):
        pv_p = f"{mat}/disease_Pv.npy"; be_p = f"{mat}/disease_Be.npy"
        print(f"  [{name}] NOTE: no Firth-final matrix; using SCORE matrices {os.path.basename(pv_p)}", flush=True)
    Pv = np.load(pv_p).astype(np.float64); Be = np.load(be_p).astype(np.float64)
    svrows = [ln.rstrip("\n").split("\t") for ln in open(f"{mat}/disease_svorder.txt") if ln.strip()]
    sv = [r[0] for r in svrows]
    tier = [(tier_override or (r[1] if len(r) > 1 else "")) for r in svrows]
    meta = list(csv.DictReader(open(f"{mat}/disease_meta.tsv"), delimiter="\t"))
    sep = {}
    cells = f"{fir}/disease_firth_cells.tsv"
    if os.path.exists(cells):
        for r in csv.DictReader(open(cells), delimiter="\t"):
            try:
                sep[(r["sv_locus"], r["phecode"])] = (int(r.get("separation", 0) or 0),
                                                      float(r["fisher_p"]) if r.get("fisher_p") not in ("", "nan", None) else np.nan)
            except (ValueError, KeyError):
                continue
    src = {}
    hits = f"{fir}/disease_hits_final.tsv"
    if os.path.exists(hits):
        for r in csv.DictReader(open(hits), delimiter="\t"):
            src[(r["sv_locus"], r["phecode"])] = r.get("p_source", "score")
    assert Pv.shape == (len(sv), len(meta)), f"{name}: matrix {Pv.shape} != ({len(sv)},{len(meta)})"
    print(f"  [{name}] {Pv.shape[0]:,} SVs x {Pv.shape[1]:,} phecodes | tiers {dict(Counter(tier))} "
          f"| sep-cells {len(sep):,}", flush=True)
    return dict(name=name, Pv=Pv, Be=Be, sv=sv, tier=tier, meta=meta, sep=sep, src=src)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", required=True, help="name:results_root[:tier_override]")
    ap.add_argument("--out", required=True)
    ap.add_argument("--validate", default=None, help="expected hits tsv (single-arm regression check)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print("=== loading arms (per-cell p/beta ONLY) ===", flush=True)
    arms = []
    for spec in a.arms:
        p = spec.split(":")
        arms.append(load_arm(p[0], p[1], p[2] if len(p) > 2 else None))

    ref = arms[0]
    ref_ph = [m["phecode"] for m in ref["meta"]]
    for arm in arms[1:]:
        ph = [m["phecode"] for m in arm["meta"]]
        assert ph == ref_ph, (f"ABORT {arm['name']}: phecode columns differ from {ref['name']} "
                              f"({len(ph)} vs {len(ref_ph)}) -- arms are not aligned; union would be meaningless")
        t_ref = [m["tested"] == "1" for m in ref["meta"]]; t_a = [m["tested"] == "1" for m in arm["meta"]]
        assert t_ref == t_a, f"ABORT {arm['name']}: `tested` flags differ from {ref['name']}"
    print(f"GATE PASS: all {len(arms)} arms share identical phecode axis ({len(ref_ph):,}) + tested flags", flush=True)

    Pv = np.concatenate([x["Pv"] for x in arms], axis=0)
    Be = np.concatenate([x["Be"] for x in arms], axis=0)
    sv = [s for x in arms for s in x["sv"]]
    tiers = np.array([t for x in arms for t in x["tier"]])
    sep = {}; src = {}
    for x in arms:
        sep.update(x["sep"]); src.update(x["src"])
    M, J = Pv.shape
    assert len(sv) == M and len(set(sv)) == M, "ABORT: duplicate sv_locus across arms -- would double-count"
    meta = ref["meta"]
    tested = np.array([m["tested"] == "1" for m in meta])
    ph_label = {m["phecode"]: m.get("label", "") for m in meta}
    ph_cases = {m["phecode"]: int(m["cases"]) for m in meta}
    ph_n = {m["phecode"]: int(m["n"]) for m in meta}
    print(f"UNION: {M:,} SVs x {J:,} phecodes | tiers {dict(Counter(tiers.tolist()))}", flush=True)

    sv_of = np.repeat(np.arange(M), J); ph_of = np.tile(np.arange(J), M)
    flat = Pv.ravel(); Bflat = Be.ravel()
    valid = np.isfinite(flat) & tested[ph_of]
    q = np.full(M * J, np.nan)
    for tn in sorted(set(tiers.tolist())):
        if not tn:
            continue
        ix = np.nonzero(valid & (tiers[sv_of] == tn))[0]
        if ix.size:
            q[ix] = bh(flat[ix])
            print(f"  BH pool '{tn}': {ix.size:,} valid cells", flush=True)
    ntest = int(valid.sum())
    with np.errstate(all="ignore"):
        from scipy.stats import chi2 as c2
        lam = {tn: float(c2.isf(np.median(flat[valid & (tiers[sv_of] == tn)]), 1) / 0.4549)
               for tn in sorted(set(tiers.tolist())) if tn and (valid & (tiers[sv_of] == tn)).any()}
    print(f"  recomputed: ntest={ntest:,} (Bonferroni denom)  lambda_GC={ {k: round(v,3) for k,v in lam.items()} }",
          flush=True)

    keep = np.nonzero(valid & ((np.nan_to_num(q, nan=1.0) < 0.10) | (flat < 1e-4)))[0]
    keep = keep[np.argsort(flat[keep])]
    n_sig = n_ind = n_sep = 0
    out = f"{a.out}/disease_hits_union.tsv"
    with open(out, "w") as f:
        f.write("sv_locus\ttier\tphecode\tlabel\tcases\tn\tbeta_final_logOR\tp_final\tp_source\t"
                "q_bh_tier\tseparation\tfisher_p\ttier_label\n")
        for idx in keep:
            i = idx // J; j = idx % J; qi = q[idx]; ph = meta[j]["phecode"]
            s, fp = sep.get((sv[i], ph), (0, np.nan))
            sig = bool(np.isfinite(qi) and qi < 0.05)
            base = "significant" if sig else ("indicative" if (np.isfinite(qi) and qi < 0.10) else "sub-threshold")
            lab = "indicative (separation)" if s else base
            if s and (sig or base == "indicative"):
                n_sep += 1
            elif sig:
                n_sig += 1
            elif base == "indicative":
                n_ind += 1
            f.write(f"{sv[i]}\t{tiers[i]}\t{ph}\t{ph_label.get(ph,'')}\t{ph_cases.get(ph,0)}\t{ph_n.get(ph,0)}\t"
                    f"{Bflat[idx]:.4g}\t{flat[idx]:.3g}\t{src.get((sv[i],ph),'score')}\t"
                    f"{'' if not np.isfinite(qi) else f'{qi:.4g}'}\t{s}\t"
                    f"{'' if not np.isfinite(fp) else f'{fp:.3g}'}\t{lab}\n")
    np.save(f"{a.out}/disease_Pv_union.npy", Pv.astype(np.float32))
    np.save(f"{a.out}/disease_Be_union.npy", Be.astype(np.float32))
    with open(f"{a.out}/disease_svorder_union.txt", "w") as f:
        f.write("\n".join(f"{sv[i]}\t{tiers[i]}" for i in range(M)) + "\n")
    print(f"UNION_DONE: {n_sig} significant / {n_ind} indicative / {n_sep} separation-flagged "
          f"({len(keep):,} hit-band rows) -> {out}", flush=True)

    if a.validate:
        exp = list(csv.DictReader(open(a.validate), delimiter="\t"))
        e = Counter(r["tier_label"] for r in exp)
        got = Counter(r["tier_label"] for r in csv.DictReader(open(out), delimiter="\t"))
        print("\n=== VALIDATION vs " + os.path.basename(a.validate) + " ===", flush=True)
        gotrows = list(csv.DictReader(open(out), delimiter="\t"))
        for k in sorted(set(e) | set(got)):
            m = "OK " if e.get(k, 0) == got.get(k, 0) else "DIFF"
            print(f"  {m}  {k:26} expected {e.get(k,0):5}  got {got.get(k,0):5}")
        ok = True
        for lab in ("significant", "indicative"):
            es = {(r["sv_locus"], r["phecode"]) for r in exp if r["tier_label"] == lab}
            gs = {(r["sv_locus"], r["phecode"]) for r in gotrows if r["tier_label"] == lab}
            same = (es == gs); ok &= same
            print(f"  {'OK ' if same else 'FAIL'}  {lab.upper()} CELL SET identical: {same} "
                  f"(expected {len(es)}, got {len(gs)}, symmetric diff {len(es ^ gs)})")
        exp_cell = {(r["sv_locus"], r["phecode"]): r for r in exp}
        got_cell = {(r["sv_locus"], r["phecode"]): r for r in gotrows}
        extra = set(got_cell) - set(exp_cell); missing = set(exp_cell) - set(got_cell)
        edge = [k for k in extra if abs(float(got_cell[k]["p_final"]) - 1e-4) / 1e-4 < 1e-2]
        nonhit_edge = all(got_cell[k]["tier_label"] == "sub-threshold" for k in extra)
        print(f"  INFO  reporting-band delta: +{len(extra)} / -{len(missing)} rows; {len(edge)}/{len(extra)} "
              f"sit on the p=1e-4 float32 edge; all extras are sub-threshold (non-hit): {nonhit_edge}")
        if extra and not (len(edge) == len(extra) and nonhit_edge and not missing):
            ok = False
            print("  FAIL  band delta is NOT explained by the float32 edge -- investigate")
        print(f"\n[VALIDATION] {'PASS -- all FDR calls reproduced exactly (band delta = float32 edge only)' if ok else 'FAIL'}",
              flush=True)
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
