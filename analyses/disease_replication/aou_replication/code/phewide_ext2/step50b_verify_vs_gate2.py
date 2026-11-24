#!/usr/bin/env python3
"""Step 50b — GATE on the full rescore: does step50's pgen reproduce Gate 2's dosages?"""
import os, sys, csv, json
import numpy as np
import pgenlib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phewide_ext"))
import step42_gate2_headline_rescore as G
import step51_disease_scoretest_ext2 as S

R_MIN = 0.9999
MAXDELTA = 1e-3
OUT = f"{G.RV}/step28_phewide_ext2"


def main():
    ref_path = f"{G.OUT}/rowlevel/dosages.npz"
    if not os.path.exists(ref_path):
        raise SystemExit(f"ABORT: Gate-2 reference not found at {ref_path}")
    z = np.load(ref_path)
    want = {k.split("::", 1)[1] for k in z.files if k.startswith("new::")}
    print(f"Gate-2 reference: {len(want)} SVs", flush=True)

    pre = S.manifest_prefixes()
    got = {}
    for pfx, nv in pre:
        ids = [l.split("\t")[2] for l in open(f"{pfx}.pvar") if not l.startswith("#")]
        local = {sid: i for i, sid in enumerate(ids) if sid in want}
        if not local:
            continue
        n = len(S.psam_iids(pfx)); buf = np.empty(n, np.float64)
        r = pgenlib.PgenReader(f"{pfx}.pgen".encode())
        for sid, v in local.items():
            r.read_dosages(v, buf); got[sid] = buf.copy()
        r.close()
    missing = sorted(want - set(got))
    if missing:
        raise SystemExit(f"ABORT: {len(missing)} headline SVs absent from the ext2 pgen: {missing}")

    qc = {r["sv_locus"]: r for r in csv.DictReader(open(f"{OUT}/dosage_qc_ext2.tsv"), delimiter="\t")}
    g2 = {r["sv_locus"]: r for r in csv.DictReader(open(f"{G.RESULTS}/gate2_dosage_correlation.tsv"),
                                                   delimiter="\t")}

    rows = []; bad = []
    for sv in sorted(want):
        a = z[f"new::{sv}"].astype(np.float64); b = got[sv]
        if a.shape != b.shape:
            raise SystemExit(f"ABORT: {sv} length {b.shape} vs reference {a.shape}")
        r = G.pearson(a, b)
        md = float(np.max(np.abs(a - b)))
        nu_new = int(qc[sv]["n_used"]); nu_g2 = int(g2[sv]["n_used_new"])
        ok = (np.isfinite(r) and r > R_MIN and md < MAXDELTA and nu_new == nu_g2)
        if not ok:
            bad.append(sv)
        rows.append(dict(sv_locus=sv, r_vs_gate2=r, max_abs_delta=md,
                         n_used_step50=nu_new, n_used_gate2=nu_g2, pass_=int(ok)))
        print(f"  {'OK  ' if ok else 'FAIL'} r={r:.8f}  max|d|={md:.2e}  "
              f"n_used {nu_new} vs {nu_g2}  {sv}", flush=True)

    ks = list(rows[0].keys())
    with open(G.guard(f"{OUT}/step50b_vs_gate2.tsv"), "w") as f:
        f.write("\t".join(ks) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[k]:.10g}" if isinstance(r[k], float) else str(r[k]) for k in ks) + "\n")

    if bad:
        raise SystemExit(f"\nGATE FAILED for {len(bad)} SVs: {bad}\n"
                         f"The full-rescore pipeline does NOT reproduce Gate 2's dosages — stop and fix.")
    print(f"\nGATE PASSED: all {len(rows)} headline SVs reproduce Gate 2 "
          f"(min r {min(r['r_vs_gate2'] for r in rows):.8f}, "
          f"max |delta| {max(r['max_abs_delta'] for r in rows):.2e})", flush=True)
    with open(G.guard(f"{OUT}/step50b_summary.json"), "w") as f:
        json.dump(dict(n=len(rows), min_r=min(r["r_vs_gate2"] for r in rows),
                       max_delta=max(r["max_abs_delta"] for r in rows), all_pass=True), f, indent=1)


if __name__ == "__main__":
    main()
