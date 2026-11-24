#!/usr/bin/env python3
"""Aggregate the joint-arm full-panel fanout diags -> real-vs-artifact table for one cv2 version."""
import glob, math, sys
from collections import Counter, defaultdict
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
VER = sys.argv[1] if len(sys.argv) > 1 else "0p3"
FAN = CENSV_SCRATCH + "/ssv_joint"
DIAG_DIR = f"{FAN}/fanout_out/{VER}"
MANIFEST = f"{FAN}/manifest_joint_{VER}.tsv"
OUTBASE = f"{FAN}/fanout_out/agg_{VER}"

def ff(x):
    try: return float(x)
    except: return float("nan")

ref = {}
with open(MANIFEST) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        ref[(p[ix["region"]], p[ix["sv_locus"]], p[ix["trait"]])] = (
            ff(p[ix["cv"]]), ff(p[ix["pip_matched"]]), ff(p[ix["pip_impute"]]), ff(p[ix["svSNP_r"]]))

cells = defaultdict(dict)
diag_files = sorted(glob.glob(DIAG_DIR + "/*_fanout_diag.tsv"))
for f in diag_files:
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr): continue
            r = dict(zip(hdr, p))
            cells[(r["region"], r["sv_locus"], r["trait"])][r["arm"]] = r
print(f"[{VER}] diag files: {len(diag_files)}  cells: {len(cells)}")

recs = []
for (region, sv, trait), arms in cells.items():
    m = arms.get("matched", {}); im = arms.get("impute", {})
    cv, jmargz, jpip, mxz = ref.get((region, sv, trait), (float("nan"),)*4)
    pm12 = ff(m.get("sv_pip", "nan")); pi12 = ff(im.get("sv_pip", "nan"))
    best12 = max([x for x in (pm12, pi12) if not math.isnan(x)], default=float("nan"))
    survivor = (not math.isnan(best12)) and best12 > 0.5
    impute_surv = (not math.isnan(pi12)) and pi12 > 0.5
    matched_surv = (not math.isnan(pm12)) and pm12 > 0.5
    recs.append(dict(region=region, sv=sv, trait=trait, cv=cv, joint_margz=jmargz, joint_pip=jpip,
                     max_abs_z_snp=mxz, pip_m_12k=pm12, pip_i_12k=pi12, best_12k=best12,
                     max_snp_m=ff(m.get("max_snp_pip", "nan")), max_snp_i=ff(im.get("max_snp_pip", "nan")),
                     conv_m=m.get("converged", "NA"), conv_i=im.get("converged", "NA"),
                     survivor=survivor, impute_surv=impute_surv, matched_surv=matched_surv))

recs.sort(key=lambda r: (-r["best_12k"] if not math.isnan(r["best_12k"]) else 1, r["region"]))
surv = [r for r in recs if r["survivor"]]
imp = [r for r in recs if r["impute_surv"]]

cols = ["region","sv","trait","cv","joint_margz","joint_pip","max_abs_z_snp","pip_m_12k","pip_i_12k",
        "best_12k","max_snp_m","max_snp_i","conv_m","conv_i","survivor","impute_surv","matched_surv"]
def w(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
import os; os.makedirs(OUTBASE, exist_ok=True)
w(f"{OUTBASE}/fanout_all_cells_{VER}.tsv", recs)
w(f"{OUTBASE}/fanout_survivors_{VER}.tsv", [r for r in recs if r["survivor"]])

n_matched_only = sum(1 for r in recs if r["matched_surv"] and not r["impute_surv"])
print(f"\n=== [{VER}] VERDICT: {len(recs)} cells ===")
print(f"  COLLAPSE (both arms lose): {sum(1 for r in recs if not r['survivor'])}")
print(f"  survive MATCHED-only (finite-z artifact): {n_matched_only}")
print(f"  survive IMPUTE arm (DECISIVE / candidate-real): {len(imp)} cells = {len(set(r['sv'] for r in imp))} SVs")
print(f"  survive EITHER arm: {len(surv)} cells = {len(set(r['sv'] for r in surv))} SVs")
print("\nimpute-survivors by cv2 bin:")
ab = Counter(min(int(r["cv"]*10),9) for r in recs if not math.isnan(r["cv"]))
sb = Counter(min(int(r["cv"]*10),9) for r in imp if not math.isnan(r["cv"]))
for k in sorted(ab): print(f"  cv~{k/10:.1f}: impute-survive {sb.get(k,0)}/{ab[k]}")
print("\nTOP impute-survivors (region sv trait cv joint_pip pip_i_12k max_snp_i):")
for r in imp[:30]:
    print(f"  {r['region'][:20]:20s} {r['sv'][:32]:32s} {r['trait'][:12]:12s} cv={r['cv']:.2f} jpip={r['joint_pip']:.3f} pi12k={r['pip_i_12k']:.3f} maxSNP={r['max_snp_i']:.2f}")
