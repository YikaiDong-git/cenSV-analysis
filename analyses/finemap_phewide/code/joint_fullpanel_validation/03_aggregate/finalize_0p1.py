#!/usr/bin/env python3
"""Merged cv2>=0.1 fanout verdict over the FULL 2,978 survivor_gc cells."""
import glob, math, os
from collections import Counter, defaultdict
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
FAN = CENSV_SCRATCH + "/ssv_joint"
FULL = f"{FAN}/manifest_joint_0p1_FULL.tsv"
OUT = f"{FAN}/fanout_out/agg_0p1_FULL"
def ff(x):
    try: return float(x)
    except: return float("nan")

res = {}; both = defaultdict(dict)
for d, tag in [(f"{FAN}/fanout_out/0p3", "0p3"), (f"{FAN}/fanout_out/0p1", "0p1")]:
    for f in glob.glob(d + "/*_fanout_diag.tsv"):
        with open(f) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) < len(hdr): continue
                r = dict(zip(hdr, p))
                k = (r["region"], r["sv_locus"], r["trait"], r["arm"])
                res[k] = r
                both[k][tag] = r.get("sv_pip", "NA")

disagree = 0; overlap = 0
for k, v in both.items():
    if "0p3" in v and "0p1" in v:
        overlap += 1
        if v["0p3"] != v["0p1"]: disagree += 1

ref = {}
with open(FULL) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    order = []
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        key = (p[ix["region"]], p[ix["sv_locus"]], p[ix["trait"]])
        ref[key] = (ff(p[ix["cv"]]), ff(p[ix["pip_matched"]]), ff(p[ix["pip_impute"]]), ff(p[ix["svSNP_r"]]))
        order.append(key)

recs = []; missing = 0
for (region, sv, trait) in order:
    m = res.get((region, sv, trait, "matched")); im = res.get((region, sv, trait, "impute"))
    if m is None and im is None:
        missing += 1; continue
    cv, jmargz, jpip, mxz = ref[(region, sv, trait)]
    pm12 = ff(m.get("sv_pip", "nan")) if m else float("nan")
    pi12 = ff(im.get("sv_pip", "nan")) if im else float("nan")
    best12 = max([x for x in (pm12, pi12) if not math.isnan(x)], default=float("nan"))
    recs.append(dict(region=region, sv=sv, trait=trait, cv=cv, joint_margz=jmargz, joint_pip=jpip,
                     max_abs_z_snp=mxz, pip_m_12k=pm12, pip_i_12k=pi12, best_12k=best12,
                     max_snp_m=ff(m.get("max_snp_pip", "nan")) if m else float("nan"),
                     max_snp_i=ff(im.get("max_snp_pip", "nan")) if im else float("nan"),
                     survivor=(not math.isnan(best12)) and best12 > 0.5,
                     impute_surv=(not math.isnan(pi12)) and pi12 > 0.5,
                     matched_surv=(not math.isnan(pm12)) and pm12 > 0.5))
recs.sort(key=lambda r: (-r["best_12k"] if not math.isnan(r["best_12k"]) else 1, r["region"]))
imp = [r for r in recs if r["impute_surv"]]
cols = ["region","sv","trait","cv","joint_margz","joint_pip","max_abs_z_snp","pip_m_12k","pip_i_12k",
        "best_12k","max_snp_m","max_snp_i","survivor","impute_surv","matched_surv"]
os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/fanout_all_cells_0p1_FULL.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in recs:
        fh.write("\t".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
with open(f"{OUT}/fanout_survivors_0p1_FULL.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in [x for x in recs if x["survivor"]]:
        fh.write("\t".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

print(f"[0p1 FULL] cells resolved: {len(recs)}/{len(order)}  (missing fanout: {missing})")
print(f"[reuse consistency] overlap-in-both-dirs cells(x arm): {overlap}  disagreements: {disagree}")
print(f"  COLLAPSE (both arms lose): {sum(1 for r in recs if not r['survivor'])}")
print(f"  survive MATCHED-only (finite-z artifact): {sum(1 for r in recs if r['matched_surv'] and not r['impute_surv'])}")
print(f"  survive IMPUTE arm (DECISIVE): {len(imp)} cells = {len(set(r['sv'] for r in imp))} SVs")
print(f"  survive EITHER arm: {sum(1 for r in recs if r['survivor'])} cells = {len(set(r['sv'] for r in recs if r['survivor']))} SVs")
print("\nTOP impute-survivors (region sv trait cv joint_pip pip_i_12k maxSNP_i):")
for r in imp[:30]:
    print(f"  {r['region'][:20]:20s} {r['sv'][:30]:30s} {r['trait'][:12]:12s} cv={r['cv']:.2f} jpip={r['joint_pip']:.3f} pi12k={r['pip_i_12k']:.3f} maxSNP={r['max_snp_i']:.2f}")
