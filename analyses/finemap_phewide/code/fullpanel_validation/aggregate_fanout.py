#!/usr/bin/env python3
"""Aggregate the per-region fan-out diags into the definitive real-vs-artifact table."""
import glob, math
from collections import Counter, defaultdict
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
FAN = CENSV_SCRATCH + "/ssv_fanout"

def ff(x):
    try: return float(x)
    except: return float("nan")

cells = defaultdict(dict)
diag_files = sorted(glob.glob(FAN + "/out/*_fanout_diag.tsv"))
for f in diag_files:
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            r = dict(zip(hdr, p))
            cells[(r["sv_locus"], r["trait"])][r["arm"]] = r

print(f"regions aggregated: {len(diag_files)}  cells: {len(cells)}")

recs = []
for (sv, trait), arms in cells.items():
    m = arms.get("matched", {}); im = arms.get("impute", {})
    any_r = m or im
    region = any_r.get("region", ""); cv = ff(any_r.get("cv", "nan")); rr = ff(any_r.get("svSNP_r_union", "nan"))
    pm_u = ff(any_r.get("pip_matched_union", "nan")); pi_u = ff(any_r.get("pip_impute_union", "nan"))
    pm12 = ff(m.get("sv_pip", "nan")); pi12 = ff(im.get("sv_pip", "nan"))
    best12 = max([x for x in (pm12, pi12) if not math.isnan(x)], default=float("nan"))
    survivor = (not math.isnan(best12)) and best12 > 0.5
    recs.append(dict(region=region, sv=sv, trait=trait, cv=cv, svSNP_r=rr,
                     pip_m_union=pm_u, pip_i_union=pi_u, pip_m_12k=pm12, pip_i_12k=pi12,
                     best_12k=best12, survivor=survivor,
                     max_snp_m=ff(m.get("max_snp_pip", "nan")), max_snp_i=ff(im.get("max_snp_pip", "nan"))))

recs.sort(key=lambda r: (-r["best_12k"] if not math.isnan(r["best_12k"]) else 1, r["region"]))
surv = [r for r in recs if r["survivor"]]
coll = [r for r in recs if not r["survivor"]]

def w(path, rows):
    cols = ["region", "sv", "trait", "cv", "svSNP_r", "pip_m_union", "pip_i_union", "pip_m_12k", "pip_i_12k", "best_12k", "max_snp_m", "max_snp_i", "survivor"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
w(FAN + "/fanout_all_cells.tsv", recs)
w(FAN + "/fanout_survivors.tsv", surv)

print(f"\n=== VERDICT: {len(surv)} SURVIVE (12k sv_pip>0.5 either arm) / {len(recs)} cells ; {len(coll)} COLLAPSE ===")
print(f"unique survivor SVs: {len(set(r['sv'] for r in surv))} / {len(set(r['sv'] for r in recs))}")
print("\nsurvivors by cv2 bin:")
sb = Counter(min(int(r["cv"] * 10), 9) for r in surv if not math.isnan(r["cv"]))
ab = Counter(min(int(r["cv"] * 10), 9) for r in recs if not math.isnan(r["cv"]))
for k in sorted(ab):
    print(f"  cv~{k/10:.1f}: survive {sb.get(k,0)}/{ab[k]}")
print("\nsurvivor arm: impute-only", sum(1 for r in surv if ff(str(r['pip_i_12k']))>0.5 and not(ff(str(r['pip_m_12k']))>0.5)),
      "| matched-only", sum(1 for r in surv if ff(str(r['pip_m_12k']))>0.5 and not(ff(str(r['pip_i_12k']))>0.5)),
      "| both", sum(1 for r in surv if ff(str(r['pip_m_12k']))>0.5 and ff(str(r['pip_i_12k']))>0.5))
print("\nTOP survivors (region sv trait cv svSNP_r best_12k):")
for r in surv[:25]:
    print(f"  {r['region'][:22]:22s} {r['sv'][:34]:34s} {r['trait'][:12]:12s} cv={r['cv']:.2f} r={r['svSNP_r']:.2f} pip12k={r['best_12k']:.3f}")
