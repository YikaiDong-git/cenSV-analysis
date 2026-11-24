#!/usr/bin/env python3
"""Compare with-n+check_prior=T survival vs the no-n classification, per no-n class."""
import glob, math
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
FAN = CENSV_SCRATCH + "/ssv_joint"
def ff(x):
    try: return float(x)
    except: return float("nan")

noN = {}
with open(f"{FAN}/fanout_out/agg_0p3/fanout_all_cells_0p3.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        key = (p[ix["region"]], p[ix["sv"]], p[ix["trait"]])
        pi = ff(p[ix["pip_i_12k"]]); mx = ff(p[ix["max_snp_i"]]); conv = p[ix["conv_i"]]
        if p[ix["impute_surv"]] == "True": c = "survive"
        elif conv != "TRUE": c = "nonconv"
        elif mx > 0.5: c = "snpwin"
        else: c = "degen"
        noN[key] = (c, pi)

withN = {}
for f in glob.glob(f"{FAN}/nsens_out/*_nsens_diag.tsv"):
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); jx = {c: i for i, c in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr) or p[jx["arm"]] != "impute": continue
            key = (p[jx["region"]], p[jx["sv_locus"]], p[jx["trait"]])
            withN[key] = (ff(p[jx["sv_pip"]]), p[jx["converged"]], ff(p[jx["n_gwas"]]))

from collections import Counter, defaultdict
tab = defaultdict(lambda: [0, 0, 0])
flips = []
for key, (c, pi) in noN.items():
    if key not in withN: continue
    wp, wconv, ng = withN[key]
    tab[c][0] += 1
    surv = (not math.isnan(wp)) and wp > 0.5
    if surv: tab[c][1] += 1
    if wconv == "TRUE": tab[c][2] += 1
    if c != "survive" and surv:
        flips.append((key, c, pi, wp))

print("with-n (n=GWAS sample size, check_prior=TRUE) vs no-n class:")
print(f"{'no-n class':12s} {'n':>4s} {'with-n survive':>15s} {'with-n converged':>17s}")
for c in ("survive", "snpwin", "degen", "nonconv"):
    n, s, cv = tab[c]
    if n: print(f"{c:12s} {n:>4d} {s:>15d} {cv:>17d}")
print(f"\nCOLLAPSE->SURVIVE flips with n: {len(flips)} (of {sum(tab[c][0] for c in ('snpwin','degen','nonconv'))} collapse cells tested)")
for key, c, pi, wp in flips[:15]:
    print(f"  FLIP {key[1][:28]:28s} {key[2][:10]:10s} class={c} no-n_pip={pi:.3f} with-n_pip={wp:.3f}")
