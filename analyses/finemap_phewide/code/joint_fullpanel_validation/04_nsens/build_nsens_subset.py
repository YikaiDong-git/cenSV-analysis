#!/usr/bin/env python3
import math
import os as _censv_os
CENSV_SCRATCH = _censv_os.environ["CENSV_SCRATCH"]
FAN = CENSV_SCRATCH + "/ssv_joint"
ALL = f"{FAN}/fanout_out/agg_0p3/fanout_all_cells_0p3.tsv"
MAN = f"{FAN}/manifest_joint_0p3.tsv"
K = 50

def ff(x):
    try: return float(x)
    except: return float("nan")

classes = {"survive": [], "snpwin": [], "degen": [], "nonconv": []}
with open(ALL) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(h)}
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        key = (p[ix["region"]], p[ix["sv"]], p[ix["trait"]])
        pi = ff(p[ix["pip_i_12k"]]); mx = ff(p[ix["max_snp_i"]]); conv = p[ix["conv_i"]]
        if p[ix["impute_surv"]] == "True": classes["survive"].append(key)
        elif conv != "TRUE": classes["nonconv"].append(key)
        elif mx > 0.5: classes["snpwin"].append(key)
        else: classes["degen"].append(key)

pick = set(classes["survive"])
for c in ("snpwin", "degen", "nonconv"):
    pick |= set(sorted(classes[c])[:K])

out = []
with open(MAN) as fh:
    hdr = fh.readline()
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if (p[0], p[4], p[6]) in pick:
            out.append(ln.rstrip("\n"))
out.sort()
with open(f"{FAN}/manifest_nsens.tsv", "w") as fh:
    fh.write(hdr); [fh.write(l + "\n") for l in out]
print(f"class sizes: " + ", ".join(f"{c}={len(v)}" for c, v in classes.items()))
print(f"subset: {len(out)} cells  (all {len(classes['survive'])} survivors + up to {K} each of snpwin/degen/nonconv)")
