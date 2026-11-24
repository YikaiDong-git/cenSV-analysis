#!/usr/bin/env python3
import pandas as pd, numpy as np, re, os, glob
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]
ROOT = CENSV_ROOT + ""
ADMIT = f"{ROOT}/analyses/sv_catalog/output/ref_ukbb/admit/admit_decompose_cv01.tsv"
DBASE = f"{ROOT}/sv_extraction/decompose/mapq20_minsvlen10_maxpart10_partdist1000_maxedit200_maxsv500000_qgap200_rgap200"
PERCH = f"{ROOT}/analyses/cdr_qtl/results/cdr_qtl_corrected_1e6/perchrom"
EN    = f"{ROOT}/results/elastic_net/decompose/cosine_nj"

def parse_ids(s):
    if not isinstance(s, str): return set()
    return set(x for x in re.split(r"[;,]", s.replace("_B", ";_B")) if x.startswith("svim_asm."))

adf = pd.read_csv(ADMIT, sep="\t", dtype=str)
adf = adf[adf["status"] == "ADMIT"].copy()
adf["cv_r2"] = adf["cv_r2"].astype(float)
adf["chr"] = adf["sv_locus"].str.extract(r"^(chr[0-9XY]+)_")
inv_by_chr = {}
for ch, g in adf.groupby("chr"):
    inv = {}
    for _, r in g.iterrows():
        for i in parse_ids(r["sv_id"]):
            inv.setdefault(i, []).append((r["cv_r2"], r["sv_locus"]))
    inv_by_chr[ch] = inv
print(f"[admit] {len(adf)} ADMIT rows across {adf['chr'].nunique()} chroms")
print(f"[admit] chr10 ADMIT loci: {(adf['chr']=='chr10').sum()}")
print(f"[base decompose] chr10 present: {os.path.isdir(f'{DBASE}/chr10')} ; chr10 npz: {os.path.exists(f'{DBASE}/chr10/chr10_sv_genotypes.npz')}")

files = sorted(glob.glob(f"{PERCH}/*_perm.tsv"))
chroms = sorted({os.path.basename(f).split("_cdr")[0].split("_mean")[0] for f in files},
                key=lambda c: (len(c), c))
print(f"\n[interim] result chroms: {chroms}")

print("\n=== per-chrom: interim FLOOR-hits crosswalk to ADMIT ===")
print(f"{'chr':6} {'floor_hits':>10} {'mapped':>7} {'drop':>5} {'med_overlapfrac':>15} {'min':>5} {'has_var':>7}")
rows = []
var_cache = {}
for ch in chroms:
    vfile = f"{DBASE}/{ch}/{ch}_sv_variants.tsv"
    has_var = os.path.exists(vfile)
    if not has_var:
        print(f"{ch:6} {'--':>10} {'--':>7} {'--':>5} {'NO base var table':>15}")
        continue
    var = pd.read_csv(vfile, sep="\t", dtype=str)
    var["ids"] = var["id"].map(parse_ids)
    var["key"] = var["pos"].astype(str)+"|"+var["svtype"]+"|"+var["svlen"].astype(str)
    key2ids = dict(zip(var["key"], var["ids"]))
    inv = inv_by_chr.get(ch, {})
    hits = []
    for f in glob.glob(f"{PERCH}/{ch}_*_perm.tsv"):
        d = pd.read_csv(f, sep="\t")
        fh = d[d["p_emp_perm"] <= 2.0e-6]
        hits.append(fh)
    if not hits:
        continue
    H = pd.concat(hits, ignore_index=True).drop_duplicates(subset=["pos","svtype","svlen"])
    mapped, drop, ofracs = 0, 0, []
    for _, h in H.iterrows():
        key = f"{int(h['pos'])}|{h['svtype']}|{int(h['svlen'])}"
        members = key2ids.get(key, set())
        if not members:
            continue
        nadm = sum(1 for i in members if i in inv)
        if nadm > 0:
            mapped += 1; ofracs.append(nadm/len(members))
        else:
            drop += 1
    med = np.median(ofracs) if ofracs else float("nan")
    mn  = np.min(ofracs) if ofracs else float("nan")
    print(f"{ch:6} {len(H):>10} {mapped:>7} {drop:>5} {med:>15.3f} {mn:>5.2f} {str(has_var):>7}")
    rows.append((ch, len(H), mapped, drop, med))

tot_h = sum(r[1] for r in rows); tot_m = sum(r[2] for r in rows); tot_d = sum(r[3] for r in rows)
print(f"\n[TOTAL] interim floor-hits={tot_h}  map-to-ADMIT={tot_m}  would-drop(non-ADMIT)={tot_d}")
