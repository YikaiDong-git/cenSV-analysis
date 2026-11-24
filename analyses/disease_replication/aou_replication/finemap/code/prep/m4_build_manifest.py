#!/usr/bin/env python3
"""M4 — build the SuSiE fine-map target manifest."""
import csv, os, collections

REP = "/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication/aou_replication/results_vds"
MCC = {"MCC2": REP + "/step30_phewide", "MCC1": REP + "/step30_phewide_mcc1"}
STRATA = {
    "all-anc": "disease_matrices_firth/disease_hits_final.tsv",
    "EUR": "disease_matrices_firth/per_ancestry/disease_hits_final_EUR.tsv",
    "AFR": "disease_matrices_firth/per_ancestry/disease_hits_final_AFR.tsv",
    "AMR": "disease_matrices_firth/per_ancestry/disease_hits_final_AMR.tsv",
}
Q_MAX, CV_MIN = 0.10, 0.10
OUTDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/targets"
os.makedirs(OUTDIR, exist_ok=True)

def load_cv(path):
    cv = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try: cv[r["sv_locus"]] = float(r["cv_r2"])
            except Exception: pass
    return cv

cv = load_cv(MCC["MCC2"] + "/sv_cv_r2_all.tsv")

Fit = lambda: {"svs": {}, "mcc": set(), "label": "", "cases": 0, "cvmax": 0.0, "qmin": 1.0, "tiers": set()}
fits = collections.defaultdict(Fit)
cell_count = collections.Counter()

for mcc, tree in MCC.items():
    for stratum, rel in STRATA.items():
        path = os.path.join(tree, rel)
        if not os.path.exists(path):
            print(f"  MISSING {mcc}/{stratum}: {path}"); continue
        with open(path) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                q = r.get("q_bh_tier", "")
                if q in ("", "NA", "nan"): continue
                q = float(q)
                if q >= Q_MAX: continue
                sv = r["sv_locus"]; c = cv.get(sv)
                if c is None or c < CV_MIN: continue
                cell_count[(mcc, stratum)] += 1
                chrom = sv.split("_")[0]
                phe = r["phecode"]
                k = (chrom, phe, stratum)
                fit = fits[k]
                fit["label"] = r.get("label", "")
                fit["mcc"].add(mcc)
                fit["svs"][sv] = max(fit["svs"].get(sv, 0.0), c)
                try: fit["cases"] = max(fit["cases"], int(float(r.get("cases", 0))))
                except Exception: pass
                fit["cvmax"] = max(fit["cvmax"], c)
                fit["qmin"] = min(fit["qmin"], q)
                fit["tiers"].add(r.get("tier_label", ""))

rows = []
for (chrom, phe, stratum), fit in fits.items():
    svs = sorted(fit["svs"], key=lambda s: -fit["svs"][s])
    rows.append({
        "fit_id": f"{chrom}.{phe}.{stratum}",
        "chrom": chrom, "phecode": phe, "label": fit["label"], "stratum": stratum,
        "primary_mcc": "MCC2" if "MCC2" in fit["mcc"] else "MCC1",
        "mcc_hits": ";".join(sorted(fit["mcc"])),
        "n_sv": len(svs), "cases_max": fit["cases"],
        "cv_r2_max": round(fit["cvmax"], 4), "q_min": f'{fit["qmin"]:.3g}',
        "tier_labels": ";".join(sorted(t for t in fit["tiers"] if t)),
        "sv_ids": ",".join(svs),
        "window_start": "", "window_end": "",
    })
def chr_key(c):
    v = c.replace("chr", ""); return (int(v), "") if v.isdigit() else (99, v)
rows.sort(key=lambda r: (chr_key(r["chrom"]), r["stratum"], r["phecode"]))
fields = ["fit_id","chrom","phecode","label","stratum","primary_mcc","mcc_hits","n_sv",
          "cases_max","cv_r2_max","q_min","tier_labels","sv_ids","window_start","window_end"]
mpath = os.path.join(OUTDIR, "target_manifest.tsv")
with open(mpath, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader()
    for r in rows: w.writerow(r)

by_stratum = collections.Counter(r["stratum"] for r in rows)
by_stratum_mcc = collections.Counter((r["stratum"], r["primary_mcc"]) for r in rows)
by_chrom = collections.Counter(r["chrom"] for r in rows)
phecodes = {r["phecode"] for r in rows}
allsv = set();
for r in rows: allsv.update(r["sv_ids"].split(","))
lines = []
lines.append(f"FITS = {len(rows)}  | phecodes = {len(phecodes)} | distinct SVs = {len(allsv)} | chrom-windows = {len(by_chrom)}")
lines.append(f"loose CELL counts (mcc,stratum) [cross-check vs verification]: " +
             ", ".join(f"{m}/{s}={n}" for (m,s),n in sorted(cell_count.items())))
lines.append("FITS by stratum: " + ", ".join(f"{s}={n}" for s,n in sorted(by_stratum.items())))
lines.append("FITS by stratum x primary_mcc: " + ", ".join(f"{s}:{m}={n}" for (s,m),n in sorted(by_stratum_mcc.items())))
lines.append("FITS per chrom: " + ", ".join(f"{c}={by_chrom[c]}" for c in sorted(by_chrom, key=chr_key)))
spath = os.path.join(OUTDIR, "manifest_summary.txt")
with open(spath, "w") as f: f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwrote {mpath} ({len(rows)} fits)\nwrote {spath}")
