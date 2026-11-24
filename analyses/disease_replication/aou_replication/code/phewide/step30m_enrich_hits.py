#!/usr/bin/env python3
import os, csv
from collections import defaultdict

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = f"{BASE}/aou_replication/results_vds/step30_phewide"
ANN = f"{BASE}/phewide_all_handoff/sv_annotation/sv_annotation_full.tsv"
HITS = f"{RV}/phewide_disease_full.tsv"
FIRTH = f"{RV}/phewide_disease_full_firthrefit.tsv"
OUT = f"{RV}/phewide_disease_hits_enriched.tsv"


def main():
    ann = {r["sv_locus"]: r for r in csv.DictReader(open(ANN), delimiter="\t")}
    firth = {(r["sv_locus"], r["phecode"]): r
             for r in csv.DictReader(open(FIRTH), delimiter="\t")}
    hits = []
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        q = r["q_bh_tier"]
        if q == "" or float(q) >= 0.05:
            continue
        hits.append(r)
    prom = defaultdict(int)
    for r in hits:
        a = ann.get(r["sv_locus"], {})
        prom[a.get("sv_key", r["sv_locus"])] += 1

    cols = ["sv_locus", "sv_key", "chrom", "svtype", "phecode", "label", "cases", "n",
            "beta_score", "p_score", "q_bh_tier", "OR_firth", "p_firth",
            "firth_converged", "firth_confirmed", "cv_r2", "n_svim_primitives",
            "block_span_mb", "discovery_set", "n_hits_this_sv"]
    hits.sort(key=lambda r: float(r["chi2_p"]))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in hits:
            sv = r["sv_locus"]
            a = ann.get(sv, {})
            fr = firth.get((sv, r["phecode"]), {})
            cvr2 = a.get("cv_r2", "")
            w.writerow([
                sv, a.get("sv_key", ""), sv.split("_")[0], a.get("svtype", ""),
                r["phecode"], r["label"], r["cases"], r["n"],
                r["beta_logOR"], r["chi2_p"], r["q_bh_tier"],
                fr.get("OR_firth", ""), fr.get("p_firth", ""),
                fr.get("firth_converged", ""), fr.get("firth_confirmed", ""),
                (f"{float(cvr2):.3f}" if cvr2 not in ("", None) else ""),
                a.get("n_svim_primitives", ""), a.get("block_span_mb", ""),
                a.get("discovery_set", ""), prom[a.get("sv_key", sv)],
            ])
    n_conv = sum(1 for r in hits
                 if firth.get((r["sv_locus"], r["phecode"]), {}).get("firth_converged") == "1")
    print(f"wrote {len(hits)} enriched hits -> {OUT}")
    print(f"converged (exact Firth): {n_conv}/{len(hits)} = {100*n_conv/len(hits):.1f}%")
    sz = os.path.getsize(OUT)
    print(f"size: {sz/1e6:.3f} MB ({'OK' if sz < 4.9e6 else 'TOO BIG'})")


if __name__ == "__main__":
    main()
