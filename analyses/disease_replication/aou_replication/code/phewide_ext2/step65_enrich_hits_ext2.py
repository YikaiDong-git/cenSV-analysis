#!/usr/bin/env python3
"""Step 65 (ext2 arm) — enrich the score-BH significant disease hits with SV annotation (← step30m)."""
import os, sys, csv
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import step51_disease_scoretest_ext2 as S

HOME = os.path.expanduser("~")
BASE = f"{HOME}/workspace/working/ydon268/Project1_Centromere/analyses/disease_replication"
RV = S.OUTA
META2 = f"{BASE}/phewide_all_handoff/full_11129_ext2/sv_metadata_ext2.tsv"
ANN_OLD = f"{BASE}/phewide_all_handoff/sv_annotation/sv_annotation_full.tsv"
HITS = f"{RV}/phewide_disease_full.tsv"
FIRTH = f"{RV}/phewide_disease_full_firthrefit.tsv"
OUT = f"{RV}/phewide_disease_hits_enriched.tsv"


def main():
    ext = {r["sv_locus"]: r for r in csv.DictReader(open(META2), delimiter="\t")}
    old = {}
    if os.path.exists(ANN_OLD):
        old = {r["sv_locus"]: r for r in csv.DictReader(open(ANN_OLD), delimiter="\t")}
    print(f"annotation: ext2 metadata {len(ext)} loci; legacy annotation {len(old)} loci", flush=True)

    mism = [k for k, r in old.items()
            if k in ext and f"{ext[k]['chrom']}_sv{ext[k]['sv_index']}" != r.get("sv_key")]
    if mism:
        raise SystemExit(f"ABORT: sv_key derivation disagrees with the legacy annotation for {len(mism)} "
                         f"loci (e.g. {mism[:3]}) — promiscuity would be wrong.")
    print(f"  sv_key derivation f'{{chrom}}_sv{{sv_index}}' verified against {len(old)} legacy loci: "
          f"0 mismatches", flush=True)

    firth = {(r["sv_locus"], r["phecode"]): r for r in csv.DictReader(open(FIRTH), delimiter="\t")}
    hits = []
    for r in csv.DictReader(open(HITS), delimiter="\t"):
        q = r["q_bh_tier"]
        if q == "" or float(q) >= 0.05:
            continue
        hits.append(r)

    def svkey(sv):
        e = ext.get(sv)
        return f"{e['chrom']}_sv{e['sv_index']}" if e else old.get(sv, {}).get("sv_key", sv)

    prom = defaultdict(int)
    for r in hits:
        prom[svkey(r["sv_locus"])] += 1

    cols = ["sv_locus", "sv_key", "chrom", "svtype", "phecode", "label", "cases", "n",
            "beta_score", "p_score", "q_bh_tier", "OR_firth", "p_firth",
            "firth_converged", "firth_confirmed", "cv_r2", "n_weight_snps", "sv_category",
            "sv_strategy", "sv_arm", "tier", "n_svim_primitives", "block_span_mb",
            "discovery_set", "ann_source", "n_hits_this_sv"]
    hits.sort(key=lambda r: float(r["chi2_p"]))
    n_ext_only = 0
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in hits:
            sv = r["sv_locus"]
            e = ext.get(sv, {})
            a = old.get(sv, {})
            fr = firth.get((sv, r["phecode"]), {})
            src = "ext2+legacy" if a else "ext2_only"
            if not a:
                n_ext_only += 1
            cvr2 = e.get("cv_r2", a.get("cv_r2", ""))
            parts = sv.split("_")
            svtype = a.get("svtype", "") or (parts[2] if len(parts) > 2 else "")
            w.writerow([
                sv, svkey(sv), e.get("chrom", sv.split("_")[0]), svtype,
                r["phecode"], r["label"], r["cases"], r["n"],
                r["beta_logOR"], r["chi2_p"], r["q_bh_tier"],
                fr.get("OR_firth", ""), fr.get("p_firth", ""),
                fr.get("firth_converged", ""), fr.get("firth_confirmed", ""),
                (f"{float(cvr2):.3f}" if cvr2 not in ("", None) else ""),
                e.get("n_weight_snps", ""), e.get("sv_category", ""), e.get("sv_strategy", ""),
                e.get("sv_arm", ""), e.get("tier", ""),
                a.get("n_svim_primitives", ""), a.get("block_span_mb", ""), a.get("discovery_set", ""),
                src, prom[svkey(sv)],
            ])
    n_conv = sum(1 for r in hits
                 if firth.get((r["sv_locus"], r["phecode"]), {}).get("firth_converged") == "1")
    n_nofirth = sum(1 for r in hits if (r["sv_locus"], r["phecode"]) not in firth)
    print(f"\nwrote {len(hits)} enriched score-BH hits -> {OUT}")
    print(f"  exact-Firth converged : {n_conv}/{len(hits)} = {100*n_conv/max(len(hits),1):.1f}%")
    print(f"  no exact-Firth row    : {n_nofirth}  (outside step51b's refit band)")
    print(f"  ann_source=ext2_only  : {n_ext_only}/{len(hits)} "
          f"({100*n_ext_only/max(len(hits),1):.1f}%) — no legacy n_svim_primitives/block_span_mb")
    sz = os.path.getsize(OUT)
    print(f"  size: {sz/1e6:.3f} MB ({'OK' if sz < 4.9e6 else 'TOO BIG — split before committing'})")


if __name__ == "__main__":
    main()
