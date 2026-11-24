#!/usr/bin/env python3
"""M4-ext2 — the SuSiE fine-map target manifest for the 11,129-predictor arm."""
import os
import csv
import sys
import argparse
import collections

REP = ("/home/dataproc/workspace/working/ydon268/Project1_Centromere/analyses/"
       "disease_replication/aou_replication/results_vds")
FM = REP.replace("/results_vds", "/finemap")

PUB = {"MCC2": f"{REP}/step30_phewide", "MCC1": f"{REP}/step30_phewide_mcc1"}
PUB_CV = f"{REP}/step30_phewide/sv_cv_r2_all.tsv"
EXT = {"MCC2": f"{REP}/step30_phewide_ext2", "MCC1": f"{REP}/step30_phewide_mcc1_ext2"}
EXT_CV = f"{REP}/step28_phewide_ext2/sv_cv_r2_all.tsv"

STRATA = {
    "all-anc": "disease_matrices_firth/disease_hits_final.tsv",
    "EUR": "disease_matrices_firth/per_ancestry/disease_hits_final_EUR.tsv",
    "AFR": "disease_matrices_firth/per_ancestry/disease_hits_final_AFR.tsv",
    "AMR": "disease_matrices_firth/per_ancestry/disease_hits_final_AMR.tsv",
}
Q_MAX, CV_MIN = 0.10, 0.10

GATE_EXPECT = dict(fits=170, phecodes=120, svs=304, chrom_windows=19)


def load_cv(path):
    cv = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                cv[r["sv_locus"]] = float(r["cv_r2"])
            except Exception:
                pass
    return cv


def chr_key(c):
    v = c.replace("chr", "")
    return (int(v), "") if v.isdigit() else (99, v)


def build(mcc_trees, cv_path):
    cv = load_cv(cv_path)
    Fit = (lambda: {"svs": {}, "mcc": set(), "label": "", "cases": 0,
                    "cvmax": 0.0, "qmin": 1.0, "tiers": set()})
    fits = collections.defaultdict(Fit)
    cell_count = collections.Counter()
    missing = []

    for mcc, tree in mcc_trees.items():
        for stratum, rel in STRATA.items():
            path = os.path.join(tree, rel)
            if not os.path.exists(path):
                missing.append(f"{mcc}/{stratum}")
                continue
            with open(path) as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    q = r.get("q_bh_tier", "")
                    if q in ("", "NA", "nan"):
                        continue
                    q = float(q)
                    if q >= Q_MAX:
                        continue
                    sv = r["sv_locus"]
                    c = cv.get(sv)
                    if c is None or c < CV_MIN:
                        continue
                    cell_count[(mcc, stratum)] += 1
                    k = (sv.split("_")[0], r["phecode"], stratum)
                    fit = fits[k]
                    fit["label"] = r.get("label", "")
                    fit["mcc"].add(mcc)
                    fit["svs"][sv] = max(fit["svs"].get(sv, 0.0), c)
                    try:
                        fit["cases"] = max(fit["cases"], int(float(r.get("cases", 0))))
                    except Exception:
                        pass
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
    rows.sort(key=lambda r: (chr_key(r["chrom"]), r["stratum"], r["phecode"]))
    return rows, cell_count, missing


def stats(rows):
    allsv = set()
    for r in rows:
        allsv.update(r["sv_ids"].split(","))
    return dict(fits=len(rows), phecodes=len({r["phecode"] for r in rows}),
                svs=len(allsv), chrom_windows=len({r["chrom"] for r in rows}))


FIELDS = ["fit_id", "chrom", "phecode", "label", "stratum", "primary_mcc", "mcc_hits", "n_sv",
          "cases_max", "cv_r2_max", "q_min", "tier_labels", "sv_ids", "window_start", "window_end"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-gate", action="store_true")
    a = ap.parse_args()

    if not a.no_gate:
        print("=" * 96)
        print("GATE — reproduce the PUBLISHED manifest with this code before trusting any ext2 number")
        print("=" * 96)
        prow, _, pmiss = build(PUB, PUB_CV)
        got = stats(prow)
        ok = True
        for k, want in GATE_EXPECT.items():
            good = got[k] == want
            ok &= good
            print(f"  {k:<16} expect {want:>5}   got {got[k]:>5}   {'PASS' if good else '*** FAIL ***'}")
        if pmiss:
            print(f"  missing published tables: {pmiss}")
        if not ok:
            print("\nGATE FAILED — the port does not reproduce the published manifest. "
                  "Do not trust the ext2 counts.")
            sys.exit(1)
        print("  GATE PASS — selection logic is faithful\n")

    print("=" * 96)
    print("BUILD — ext2 (11,129) target manifest")
    print("=" * 96)
    rows, cell_count, missing = build(EXT, EXT_CV)
    if missing:
        print(f"  ⚠️ MISSING ext2 hit tables: {missing}")
    st = stats(rows)

    outdir = f"{FM}/targets_ext2"
    os.makedirs(outdir, exist_ok=True)
    assert not os.path.islink(outdir), f"refusing to write into symlinked dir: {outdir}"
    mpath = f"{outdir}/target_manifest.tsv"
    with open(mpath, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    by_stratum = collections.Counter(r["stratum"] for r in rows)
    by_sm = collections.Counter((r["stratum"], r["primary_mcc"]) for r in rows)
    by_chrom = collections.Counter(r["chrom"] for r in rows)
    xrows = [r for r in rows if r["chrom"] == "chrX"]

    L = []
    L.append(f"FITS = {st['fits']}  | phecodes = {st['phecodes']} | distinct SVs = {st['svs']} "
             f"| chrom-windows = {st['chrom_windows']}")
    L.append(f"  of which chrX (NO PLOIDY HANDLING in step50 — hypothesis-only): {len(xrows)} fits, "
             f"{len({s for r in xrows for s in r['sv_ids'].split(',')})} SVs")
    L.append(f"  autosomal: {st['fits'] - len(xrows)} fits")
    L.append("loose CELL counts (mcc,stratum): " +
             ", ".join(f"{m}/{s}={n}" for (m, s), n in sorted(cell_count.items())))
    L.append("FITS by stratum: " + ", ".join(f"{s}={n}" for s, n in sorted(by_stratum.items())))
    L.append("FITS by stratum x primary_mcc: " +
             ", ".join(f"{s}:{m}={n}" for (s, m), n in sorted(by_sm.items())))
    L.append("FITS per chrom: " + ", ".join(f"{c}={by_chrom[c]}" for c in sorted(by_chrom, key=chr_key)))
    with open(f"{outdir}/manifest_summary.txt", "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))

    old_sum = f"{FM}/targets/manifest_summary.txt"
    H = []
    H.append("HEAD-TO-HEAD — fine-map target manifest, 6,513 vs 11,129")
    H.append("")
    H.append(f"{'metric':<24}{'6,513':>10}{'11,129':>10}{'ratio':>10}")
    for k, lab in (("fits", "fits"), ("phecodes", "phecodes"),
                   ("svs", "distinct SVs"), ("chrom_windows", "chrom-windows")):
        o, n = GATE_EXPECT[k], st[k]
        H.append(f"{lab:<24}{o:>10}{n:>10}{n / o:>9.2f}x")
    H.append("")
    if os.path.exists(old_sum):
        oldrows = list(csv.DictReader(open(f"{FM}/targets/target_manifest.tsv"), delimiter="\t"))
        oldk = {r["fit_id"] for r in oldrows}
        newk = {r["fit_id"] for r in rows}
        H.append(f"fit_ids in BOTH        : {len(oldk & newk)}")
        H.append(f"fit_ids only in 6,513  : {len(oldk - newk)}")
        H.append(f"fit_ids only in 11,129 : {len(newk - oldk)}")
        H.append("")
        H.append("⚠️ A fit_id present in both is NOT the same fit — its SV set is drawn from a different")
        H.append("   panel, so the competitor list and the LD differ. Overlap counts the question asked,")
        H.append("   not the answer given.")
    H.append("")
    H.append("Cost drivers for the fit stage (all_window competitor panel):")
    H.append("  6,513 arm: 5,143 distinct SVs across 21 chroms; median 242 window SVs, chr19 max 2,177")
    H.append("  11,129 arm: 7,633 distinct SVs across 22 chroms (+48% SVs, chrX new)")
    with open(f"{outdir}/headtohead_vs_6513.txt", "w") as f:
        f.write("\n".join(H) + "\n")
    print()
    print("\n".join(H))
    print(f"\nwrote {mpath}")
    print("M4_EXT2_DONE")


if __name__ == "__main__":
    main()
