#!/usr/bin/env python3
"""run_finemap_target — orchestrate ONE fit (chrom x phecode x stratum) end to end."""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time

import finemap_core as core

HERE = os.path.dirname(os.path.abspath(__file__))
FM = core.FM
MANIFEST = f"{FM}/targets_ext2/target_manifest.tsv"
RESULTS = f"{FM}/results_ext2"
FIT_R = f"{HERE}/finemap_fit.R"
import _ext2_guard  # noqa: F401,E402 - fail-closed isolation assertions (side effects on import)


def manifest_row(fit_id=None, chrom=None, phecode=None, stratum=None):
    with open(MANIFEST) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if fit_id and r["fit_id"] == fit_id:
                return r
            if (not fit_id) and r["chrom"] == chrom and r["phecode"] == phecode and r["stratum"] == stratum:
                return r
    raise SystemExit(f"fit not found in manifest: {fit_id or (chrom, phecode, stratum)}")


def run_one(row, panels, force=False, L=10):
    fit_id = row["fit_id"]
    chrom, phecode, stratum = row["chrom"], row["phecode"], row["stratum"]
    mcc = "mcc2" if row["primary_mcc"] == "MCC2" else "mcc1"
    hit_svs = [s for s in row["sv_ids"].split(",") if s]
    out_root = f"{RESULTS}/fits/{row['stratum']}/{fit_id}"
    os.makedirs(out_root, exist_ok=True)
    rows = []
    for mode in panels:
        od = f"{out_root}/{mode}"
        t0 = time.time()
        core.build_fit(chrom, phecode, stratum, mcc, hit_svs, od, panel_mode=mode, force=force)
        subprocess.run(["Rscript", FIT_R, "--dir", od, "--L", str(L)], check=True)
        summ = f"{od}/susie_fit_summary.tsv"
        s = {}
        if os.path.exists(summ):
            with open(summ) as fh:
                s = list(csv.DictReader(fh, delimiter="\t"))[0]
        s.update({"fit_id": fit_id, "panel_mode": mode, "primary_mcc": row["primary_mcc"],
                  "label": row.get("label", ""), "elapsed_s": round(time.time() - t0, 1)})
        rows.append(s)
        print(f"[{fit_id}/{mode}] done ({time.time()-t0:.0f}s)", flush=True)
    if rows:
        keys = ["fit_id", "panel_mode", "primary_mcc", "label", "n_variants", "n_sv", "n_snp",
                "n_present", "n_eff", "L", "n_cs_indiv", "n_cs_rss_lpm", "n_cs_rss_logit",
                "jaccard_indiv_vs_rsslpm", "pipcor_indiv_vs_rsslpm",
                "jaccard_indiv_vs_rsslogit", "pipcor_indiv_vs_rsslogit",
                "value_indiv", "value_rss_lpm", "value_rss_logit",
                "sum_sv_pip_indiv", "sum_sv_pip_rss_lpm", "sum_sv_pip_rss_logit",
                "top_sv_pip_indiv", "top_sv_pip_rss_logit", "max_abs_sv_snp_r",
                "converged_indiv", "converged_rss_lpm", "converged_rss_logit", "elapsed_s"]
        with open(f"{out_root}/fit_summary_bymode.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit-id")
    ap.add_argument("--chrom"); ap.add_argument("--phecode"); ap.add_argument("--stratum")
    ap.add_argument("--panels", default="hit_only,all_window")
    ap.add_argument("--L", type=int, default=10)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    row = manifest_row(a.fit_id, a.chrom, a.phecode, a.stratum)
    panels = [p for p in a.panels.split(",") if p]
    rows = run_one(row, panels, force=a.force, L=a.L)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
