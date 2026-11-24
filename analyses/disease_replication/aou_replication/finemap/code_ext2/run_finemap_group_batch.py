#!/usr/bin/env python3
"""run_finemap_group_batch — BATCHED-SCORE production runner."""
from __future__ import annotations
import argparse
import csv
import os
import subprocess
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

_W = int(os.environ.get("FINEMAP_GROUP_WORKERS", "3"))
_BLAS = os.environ.get("FINEMAP_BLAS_THREADS", str(max(1, 16 // max(_W, 1))))
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, _BLAS)

import finemap_core as core  # noqa: E402
import finemap_group as grp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = f"{core.FM}/targets_ext2/target_manifest.tsv"
RESULTS = f"{core.FM}/results_ext2"
FIT_R = f"{HERE}/finemap_fit.R"
import _ext2_guard  # noqa: F401,E402 - fail-closed isolation assertions (side effects on import)
BATCH_LOG = f"{RESULTS}/_group_batch_status.tsv"
STRATUM_ORDER = {"AMR": 0, "AFR": 1, "EUR": 2, "all-anc": 3}


def load_groups(strata, panels, chroms=None):
    groups = defaultdict(list)
    with open(MANIFEST) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if strata and r["stratum"] not in strata:
                continue
            if chroms and r["chrom"] not in chroms:
                continue
            mcc = "mcc2" if r["primary_mcc"] == "MCC2" else "mcc1"
            hit = [s for s in r["sv_ids"].split(",") if s]
            for mode in panels:
                groups[(r["chrom"], r["stratum"])].append({
                    "phecode": r["phecode"], "mcc": mcc, "hit_sv_ids": hit, "panel_mode": mode,
                    "snp_select": "assoc_tag", "fit_id": r["fit_id"],
                    "out_dir": f"{RESULTS}/fits/{r['stratum']}/{r['fit_id']}/{mode}"})
    return groups


def _susie(od):
    cp = subprocess.run(["Rscript", FIT_R, "--dir", od, "--L", "10"],
                        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if cp.returncode != 0:
        open(f"{od}/ERROR.txt", "w").write("susie failed:\n" + cp.stderr.decode("utf-8", "replace")[-800:])
        return False
    for b in ("Mres.bin", "yres.bin"):
        try:
            os.remove(f"{od}/{b}")
        except OSError:
            pass
    open(f"{od}/_SUCCESS", "w").write("ok\n")
    return True


def run_group(item):
    (chrom, stratum), fits = item
    todo = [f for f in fits if not os.path.exists(f"{f['out_dir']}/_SUCCESS")
            and not os.path.exists(f"{f['out_dir']}/_SKIP")]
    if not todo:
        return {"chrom": chrom, "stratum": stratum, "built": 0, "susie": 0, "skip": 0,
                "cached": len(fits), "fail": 0, "sec": 0.0}
    t0 = time.time()
    built = susie_ok = skipped = failed = 0
    try:
        res = grp.build_group(chrom, stratum, todo)
    except Exception as e:  # noqa: BLE001 - a group-level failure marks its fits failed, batch continues
        for f in todo:
            os.makedirs(f["out_dir"], exist_ok=True)
            open(f"{f['out_dir']}/ERROR.txt", "w").write(f"build_group failed: {e}\n")
        return {"chrom": chrom, "stratum": stratum, "built": 0, "susie": 0, "skip": 0,
                "cached": len(fits) - len(todo), "fail": len(todo), "sec": round(time.time() - t0, 1)}
    by_key = {}
    for r in res:
        by_key[(r.get("phecode"), r.get("panel_mode"))] = r
    for f in todo:
        r = by_key.get((f["phecode"], f["panel_mode"]), {})
        st = r.get("status", "MISSING")
        if st == "built":
            built += 1
            if _susie(f["out_dir"]):
                susie_ok += 1
            else:
                failed += 1
        elif st.startswith("SKIP"):
            os.makedirs(f["out_dir"], exist_ok=True)
            open(f"{f['out_dir']}/_SKIP", "w").write(st + "\n")
            skipped += 1
        else:
            failed += 1
    return {"chrom": chrom, "stratum": stratum, "built": built, "susie": susie_ok, "skip": skipped,
            "cached": len(fits) - len(todo), "fail": failed, "sec": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strata", default="")
    ap.add_argument("--panels", default="hit_only,all_window")
    ap.add_argument("--limit-groups", type=int, default=0, help="cap number of groups (smoke)")
    ap.add_argument("--chroms", default="", help="comma list, e.g. chr19 (scheduling only)")
    ap.add_argument("--exclude-chroms", default="", help="comma list to skip (scheduling only)")
    a = ap.parse_args()
    strata = [s for s in a.strata.split(",") if s] or None
    panels = tuple(p for p in a.panels.split(",") if p)
    chroms = [c for c in a.chroms.split(",") if c] or None
    groups = load_groups(strata, panels, chroms)
    if a.exclude_chroms:
        drop = {c for c in a.exclude_chroms.split(",") if c}
        groups = {k: v for k, v in groups.items() if k[0] not in drop}
    items = sorted(groups.items(), key=lambda kv: (STRATUM_ORDER.get(kv[0][1], 9), -len(kv[1])))
    if a.limit_groups:
        items = items[:a.limit_groups]
    nfits = sum(len(v) for _, v in items)
    print(f"GROUP BATCH: {len(items)} groups, {nfits} fit-modes, {_W} workers x {_BLAS} BLAS threads", flush=True)

    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=_W) as ex:
        futs = {ex.submit(run_group, it): it[0] for it in items}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            rows.append(r)
            print(f"  [{i}/{len(items)}] {r['chrom']}.{r['stratum']}: built={r['built']} susie={r['susie']} "
                  f"skip={r['skip']} cached={r['cached']} fail={r['fail']} ({r['sec']:.0f}s) "
                  f"| elapsed {time.time()-t0:.0f}s", flush=True)
            with open(BATCH_LOG + ".tmp", "w") as fh:
                w = csv.DictWriter(fh, fieldnames=["chrom", "stratum", "built", "susie", "skip", "cached", "fail", "sec"],
                                   delimiter="\t")
                w.writeheader()
                for rr in rows:
                    w.writerow(rr)
            os.replace(BATCH_LOG + ".tmp", BATCH_LOG)
    tot = {k: sum(r[k] for r in rows) for k in ("built", "susie", "skip", "cached", "fail")}
    print(f"GROUP BATCH DONE: {tot} of {nfits} fit-modes ({time.time()-t0:.0f}s) -> {BATCH_LOG}", flush=True)


if __name__ == "__main__":
    main()
