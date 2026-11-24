#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

_BLAS = os.environ.get("FINEMAP_BLAS_THREADS", "1")
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, _BLAS)

import finemap_core as core  # noqa: E402 - must follow the BLAS-thread pin above

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = f"{core.FM}/targets_ext2/target_manifest.tsv"
RESULTS = f"{core.FM}/results_ext2"
FIT_R = f"{HERE}/finemap_fit.R"
import _ext2_guard  # noqa: F401,E402 - fail-closed isolation assertions (side effects on import)
BATCH_LOG = f"{RESULTS}/_batch_status.tsv"

WORKERS = {"AMR": 10, "AFR": 8, "EUR": 5, "all-anc": 3}
STRATUM_ORDER = {"AMR": 0, "AFR": 1, "EUR": 2, "all-anc": 3}


def load_fits(strata, panels):
    out = []
    with open(MANIFEST) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if strata and r["stratum"] not in strata:
                continue
            for mode in panels:
                out.append((r, mode))
    return out


def run_one(task):
    r, mode = task
    fit_id = r["fit_id"]
    od = f"{RESULTS}/fits/{r['stratum']}/{fit_id}/{mode}"
    succ = f"{od}/_SUCCESS"
    if os.path.exists(succ):
        return (fit_id, mode, r["stratum"], "cached", 0.0)
    t0 = time.time()
    mcc = "mcc2" if r["primary_mcc"] == "MCC2" else "mcc1"
    hit = [s for s in r["sv_ids"].split(",") if s]
    try:
        core.build_fit(r["chrom"], r["phecode"], r["stratum"], mcc, hit, od, panel_mode=mode)
        cp = subprocess.run(["Rscript", FIT_R, "--dir", od, "--L", "10"],
                            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for b in ("Mres.bin", "yres.bin"):
            try:
                os.remove(f"{od}/{b}")
            except OSError:
                pass
        with open(succ, "w") as fh:
            fh.write(f"ok {time.time()-t0:.0f}s\n")
        return (fit_id, mode, r["stratum"], "done", time.time() - t0)
    except subprocess.CalledProcessError as e:
        msg = (e.stderr or b"").decode("utf-8", "replace")[-500:]
        if os.path.isdir(od):
            open(f"{od}/ERROR.txt", "w").write("Rscript failed:\n" + msg + "\n")
        return (fit_id, mode, r["stratum"], f"FAIL_R", time.time() - t0)
    except Exception as e:  # noqa: BLE001 - capture any build failure per-fit, keep the batch going
        if os.path.isdir(od):
            open(f"{od}/ERROR.txt", "w").write(f"build failed: {e}\n")
        return (fit_id, mode, r["stratum"], f"FAIL: {str(e)[:120]}", time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strata", default="", help="comma list (AMR,EUR,all-anc); empty = all")
    ap.add_argument("--panels", default="hit_only,all_window")
    ap.add_argument("--limit", type=int, default=0, help="cap total fit-modes (smoke test)")
    a = ap.parse_args()
    strata = [s for s in a.strata.split(",") if s] or None
    panels = tuple(p for p in a.panels.split(",") if p)
    tasks = load_fits(strata, panels)
    if a.limit:
        tasks = tasks[:a.limit]

    by_strat = defaultdict(list)
    for t in tasks:
        by_strat[t[0]["stratum"]].append(t)
    total = len(tasks)
    print(f"BATCH: {total} fit-modes over strata {dict((s, len(v)) for s, v in by_strat.items())}", flush=True)

    t0 = time.time()
    done = cached = failed = 0
    rows = []
    for strat in sorted(by_strat, key=lambda s: STRATUM_ORDER.get(s, 9)):
        ts = by_strat[strat]
        w = WORKERS.get(strat, 3)
        print(f"=== stratum {strat}: {len(ts)} fit-modes, {w} workers ===", flush=True)
        with ProcessPoolExecutor(max_workers=w) as ex:
            futs = {ex.submit(run_one, t): t for t in ts}
            for i, fut in enumerate(as_completed(futs), 1):
                fid, mode, st, status, dt = fut.result()
                rows.append({"fit_id": fid, "panel_mode": mode, "stratum": st, "status": status, "seconds": round(dt, 1)})
                if status == "done":
                    done += 1
                elif status == "cached":
                    cached += 1
                else:
                    failed += 1
                if i % 5 == 0 or status.startswith("FAIL"):
                    print(f"  [{strat} {i}/{len(ts)}] {fid}/{mode}: {status} ({dt:.0f}s) | "
                          f"tot done={done} cached={cached} fail={failed} ({time.time()-t0:.0f}s)", flush=True)
                with open(BATCH_LOG + ".tmp", "w") as fh:
                    wtr = csv.DictWriter(fh, fieldnames=["fit_id", "panel_mode", "stratum", "status", "seconds"], delimiter="\t")
                    wtr.writeheader()
                    for rr in rows:
                        wtr.writerow(rr)
                os.replace(BATCH_LOG + ".tmp", BATCH_LOG)
    print(f"BATCH DONE: {done} done, {cached} cached, {failed} failed of {total} ({time.time()-t0:.0f}s) -> {BATCH_LOG}", flush=True)


if __name__ == "__main__":
    main()
