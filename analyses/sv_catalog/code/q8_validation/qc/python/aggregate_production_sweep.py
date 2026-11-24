"""aggregate_production_sweep.py — R6."""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def gather_diag_files(sweep_root: Path) -> list[Path]:
    return sorted(sweep_root.rglob("*__diag.tsv"))


def read_diag(p: Path) -> dict[str, str]:
    with open(p, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
    assert len(rows) == 1, f"expected 1 row in {p}, got {len(rows)}"
    return rows[0]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-root", required=True,
                    help="Path to the sweep output dir (e.g. output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3)")
    ap.add_argument("--out-aggregate", required=True)
    ap.add_argument("--out-pip-summary", required=True)
    args = ap.parse_args()

    sweep_root = Path(args.sweep_root)
    if not sweep_root.is_dir():
        print(f"ERROR: sweep root not a dir: {sweep_root}", file=sys.stderr)
        return 2

    diag_paths = gather_diag_files(sweep_root)
    print(f"[1] found {len(diag_paths)} diag.tsv files under {sweep_root}")

    rows: list[dict[str, str]] = []
    all_cols: list[str] = []
    seen_cols: set[str] = set()
    for p in diag_paths:
        r = read_diag(p)
        r["_source_path"] = str(p.relative_to(sweep_root))
        for c in r.keys():
            if c not in seen_cols:
                all_cols.append(c)
                seen_cols.add(c)
        rows.append(r)

    if not rows:
        print("WARN: no rows", file=sys.stderr)
        return 0

    out_ag = Path(args.out_aggregate)
    tmp = out_ag.parent / f"{out_ag.name}.tmp.{os.getpid()}"
    out_ag.parent.mkdir(parents=True, exist_ok=True)
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=all_cols, delimiter="\t",
                           quoting=csv.QUOTE_NONE, escapechar="\\")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in all_cols})
    os.replace(tmp, out_ag)
    print(f"[2] wrote aggregate: {out_ag}  ({len(rows)} rows, {len(all_cols)} cols)")

    has_tier = "s_rss_tier" in all_cols
    pip_by_strata: dict[tuple[str, str], list[float]] = defaultdict(list)
    for r in rows:
        pip = r.get("sv_pip", "")
        if pip in ("", "NA", "NaN", "nan"):
            continue
        arm = r.get("arm", "unknown")
        tier = r.get("s_rss_tier", "ALL") if has_tier else "ALL"
        if tier in ("", "NA", "NaN", "nan"):
            tier = "ALL"
        pip_by_strata[(tier, arm)].append(float(pip))

    summ_rows = []
    for (tier, arm), vec in sorted(pip_by_strata.items()):
        n = len(vec)
        if n == 0:
            continue
        arr = np.asarray(vec, dtype=np.float64)
        mean_p = float(arr.mean())
        median_p = float(np.median(arr))
        q90 = float(np.quantile(arr, 0.90))
        q99 = float(np.quantile(arr, 0.99))
        max_p = float(arr.max())
        n_ge_05 = int((arr >= 0.5).sum())
        n_ge_09 = int((arr >= 0.9).sum())
        summ_rows.append({
            "qc_tier": tier, "arm": arm, "n": n,
            "mean_sv_pip": f"{mean_p:.6g}",
            "median_sv_pip": f"{median_p:.6g}",
            "q90_sv_pip": f"{q90:.6g}",
            "q99_sv_pip": f"{q99:.6g}",
            "max_sv_pip": f"{max_p:.6g}",
            "n_sv_pip_ge_0.5": n_ge_05,
            "n_sv_pip_ge_0.9": n_ge_09,
        })

    out_summ = Path(args.out_pip_summary)
    summ_cols = ["qc_tier", "arm", "n", "mean_sv_pip", "median_sv_pip",
                 "q90_sv_pip", "q99_sv_pip", "max_sv_pip",
                 "n_sv_pip_ge_0.5", "n_sv_pip_ge_0.9"]
    tmp = out_summ.parent / f"{out_summ.name}.tmp.{os.getpid()}"
    out_summ.parent.mkdir(parents=True, exist_ok=True)
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=summ_cols, delimiter="\t",
                           quoting=csv.QUOTE_NONE, escapechar="\\")
        w.writeheader()
        for r in summ_rows:
            w.writerow(r)
    os.replace(tmp, out_summ)
    print(f"[3] wrote PIP summary: {out_summ}  ({len(summ_rows)} strata)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
