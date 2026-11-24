#!/usr/bin/env python
"""audit_rescue.py — Consolidate the Option-B rescue sweep with the production stubs."""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import os as _censv_os
CENSV_ROOT = _censv_os.environ["CENSV_ROOT"]

PROJ = Path(CENSV_ROOT + "/analyses/sv_catalog")
SWEEP_ROOT = PROJ / "output/ref_compare_proxy_vs_dense/sweep_6519_locked_v3"
DIAG_ROOT = SWEEP_ROOT / "diag"
RESCUE_ROOT = SWEEP_ROOT / "rescue"
MANIFEST_ROOT = SWEEP_ROOT / "manifest"

NULL_BAND = 0.02
SOFT_HIT = 0.1


STUB_STATUSES = {"R_AUG_TOO_LARGE_AT_64G", "ORIENT_GAP", "PHASE_A_MISSING"}
NA_SVPIP_SENTINELS = {"", "NA", "NaN", "nan", "NAN", "Inf", "-Inf", "inf", "-inf"}


def read_diag(p: Path) -> dict | None:
    with open(p, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
    if len(rows) != 1:
        print(f"[skip] malformed diag rows={len(rows)}: {p}", file=sys.stderr)
        return None
    return rows[0]


def classify_rescue_row(rescue_path: Path, rescue_group: str) -> dict:
    if not rescue_path.exists():
        return {
            "classification": "R_AUG_TOO_LARGE_AT_64G_POSTHOC",
            "rescue_group": rescue_group,
            "sv_pip": float("nan"),
        }
    d = read_diag(rescue_path)
    if d is None:
        return {
            "classification": "MALFORMED",
            "rescue_group": rescue_group,
            "sv_pip": float("nan"),
        }
    status = d.get("status", "")
    if status in STUB_STATUSES:
        return dict(classification=status,
                    rescue_group=d.get("rescue_group", rescue_group),
                    sv_pip=float("nan"))
    if status != "" and status != "OK":
        return dict(classification="UNKNOWN_STATUS",
                    rescue_group=d.get("rescue_group", rescue_group),
                    sv_pip=float("nan"))
    sv_pip_str = d.get("sv_pip", "")
    if sv_pip_str in NA_SVPIP_SENTINELS:
        return dict(classification="UNCONVERGED_RESCUE",
                    rescue_group=d.get("rescue_group", rescue_group),
                    sv_pip=float("nan"))
    converged = d.get("converged", "")
    warn = d.get("warn_unconverged", "")
    if converged != "TRUE" or warn not in ("", "FALSE"):
        return dict(classification="UNCONVERGED_RESCUE",
                    rescue_group=d.get("rescue_group", rescue_group),
                    sv_pip=float(sv_pip_str))
    g = d.get("rescue_group", rescue_group)
    cls = "RESCUED_OK_WIDE" if g == "widest_decile" else "RESCUED_OK"
    return dict(classification=cls, rescue_group=g, sv_pip=float(sv_pip_str))


def audit_pheno(pheno: str) -> pd.DataFrame:
    mani = pd.read_csv(MANIFEST_ROOT / f"MISSING_RESCUE_{pheno}.tsv", sep="\t")
    classifications = []
    for _, row in mani.iterrows():
        rescue_path = RESCUE_ROOT / pheno / row["sv_chrom"] / f"{row['sv_locus']}__{row['arm']}__rescue_diag.tsv"
        cls = classify_rescue_row(rescue_path, row["rescue_group"])
        classifications.append({
            "phenotype": pheno,
            "pair_id": row["pair_id"],
            "sv_locus": row["sv_locus"],
            "arm": row["arm"],
            "sv_chrom": row["sv_chrom"],
            "n_panel_snp": row["n_panel_snp"],
            **cls,
        })
    out = pd.DataFrame(classifications)
    return out


def emit_summary(audit: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pheno, sub in audit.groupby("phenotype"):
        cls_counts = Counter(sub["classification"])
        ok_rows = sub[sub["classification"].isin(["RESCUED_OK", "RESCUED_OK_WIDE"])]
        n_nan_sv_pip_in_ok = int(ok_rows["sv_pip"].isna().sum()) if len(ok_rows) else 0
        max_sv_pip = (
            float(ok_rows["sv_pip"].max(skipna=False)) if len(ok_rows) > 0 else float("nan")
        )
        median_sv_pip = (
            float(ok_rows["sv_pip"].median(skipna=False)) if len(ok_rows) > 0 else float("nan")
        )
        ok_pip_real = ok_rows["sv_pip"].dropna()
        n_hits_soft = int((ok_pip_real >= SOFT_HIT).sum())
        n_above_null = int((ok_pip_real >= NULL_BAND).sum())
        rows.append({
            "phenotype": pheno,
            "n_attempted": len(sub),
            "n_rescued_ok": int(cls_counts.get("RESCUED_OK", 0) + cls_counts.get("RESCUED_OK_WIDE", 0)),
            "n_rescued_ok_widest": int(cls_counts.get("RESCUED_OK_WIDE", 0)),
            "n_unconverged_rescue": int(cls_counts.get("UNCONVERGED_RESCUE", 0)),
            "n_stub_64G_preflight": int(cls_counts.get("R_AUG_TOO_LARGE_AT_64G", 0)),
            "n_stub_64G_posthoc": int(cls_counts.get("R_AUG_TOO_LARGE_AT_64G_POSTHOC", 0)),
            "n_orient_gap": int(cls_counts.get("ORIENT_GAP", 0)),
            "n_phase_a_missing": int(cls_counts.get("PHASE_A_MISSING", 0)),
            "n_malformed": int(cls_counts.get("MALFORMED", 0)),
            "n_unknown_status": int(cls_counts.get("UNKNOWN_STATUS", 0)),
            "n_nan_sv_pip_in_ok": n_nan_sv_pip_in_ok,
            "max_sv_pip_rescued": max_sv_pip,
            "median_sv_pip_rescued": median_sv_pip,
            "n_hits_soft_0_10": n_hits_soft,
            "n_above_null_0_02": n_above_null,
        })
    return pd.DataFrame(rows)


def _stratified_completion(audit_pheno: pd.DataFrame) -> dict:
    out = {}
    for g, sub in audit_pheno.groupby("rescue_group"):
        n_attempted = int(len(sub))
        n_ok = int(
            sub["classification"].isin(["RESCUED_OK", "RESCUED_OK_WIDE"]).sum()
        )
        out[g] = (n_ok, n_attempted, (n_ok / n_attempted) if n_attempted else 0.0)
    return out


def emit_verdict(
    summary: pd.DataFrame,
    audit: pd.DataFrame,
    completion_threshold: float,
) -> str:
    MNAR_FAIL_FRAC = 0.05
    lines = ["RESCUE_AUDIT_VERDICT"]
    for _, r in summary.iterrows():
        sub = audit[audit["phenotype"] == r["phenotype"]]
        per_group = _stratified_completion(sub)
        worst_group = min(per_group.items(), key=lambda kv: kv[1][2]) if per_group else ("", (0, 0, 0))
        global_completion = (
            r["n_rescued_ok"] / r["n_attempted"] if r["n_attempted"] else 0
        )
        post_hoc_fail_frac = (
            r["n_stub_64G_posthoc"] / r["n_attempted"] if r["n_attempted"] else 0
        )
        all_strata_pass = all(
            comp >= completion_threshold for (_, _, comp) in per_group.values()
        )
        nan_present = int(r["n_nan_sv_pip_in_ok"]) > 0

        if r["n_hits_soft_0_10"] > 0:
            v = "RESCUE_HIT"
        elif nan_present:
            v = "NULL_CONFIRMED_NAN_PRESENT"
        elif (
            global_completion >= completion_threshold
            and all_strata_pass
            and r["n_above_null_0_02"] == 0
        ):
            if post_hoc_fail_frac > MNAR_FAIL_FRAC:
                v = "NULL_CONFIRMED_MNAR_LARGE_PANEL_GAP"
            else:
                v = "NULL_CONFIRMED"
        elif global_completion >= completion_threshold and all_strata_pass:
            v = "WEAK_HIT"
        elif global_completion >= completion_threshold:
            v = "INCOMPLETE_STRATIFIED"
        else:
            v = "INCOMPLETE"
        per_group_str = " ".join(
            f"{g}={n_ok}/{n_att}={comp:.2f}" for g, (n_ok, n_att, comp) in per_group.items()
        )
        lines.append(
            f"{r['phenotype']}: VERDICT={v}  attempted={r['n_attempted']}  "
            f"rescued_ok={r['n_rescued_ok']}  global_completion={global_completion:.3f}  "
            f"worst_stratum={worst_group[0]}@{worst_group[1][2]:.3f}  "
            f"posthoc_fail_frac={post_hoc_fail_frac:.3f}  "
            f"max_sv_pip={r['max_sv_pip_rescued']:.3e}  "
            f"n_hits>=0.10={r['n_hits_soft_0_10']}  n>=0.02={r['n_above_null_0_02']}  "
            f"n_nan_sv_pip={r['n_nan_sv_pip_in_ok']}  "
            f"per_group=[{per_group_str}]"
        )
    return "\n".join(lines) + "\n"


def write_atomic_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    df.to_csv(tmp, sep="\t", index=False)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    print(f"[write] {path}  rows={len(df)}", file=sys.stderr)


def write_atomic_text(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    tmp.write_text(text)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    print(f"[write] {path}", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phenos", nargs="+", required=True)
    ap.add_argument("--out-audit", type=Path, required=True)
    ap.add_argument("--out-summary", type=Path, required=True)
    ap.add_argument("--out-verdict", type=Path, required=True)
    ap.add_argument("--completion-threshold", type=float, default=0.80)
    args = ap.parse_args()

    parts = [audit_pheno(p) for p in args.phenos]
    audit = pd.concat(parts, ignore_index=True)
    write_atomic_tsv(audit, args.out_audit)

    summary = emit_summary(audit)
    write_atomic_tsv(summary, args.out_summary)

    verdict = emit_verdict(summary, audit, args.completion_threshold)
    write_atomic_text(verdict, args.out_verdict)

    print(verdict, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
