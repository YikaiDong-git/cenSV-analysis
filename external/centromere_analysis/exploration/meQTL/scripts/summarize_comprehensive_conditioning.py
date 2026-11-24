#!/usr/bin/env python3
"""Summarize comprehensive conditioning without treating callsets as replication."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from robust_meqtl import bh_adjust


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def concat(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path, sep="\t") for path in paths if path.exists() and path.stat().st_size]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def adjust(frame: pd.DataFrame, p_column: str, q_column: str) -> None:
    frame[q_column] = np.nan
    for _, indices in frame.groupby(["panel", "window_bp"], sort=False).groups.items():
        index = list(indices)
        valid = pd.to_numeric(frame.loc[index, p_column], errors="coerce").notna()
        selected = np.asarray(index)[valid.to_numpy()]
        if len(selected):
            frame.loc[selected, q_column] = bh_adjust(frame.loc[selected, p_column].to_numpy(float))


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    root = Path(config["workflow_output"])
    task_ids = pd.read_csv(root / "manifests/analysis_tasks.tsv", sep="\t").analysis_task_id.astype(str).tolist()
    associations = concat([root / "analysis" / value / "associations.tsv.gz" for value in task_ids])
    variables = concat([root / "analysis" / value / "susie_variables.tsv.gz" for value in task_ids])
    failures = concat([root / "analysis" / value / "failures.tsv" for value in task_ids])
    pairs = pd.read_csv(root / "manifests/pairs.tsv.gz", sep="\t")
    primary = associations.loc[associations.pc_count.eq(10)].copy()
    if not variables.empty:
        merge_keys = ["pair_id", "panel", "window_bp"]
        sv_pip = variables.loc[variables.variable.eq("SV"), [*merge_keys, "pip"]].rename(columns={"pip": "susie_sv_pip"})
        primary = primary.merge(sv_pip, on=merge_keys, how="left", validate="one_to_one")
    else:
        primary["susie_sv_pip"] = np.nan
    adjust(primary, "rint_conditional_hc3_p", "conditional_hc3_bh_q")
    adjust(primary, "crossfit_residual_sv_hc3_p", "crossfit_residual_sv_bh_q")
    thresholds = config["robustness"]
    primary["window_supported"] = (
        primary.status.eq("ok")
        & primary.susie_sv_pip.ge(float(thresholds["susie_sv_pip_min"]))
        & primary.conditional_hc3_bh_q.le(float(thresholds["bh_q_max"]))
        & primary.crossfit_residual_sv_bh_q.le(float(thresholds["bh_q_max"]))
        & np.sign(primary.rint_conditional_beta).eq(np.sign(primary.crossfit_residual_sv_beta))
    )
    panel = (primary.groupby(["pair_id", "panel"], as_index=False)
             .agg(n_windows_available=("status", lambda value: int((value != "snp_callset_window_unavailable").sum())),
                  n_windows_tested=("status", lambda value: int((value == "ok").sum())),
                  n_windows_unavailable=("status", lambda value: int((value == "snp_callset_window_unavailable").sum())),
                  n_windows_supported=("window_supported", "sum"),
                  supported_windows=("window_bp", lambda value: ",".join(map(str, sorted(value[primary.loc[value.index, "window_supported"]])))),
                  min_conditional_q=("conditional_hc3_bh_q", "min"),
                  max_sv_pip=("susie_sv_pip", "max")))
    panel["panel_robust"] = panel.n_windows_supported.ge(int(thresholds["minimum_supported_windows"]))
    expected_panels = list(config["panels"])
    wide = panel.pivot(index="pair_id", columns="panel", values="panel_robust").reindex(columns=expected_panels).fillna(False)
    consensus = pairs[["pair_id", "sv_id", "feature_id", "discovery_scope"]].copy()
    robust_count = wide.sum(axis=1).reindex(consensus.pair_id).fillna(0).to_numpy(int)
    consensus["conditioning_class"] = np.select([robust_count == len(expected_panels), robust_count == 1],
                                                 ["cross_panel_robust", "single_panel"], default="not_robust")
    consensus["n_robust_sensitivity_callsets"] = robust_count
    window_summary = (primary.groupby(["panel", "window_bp", "discovery_scope"], dropna=False, as_index=False)
                      .agg(n_pairs=("pair_id", "size"), n_available=("status", lambda value: int((value != "snp_callset_window_unavailable").sum())),
                           n_unavailable=("status", lambda value: int((value == "snp_callset_window_unavailable").sum())),
                           n_tested=("status", lambda value: int((value == "ok").sum())),
                           n_supported=("window_supported", "sum"), n_susie_failures=("susie_sv_pip", lambda value: int(value.isna().sum()))))
    summary_dir = root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    write_frame(associations, summary_dir / "all_pc_associations.tsv.gz")
    write_frame(primary, summary_dir / "pc10_conditioning_results.tsv.gz")
    write_frame(panel, summary_dir / "panel_window_robustness.tsv")
    write_frame(consensus, summary_dir / "conditioning_consensus.tsv")
    write_frame(window_summary, summary_dir / "window_scope_summary.tsv")
    write_frame(failures, summary_dir / "failure_rows.tsv.gz")
    payload = {
        "n_pairs": len(pairs), "n_hierarchical_primary": int((pairs.discovery_scope == "hierarchical_primary").sum()),
        "n_global_only_secondary": int((pairs.discovery_scope == "global_only_secondary").sum()),
        "conditioning_class_counts": consensus.conditioning_class.value_counts().to_dict(),
        "n_failure_rows": int((failures.status != "callset_window_unavailable").sum()) if not failures.empty else 0,
        "n_unavailable_rows": int((failures.status == "callset_window_unavailable").sum()) if not failures.empty else 0,
        "panel_interpretation": "The two callsets are sensitivity analyses on the same participants, not independent replication. Structurally unavailable windows are excluded from negative-evidence counts.",
        "support_definition": thresholds,
    }
    (summary_dir / "summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
