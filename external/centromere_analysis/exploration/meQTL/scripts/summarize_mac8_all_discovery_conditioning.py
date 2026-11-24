#!/usr/bin/env python3
"""Validate and summarize source-locus SNP conditioning for all MAC8 q<0.01 pairs."""
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
    frames = [
        pd.read_csv(path, sep="\t", low_memory=False)
        for path in paths
        if path.exists() and path.stat().st_size
    ]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def write_table(frame: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def adjust(frame: pd.DataFrame, p_column: str, q_column: str, group_columns: list[str]) -> None:
    frame[q_column] = np.nan
    for _, indices in frame.groupby(group_columns, sort=False, dropna=False).groups.items():
        index = list(indices)
        values = pd.to_numeric(frame.loc[index, p_column], errors="coerce")
        valid = values.notna() & frame.loc[index, "unconditional_validated"].fillna(False)
        chosen = np.asarray(index)[valid.to_numpy()]
        if len(chosen):
            frame.loc[chosen, q_column] = bh_adjust(values.loc[valid].to_numpy(float))


def main() -> int:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    root = Path(config["workflow_output"])
    summary_dir = root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    pair_root = Path(config.get("pair_manifest_workflow_output", root))
    pairs = pd.read_csv(pair_root / "manifests" / "pairs.tsv.gz", sep="\t")
    input_roots = [Path(value) for value in config.get("input_workflow_outputs", [root])]
    association_frames = []
    variable_frames = []
    failure_frames = []
    for input_root in input_roots:
        input_tasks = pd.read_csv(input_root / "manifests" / "analysis_tasks.tsv", sep="\t")
        input_ids = input_tasks.analysis_task_id.astype(str).tolist()
        association_frames.append(concat([input_root / "analysis" / task_id / "associations.tsv.gz" for task_id in input_ids]))
        variable_frames.append(concat([input_root / "analysis" / task_id / "susie_variables.tsv.gz" for task_id in input_ids]))
        failure_frames.append(concat([input_root / "analysis" / task_id / "failures.tsv" for task_id in input_ids]))
    associations = pd.concat(association_frames, ignore_index=True, sort=False)
    associations = associations.drop_duplicates(["pair_id", "panel", "window_bp", "pc_count"], keep="last")
    variables = pd.concat(variable_frames, ignore_index=True, sort=False)
    if not variables.empty:
        variables = variables.drop_duplicates(["pair_id", "panel", "window_bp", "variable"], keep="last")
    failures = pd.concat(failure_frames, ignore_index=True, sort=False)
    if associations.empty:
        raise ValueError("No conditioning association outputs were found")

    primary = associations.loc[
        pd.to_numeric(associations["pc_count"], errors="coerce").eq(int(config["fixed_pc_count"]))
    ].copy()
    expected = pairs[[
        "pair_id", "scope", "discovery_scope", "expected_rint_p",
        "expected_rint_beta", "expected_raw_p", "expected_raw_beta",
        "target_followup_candidate", "source_target_distance_bp",
    ]]
    primary = primary.merge(expected, on=["pair_id", "discovery_scope"], how="left", validate="many_to_one")
    ok = primary.status.eq("ok")
    primary["unconditional_log10_p_error"] = np.nan
    primary.loc[ok, "unconditional_log10_p_error"] = np.abs(
        np.log10(pd.to_numeric(primary.loc[ok, "rint_unconditional_hc3_p"]))
        - np.log10(pd.to_numeric(primary.loc[ok, "expected_rint_p"]))
    )
    primary["unconditional_beta_error"] = np.nan
    primary.loc[ok, "unconditional_beta_error"] = np.abs(
        pd.to_numeric(primary.loc[ok, "rint_unconditional_beta"])
        - pd.to_numeric(primary.loc[ok, "expected_rint_beta"])
    )
    validation = config["validation"]
    primary["unconditional_validated"] = (
        ok
        & primary.unconditional_log10_p_error.le(float(validation["max_abs_log10_p_error"]))
        & primary.unconditional_beta_error.le(float(validation["max_abs_beta_error"]))
    )

    validation_by_pair = (
        primary.groupby(["pair_id", "scope"], as_index=False)
        .agg(
            n_windows=("analysis_task_id", "size"),
            n_ok=("status", lambda value: int(value.eq("ok").sum())),
            n_validated=("unconditional_validated", "sum"),
            max_abs_log10_p_error=("unconditional_log10_p_error", "max"),
            max_abs_beta_error=("unconditional_beta_error", "max"),
        )
    )
    validation_by_pair["all_ok_windows_validated"] = (
        validation_by_pair.n_ok.gt(0)
        & validation_by_pair.n_ok.eq(validation_by_pair.n_validated)
    )
    write_table(validation_by_pair, summary_dir / "unconditional_validation.tsv")

    if not variables.empty:
        keys = ["pair_id", "panel", "window_bp"]
        sv_pip = variables.loc[
            variables.variable.eq("SV"), [*keys, "pip"]
        ].rename(columns={"pip": "susie_sv_pip"})
        primary = primary.merge(sv_pip, on=keys, how="left", validate="one_to_one")
    else:
        primary["susie_sv_pip"] = np.nan

    for p_column, stem in (
        ("rint_conditional_hc3_p", "conditional"),
        ("crossfit_residual_sv_hc3_p", "crossfit"),
    ):
        adjust(primary, p_column, f"{stem}_bh_q_scope", ["scope", "panel", "window_bp"])
        adjust(primary, p_column, f"{stem}_bh_q_pooled", ["panel", "window_bp"])

    thresholds = config["robustness"]
    sign_agrees = np.sign(primary.rint_conditional_beta).eq(
        np.sign(primary.crossfit_residual_sv_beta)
    )
    common_support = (
        primary.unconditional_validated
        & primary.susie_sv_pip.ge(float(thresholds["susie_sv_pip_min"]))
        & sign_agrees
    )
    primary["window_supported_scope"] = (
        common_support
        & primary.conditional_bh_q_scope.le(float(thresholds["bh_q_max"]))
        & primary.crossfit_bh_q_scope.le(float(thresholds["bh_q_max"]))
    )
    primary["window_supported_pooled"] = (
        common_support
        & primary.conditional_bh_q_pooled.le(float(thresholds["bh_q_max"]))
        & primary.crossfit_bh_q_pooled.le(float(thresholds["bh_q_max"]))
    )

    panel = (
        primary.groupby(["pair_id", "scope", "panel"], as_index=False)
        .agg(
            n_windows=("window_bp", "size"),
            n_ok=("status", lambda value: int(value.eq("ok").sum())),
            n_validated=("unconditional_validated", "sum"),
            n_supported_scope=("window_supported_scope", "sum"),
            n_supported_pooled=("window_supported_pooled", "sum"),
            max_sv_pip=("susie_sv_pip", "max"),
            min_conditional_q_scope=("conditional_bh_q_scope", "min"),
            min_crossfit_q_scope=("crossfit_bh_q_scope", "min"),
            median_raw_effect_retention=("raw_effect_retention", "median"),
        )
    )
    panel["panel_robust_scope"] = panel.n_supported_scope.ge(
        int(thresholds["minimum_supported_windows"])
    )
    panel["panel_robust_pooled"] = panel.n_supported_pooled.ge(
        int(thresholds["minimum_supported_windows"])
    )

    consensus = (
        panel.groupby(["pair_id", "scope"], as_index=False)
        .agg(
            n_robust_panels_scope=("panel_robust_scope", "sum"),
            n_robust_panels_pooled=("panel_robust_pooled", "sum"),
            n_supported_windows_scope=("n_supported_scope", "sum"),
            n_supported_windows_pooled=("n_supported_pooled", "sum"),
            max_sv_pip=("max_sv_pip", "max"),
            min_conditional_q_scope=("min_conditional_q_scope", "min"),
            min_crossfit_q_scope=("min_crossfit_q_scope", "min"),
            median_raw_effect_retention=("median_raw_effect_retention", "median"),
        )
    )
    consensus["conditioning_class_scope"] = np.select(
        [consensus.n_robust_panels_scope.eq(len(config["panels"])), consensus.n_robust_panels_scope.eq(1)],
        ["cross_panel_robust", "single_panel_supported"],
        default="not_robust",
    )
    consensus["conditioning_class_pooled"] = np.select(
        [consensus.n_robust_panels_pooled.eq(len(config["panels"])), consensus.n_robust_panels_pooled.eq(1)],
        ["cross_panel_robust", "single_panel_supported"],
        default="not_robust",
    )
    consensus = pairs.merge(consensus, on=["pair_id", "scope"], how="left", validate="one_to_one")
    target_candidates = consensus.loc[
        consensus.target_followup_candidate
        & consensus.n_supported_windows_scope.fillna(0).gt(0)
    ].copy()

    window_summary = (
        primary.groupby(["scope", "panel", "window_bp"], as_index=False)
        .agg(
            n_pairs=("pair_id", "size"),
            n_ok=("status", lambda value: int(value.eq("ok").sum())),
            n_validated=("unconditional_validated", "sum"),
            n_supported_scope=("window_supported_scope", "sum"),
            n_supported_pooled=("window_supported_pooled", "sum"),
            n_susie_missing=("susie_sv_pip", lambda value: int(value.isna().sum())),
        )
    )

    write_table(primary, summary_dir / "pc10_conditioning_results.tsv.gz")
    write_table(panel, summary_dir / "panel_window_robustness.tsv")
    write_table(consensus, summary_dir / "conditioning_consensus.tsv")
    write_table(target_candidates, summary_dir / "target_followup_candidates.tsv")
    write_table(window_summary, summary_dir / "window_scope_summary.tsv")
    write_table(failures, summary_dir / "failure_rows.tsv.gz")

    invalid = validation_by_pair.loc[~validation_by_pair.all_ok_windows_validated]
    payload = {
        "n_pairs": len(pairs),
        "n_cis_pairs": int(pairs.scope.eq("cis").sum()),
        "n_trans_pairs": int(pairs.scope.eq("trans").sum()),
        "n_pairs_with_validated_window": int(validation_by_pair.n_validated.gt(0).sum()),
        "n_pairs_failing_unconditional_validation": len(invalid),
        "conditioning_class_scope": consensus.conditioning_class_scope.value_counts(dropna=False).to_dict(),
        "conditioning_class_pooled": consensus.conditioning_class_pooled.value_counts(dropna=False).to_dict(),
        "n_target_followup_candidates": len(target_candidates),
        "support_definition": thresholds,
        "q_families": {
            "scope": "all selected pairs within cis or trans, separately by SNP panel/window",
            "pooled": "all 465 selected pairs, separately by SNP panel/window",
        },
        "interpretation": "Post-discovery same-participant sensitivity analysis, not replication.",
    }
    (summary_dir / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    if len(invalid):
        raise ValueError(
            f"{len(invalid)} pairs failed exact unconditional-statistic validation; "
            f"see {summary_dir / 'unconditional_validation.tsv'}"
        )
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
