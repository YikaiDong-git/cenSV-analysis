#!/usr/bin/env python3
"""Condition a batch of at most ten discovered CpGs on one shared SNP matrix."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from comprehensive_conditioning import PC_COUNTS, crossfit_residual_sv, direct_rint, hc3_fit, residualize_matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--snp-matrix", type=Path, required=True)
    parser.add_argument("--snp-status", type=Path, required=True)
    parser.add_argument("--associations-output", type=Path, required=True)
    parser.add_argument("--susie-output", type=Path, required=True)
    parser.add_argument("--failures-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    return parser.parse_args()


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def estimable_condition_indices(covariates: np.ndarray, genotypes: np.ndarray, sv: np.ndarray) -> np.ndarray:
    base = np.column_stack([np.ones(len(sv)), covariates])
    chosen: list[int] = []
    for index in range(genotypes.shape[1]):
        trial = np.column_stack([base, genotypes[:, chosen + [index]], sv])
        if np.linalg.matrix_rank(trial) == trial.shape[1]:
            chosen.append(index)
    return np.asarray(chosen, dtype=int)


def model_row(y: np.ndarray, sv: np.ndarray, covariates: np.ndarray, extras: np.ndarray, transform: str) -> dict[str, object]:
    outcome = direct_rint(y) if transform == "rint" else y
    unconditioned = hc3_fit(outcome, np.column_stack([np.ones(len(y)), covariates, sv]))
    selected = estimable_condition_indices(covariates, extras, sv) if extras.size else np.array([], dtype=int)
    conditional = hc3_fit(outcome, np.column_stack([np.ones(len(y)), covariates, extras[:, selected], sv]))
    return {
        f"{transform}_unconditional_beta": unconditioned["beta"],
        f"{transform}_unconditional_hc3_se": unconditioned["hc3_se"],
        f"{transform}_unconditional_hc3_p": unconditioned["hc3_p"],
        f"{transform}_unconditional_status": unconditioned["status"],
        f"{transform}_unconditional_rank": unconditioned["rank"],
        f"{transform}_unconditional_df": unconditioned["df"],
        f"{transform}_unconditional_max_leverage": unconditioned["max_leverage"],
        f"{transform}_unconditional_mean_leverage": unconditioned["mean_leverage"],
        f"{transform}_unconditional_condition_number": unconditioned["condition_number"],
        f"{transform}_conditional_beta": conditional["beta"],
        f"{transform}_conditional_hc3_se": conditional["hc3_se"],
        f"{transform}_conditional_hc3_t": conditional["hc3_t"],
        f"{transform}_conditional_hc3_p": conditional["hc3_p"],
        f"{transform}_conditional_status": conditional["status"],
        f"{transform}_conditional_rank": conditional["rank"],
        f"{transform}_conditional_df": conditional["df"],
        f"{transform}_conditional_max_leverage": conditional["max_leverage"],
        f"{transform}_conditional_mean_leverage": conditional["mean_leverage"],
        f"{transform}_conditional_condition_number": conditional["condition_number"],
        "n_condition_snps_estimable": len(selected),
    }


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(args.config.read_text())
    root = Path(config["workflow_output"])
    tasks = pd.read_csv(root / "manifests/analysis_tasks.tsv", sep="\t").set_index("analysis_task_id")
    pairs = pd.read_csv(root / "manifests/pairs.tsv.gz", sep="\t").set_index("pair_id")
    task = tasks.loc[args.task_id]
    pair_ids = str(task.pair_ids).split(",")
    if len(pair_ids) > 10 or len(pair_ids) != int(task.n_outcomes):
        raise ValueError("Outcome batch contract violated")
    selected_pairs = pairs.loc[pair_ids].reset_index()
    design = np.load(config["canonical_design"], allow_pickle=True)
    samples = design["sample_ids"].astype(str).tolist()
    if "pc_counts" in design.files:
        cov_counts = dict(zip(design["pc_counts"].astype(int), design["covariate_counts"].astype(int)))
        analysis_pc_counts = PC_COUNTS
    else:
        fixed_pc = int(config.get("fixed_pc_count", 10))
        cov_counts = {fixed_pc: int(design["covariates"].shape[1])}
        analysis_pc_counts = (fixed_pc,)
    all_covariates = np.asarray(design["covariates"], float)
    snp_status = json.loads(args.snp_status.read_text())
    snp_data = np.load(args.snp_matrix, allow_pickle=True)
    if snp_data["sample_ids"].astype(str).tolist() != samples:
        raise ValueError("SNP matrix and canonical design sample order differ")
    genotypes = np.asarray(snp_data["genotypes"], float)
    sv = np.asarray(snp_data["sv_dosage"], float)
    condition_indices = np.asarray(snp_data["condition_indices"], int)
    condition_genotypes = genotypes[:, condition_indices] if len(condition_indices) else np.empty((len(samples), 0))

    outcomes: dict[str, np.ndarray] = {}
    for matrix_path, group in selected_pairs.groupby("matrix_path"):
        wanted = set(group.feature_id.astype(str))
        frame = pd.read_csv(matrix_path, sep="\t", usecols=lambda name: name == "feature_id" or name in samples)
        frame["feature_id"] = frame.feature_id.astype(str)
        frame = frame.loc[frame.feature_id.isin(wanted)].set_index("feature_id")
        for feature_id in wanted:
            if feature_id not in frame.index:
                continue
            outcomes[feature_id] = pd.to_numeric(frame.loc[feature_id, samples], errors="coerce").to_numpy(float)

    associations: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    susie_inputs: list[dict[str, str]] = []
    temp_root = os.environ.get("TMPDIR")
    with tempfile.TemporaryDirectory(prefix="meqtl_condition_", dir=temp_root) as temp_dir:
        temp = Path(temp_dir)
        for pair in selected_pairs.itertuples(index=False):
            base = {"analysis_task_id": args.task_id, "pair_id": pair.pair_id, "sv_id": pair.sv_id,
                    "feature_id": pair.feature_id, "discovery_scope": pair.discovery_scope,
                    "panel": task.panel, "window_bp": int(task.window_bp),
                    "panel_interpretation": "same-participant sensitivity callset; not independent replication"}
            if snp_status.get("status") != "ok":
                failures.append({**base, "stage": "snp_matrix", "status": snp_status.get("status", "unknown"), "message": ""})
                for pc in analysis_pc_counts:
                    associations.append({**base, "pc_count": pc, "status": f"snp_{snp_status.get('status', 'unknown')}"})
                continue
            if str(pair.feature_id) not in outcomes:
                failures.append({**base, "stage": "outcome", "status": "feature_missing", "message": str(pair.matrix_path)})
                for pc in analysis_pc_counts:
                    associations.append({**base, "pc_count": pc, "status": "feature_missing"})
                continue
            y_all = outcomes[str(pair.feature_id)]
            complete = np.isfinite(y_all) & np.isfinite(sv)
            n_complete = int(complete.sum())
            carriers = int(np.sum(sv[complete] > 0))
            if n_complete < int(config["min_complete_samples"]):
                failures.append({**base, "stage": "outcome", "status": "insufficient_complete_cases", "message": str(n_complete)})
                for pc in analysis_pc_counts:
                    associations.append({**base, "pc_count": pc, "status": "insufficient_complete_cases",
                                         "n_complete": n_complete, "n_carriers": carriers})
                continue
            y = y_all[complete]
            sv_complete = sv[complete]
            g_condition = condition_genotypes[complete]
            for pc in analysis_pc_counts:
                covariates = all_covariates[complete, :cov_counts[pc]]
                row = {**base, "pc_count": pc, "status": "ok", "primary_model": pc == 10,
                       "n_complete": n_complete, "n_carriers": carriers, "n_susie_snps": genotypes.shape[1],
                       "n_condition_snps_requested": g_condition.shape[1]}
                row.update(model_row(y, sv_complete, covariates, g_condition, "raw"))
                row.update(model_row(y, sv_complete, covariates, g_condition, "rint"))
                denominator = row["raw_unconditional_beta"]
                row["raw_effect_retention"] = row["raw_conditional_beta"] / denominator if np.isfinite(denominator) and denominator != 0 else np.nan
                if pc == 10:
                    crossfit_residual, crossfit_r2 = crossfit_residual_sv(sv_complete, genotypes[complete], covariates,
                                                                          int(config["ridge_folds"]), int(config["random_seed"]))
                    crossfit_fit = hc3_fit(direct_rint(y), np.column_stack([np.ones(n_complete), covariates, crossfit_residual]))
                    row.update({"crossfit_sv_oof_r2": crossfit_r2, "crossfit_residual_sv_beta": crossfit_fit["beta"],
                                "crossfit_residual_sv_hc3_se": crossfit_fit["hc3_se"], "crossfit_residual_sv_hc3_p": crossfit_fit["hc3_p"],
                                "crossfit_residual_sv_rank": crossfit_fit["rank"], "crossfit_residual_sv_df": crossfit_fit["df"],
                                "crossfit_residual_sv_max_leverage": crossfit_fit["max_leverage"]})
                    cov_resid_y = residualize_matrix(direct_rint(y)[:, None], covariates)[:, 0]
                    predictors = np.column_stack([sv_complete, genotypes[complete]])
                    cov_resid_x = residualize_matrix(predictors, covariates)
                    usable = np.std(cov_resid_x, axis=0) > 1e-10
                    usable[0] = True
                    variable_ids = np.concatenate([["SV"], snp_data["variant_ids"].astype(str)])[usable]
                    susie_frame = pd.DataFrame(cov_resid_x[:, usable], columns=variable_ids)
                    susie_frame.insert(0, "phenotype", cov_resid_y)
                    input_path = temp / f"{pair.pair_id}.tsv"
                    susie_frame.to_csv(input_path, sep="\t", index=False)
                    pair_seed = int(config["random_seed"]) + int(str(pair.pair_id).replace("ccpair", ""))
                    susie_inputs.append({"pair_id": pair.pair_id, "input_path": str(input_path), "seed": pair_seed})
                associations.append(row)

        susie_variables = pd.DataFrame(columns=["pair_id", "variable", "pip", "credible_sets"])
        if susie_inputs:
            input_manifest = temp / "susie_manifest.tsv"
            variable_output = temp / "susie_variables.tsv"
            status_output = temp / "susie_status.tsv"
            pd.DataFrame(susie_inputs).to_csv(input_manifest, sep="\t", index=False)
            command = shlex.split(str(config["rscript_cmd"])) + ["scripts/run_comprehensive_susie_batch.R", str(input_manifest), str(variable_output), str(status_output)]
            result = subprocess.run(command, text=True, capture_output=True)
            if result.returncode == 0:
                susie_variables = pd.read_csv(variable_output, sep="\t")
                susie_variables["analysis_task_id"] = args.task_id
                susie_variables["panel"] = task.panel
                susie_variables["window_bp"] = int(task.window_bp)
                susie_status = pd.read_csv(status_output, sep="\t")
                successful_pairs = set(susie_status.loc[susie_status.status.eq("ok"), "pair_id"].astype(str))
                susie_variables = susie_variables.loc[susie_variables.pair_id.astype(str).isin(successful_pairs)].copy()
                for row in susie_status.itertuples(index=False):
                    if row.status != "ok":
                        pair = selected_pairs.set_index("pair_id").loc[str(row.pair_id)]
                        failures.append({"analysis_task_id": args.task_id, "pair_id": row.pair_id, "sv_id": pair.sv_id,
                                         "feature_id": pair.feature_id, "discovery_scope": pair.discovery_scope,
                                         "panel": task.panel, "window_bp": int(task.window_bp), "stage": "susie",
                                         "status": row.status, "message": row.message})
            else:
                for item in susie_inputs:
                    pair = selected_pairs.set_index("pair_id").loc[item["pair_id"]]
                    failures.append({"analysis_task_id": args.task_id, "pair_id": item["pair_id"], "sv_id": pair.sv_id,
                                     "feature_id": pair.feature_id, "discovery_scope": pair.discovery_scope,
                                     "panel": task.panel, "window_bp": int(task.window_bp), "stage": "susie",
                                     "status": "batch_runner_failed", "message": result.stderr[-2000:]})

    association_frame = pd.DataFrame(associations)
    write_frame(association_frame, args.associations_output)
    write_frame(susie_variables, args.susie_output)
    write_frame(pd.DataFrame(failures, columns=["analysis_task_id", "pair_id", "sv_id", "feature_id", "discovery_scope", "panel", "window_bp", "stage", "status", "message"]), args.failures_output)
    summary = {"analysis_task_id": args.task_id, "snp_task_id": task.snp_task_id, "snp_status": snp_status.get("status"),
               "n_pairs": len(pair_ids), "n_association_rows": len(association_frame), "n_susie_rows": len(susie_variables),
               "n_failure_rows": len(failures), "node_local_susie_inputs_removed": True}
    args.summary_output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
