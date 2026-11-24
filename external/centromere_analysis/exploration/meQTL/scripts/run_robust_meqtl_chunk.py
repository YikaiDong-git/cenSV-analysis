#!/usr/bin/env python3
"""Run one robust observed-SV meQTL chunk with shared multiplier operators."""

from __future__ import annotations

import argparse
import heapq
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from robust_meqtl import (
    METHOD_VERSION,
    array_fingerprint,
    batched_reduced_residuals,
    complete_case_dosage_frequency,
    design_fingerprint,
    fwl_hc3_statistics,
    hc3_null_score_influence,
    multiplier_score_statistics,
    residualize_outcome_batched,
    residualize_predictor_batched,
    row_rank_inverse_normal,
    sample_order_fingerprint,
)
from run_sv_methylation_assoc import MATRIX_META_COLUMNS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-dosages", type=Path, required=True)
    parser.add_argument("--sv-metadata", type=Path, required=True)
    parser.add_argument("--sv-indices", required=True)
    parser.add_argument("--sv-qc", type=Path, required=True)
    parser.add_argument("--methylation-matrix", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--targets", type=Path)
    parser.add_argument("--min-carriers", type=int, default=0)
    parser.add_argument("--primary-min-carriers", type=int, default=0)
    parser.add_argument("--min-mac", type=float, default=8.0)
    parser.add_argument("--primary-min-mac", type=float, default=8.0)
    parser.add_argument("--min-feature-fraction", type=float, default=0.9)
    parser.add_argument("--multiplier-batch-size", type=int, default=250)
    parser.add_argument("--feature-batch-size", type=int, default=500)
    parser.add_argument("--score-thresholds", default="2,2.5,3,3.5,4,4.5,5,6,7,8")
    parser.add_argument("--pvalue-threshold", type=float, default=1e-4)
    parser.add_argument("--top-per-sv", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pvalues-output", type=Path, required=True)
    parser.add_argument("--null-output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def bool_series(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def scalar_string(value: np.ndarray) -> str:
    return str(np.asarray(value).item())


def load_design(path: Path) -> dict[str, object]:
    with np.load(path) as design:
        loaded = {name: design[name].copy() for name in design.files}
    samples = loaded["sample_ids"].astype(str).tolist()
    covariate_names = loaded["covariate_names"].astype(str).tolist()
    if scalar_string(loaded["method_version"]) != METHOD_VERSION:
        raise ValueError("Shared design method version is incompatible")
    if scalar_string(loaded["sample_order_fingerprint"]) != sample_order_fingerprint(samples):
        raise ValueError("Shared design sample-order fingerprint is invalid")
    expected = design_fingerprint(
        samples,
        covariate_names,
        loaded["covariates"],
        loaded["replicate_seeds"],
        loaded["mammen"],
        loaded["rademacher"],
    )
    if scalar_string(loaded["design_fingerprint"]) != expected:
        raise ValueError("Shared design fingerprint is invalid")
    mammen_fp = array_fingerprint(loaded["mammen"], labels=[*samples, "mammen"])
    rademacher_fp = array_fingerprint(
        loaded["rademacher"], labels=[*samples, "rademacher"]
    )
    if scalar_string(loaded["mammen_fingerprint"]) != mammen_fp:
        raise ValueError("Mammen operator fingerprint is invalid")
    if scalar_string(loaded["rademacher_fingerprint"]) != rademacher_fp:
        raise ValueError("Rademacher operator fingerprint is invalid")
    return loaded


def empty_outputs(
    args: argparse.Namespace,
    design: dict[str, object],
    metadata: pd.DataFrame,
    n_features_input: int,
    n_features_eligible: int,
    stale_count_rows: int,
    thresholds: np.ndarray,
) -> int:
    n_replicates = len(design["replicate_seeds"])
    sv_ids = metadata["id"].to_numpy(dtype=str) if not metadata.empty else np.asarray([], dtype=str)
    pd.DataFrame().to_csv(args.output, sep="\t", index=False, compression="gzip")
    np.savez_compressed(
        args.pvalues_output,
        sample_ids=design["sample_ids"],
        sample_order_fingerprint=design["sample_order_fingerprint"],
        design_fingerprint=design["design_fingerprint"],
        mammen_fingerprint=design["mammen_fingerprint"],
        rademacher_fingerprint=design["rademacher_fingerprint"],
        sv_ids=sv_ids,
        sv_position=np.asarray([], dtype=np.int16),
        feature_index=np.asarray([], dtype=np.int32),
        primary_family=np.asarray([], dtype=bool),
        rint_hc3_t=np.asarray([], dtype=np.float32),
        rint_hc3_p_t=np.asarray([], dtype=np.float64),
        rint_hc3_p_asymptotic=np.asarray([], dtype=np.float64),
        rint_score=np.asarray([], dtype=np.float32),
    )
    np.savez_compressed(
        args.null_output,
        sample_ids=design["sample_ids"],
        sample_order_fingerprint=design["sample_order_fingerprint"],
        design_fingerprint=design["design_fingerprint"],
        mammen_fingerprint=design["mammen_fingerprint"],
        rademacher_fingerprint=design["rademacher_fingerprint"],
        replicate_seeds=design["replicate_seeds"],
        sv_ids=sv_ids,
        score_thresholds=thresholds,
        primary_mammen_maxima=np.zeros((n_replicates, len(sv_ids)), dtype=np.float32),
        sensitivity_rademacher_maxima=np.zeros((n_replicates, len(sv_ids)), dtype=np.float32),
        primary_mammen_exceedances=np.zeros((n_replicates, len(sv_ids), len(thresholds)), dtype=np.int32),
        sensitivity_rademacher_exceedances=np.zeros((n_replicates, len(sv_ids), len(thresholds)), dtype=np.int32),
        target_keys=np.asarray([], dtype=str),
        target_rademacher_statistics=np.empty((n_replicates, 0), dtype=np.float32),
    )
    summary = {
        "method_version": METHOD_VERSION,
        "n_samples": len(design["sample_ids"]),
        "sample_ids": design["sample_ids"].astype(str).tolist(),
        "sample_order_fingerprint": scalar_string(design["sample_order_fingerprint"]),
        "design_fingerprint": scalar_string(design["design_fingerprint"]),
        "mammen_fingerprint": scalar_string(design["mammen_fingerprint"]),
        "rademacher_fingerprint": scalar_string(design["rademacher_fingerprint"]),
        "n_replicates": n_replicates,
        "n_svs_upstream_qc": len(metadata),
        "n_features_input": n_features_input,
        "n_features_actual_ge_min_fraction": n_features_eligible,
        "minimum_actual_observed": math.ceil(args.min_feature_fraction * len(design["sample_ids"])),
        "n_rows_stale_source_n_observed": stale_count_rows,
        "min_carriers": args.min_carriers,
        "primary_min_carriers": args.primary_min_carriers,
        "min_mac": args.min_mac,
        "primary_min_mac": args.primary_min_mac,
        "min_feature_fraction": args.min_feature_fraction,
        "pvalue_threshold": args.pvalue_threshold,
        "score_thresholds": thresholds.tolist(),
        "n_sensitivity_tests": 0,
        "n_primary_tests": 0,
        "n_retained": 0,
        "n_targets": 0,
        "observed_primary_score_exceedances": [0] * len(thresholds),
        "observed_sensitivity_score_exceedances": [0] * len(thresholds),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


def main() -> int:
    args = parse_args()
    if not (0.0 < args.min_feature_fraction <= 1.0):
        raise ValueError("--min-feature-fraction must be in (0, 1]")
    if args.min_mac < 0.0 or args.primary_min_mac < args.min_mac:
        raise ValueError("MAC thresholds must satisfy 0 <= min <= primary")
    if not (0.0 < args.pvalue_threshold <= 1.0):
        raise ValueError("--pvalue-threshold must be in (0, 1]")
    thresholds = np.asarray(
        [float(value) for value in args.score_thresholds.split(",") if value], dtype=float
    )
    if len(thresholds) == 0 or np.any(thresholds <= 0) or np.any(np.diff(thresholds) <= 0):
        raise ValueError("--score-thresholds must be positive and strictly increasing")

    design = load_design(args.design)
    samples = design["sample_ids"].astype(str).tolist()
    covariates = np.asarray(design["covariates"], dtype=float)
    mammen = np.asarray(design["mammen"], dtype=np.float32)
    rademacher = np.asarray(design["rademacher"], dtype=np.int8)
    n_replicates = len(mammen)
    if len(rademacher) != n_replicates:
        raise ValueError("Multiplier operators have different replicate counts")

    sv_indices = [int(value) for value in args.sv_indices.split(",") if value]
    dosage = pd.read_csv(args.sv_dosages, sep="\t")
    metadata = pd.read_csv(args.sv_metadata, sep="\t")
    qc = pd.read_csv(args.sv_qc, sep="\t")
    qc["qc_pass"] = bool_series(qc["qc_pass"])
    if not set(samples).issubset(dosage.columns):
        missing = sorted(set(samples) - set(dosage.columns))
        raise ValueError(f"Dosage input lacks canonical samples: {missing}")

    dosage = dosage[dosage["sv_index"].isin(sv_indices)].copy()
    metadata = metadata[metadata["sv_index"].isin(sv_indices)].copy()
    metadata = metadata.merge(
        qc[["sv_id", "qc_pass", "n_carriers", "call_rate", "dosage_r"]].rename(
            columns={
                "n_carriers": "upstream_n_carriers",
                "call_rate": "upstream_call_rate",
                "dosage_r": "upstream_dosage_r",
            }
        ),
        left_on="id",
        right_on="sv_id",
        how="left",
        validate="one_to_one",
    )
    metadata = metadata[metadata["qc_pass"].fillna(False)].copy()
    dosage = dosage[dosage["sv_index"].isin(metadata["sv_index"])].copy()
    metadata = metadata.set_index("sv_index").loc[dosage["sv_index"]].reset_index()
    x_values = dosage[samples].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    canonical_masks = np.isfinite(x_values)
    canonical_mac = np.asarray(
        [complete_case_dosage_frequency(row, mask[None, :])[0][0]
         for row, mask in zip(x_values, canonical_masks, strict=True)],
        dtype=float,
    )
    canonical_n = canonical_masks.sum(axis=1)
    canonical_maf = canonical_mac / (2.0 * canonical_n)
    metadata["canonical_dosage_mac"] = canonical_mac
    metadata["canonical_dosage_maf"] = canonical_maf
    profile_keep = canonical_mac >= args.min_mac
    metadata = metadata.loc[profile_keep].reset_index(drop=True)
    x_values = x_values[profile_keep]

    methylation = pd.read_csv(args.methylation_matrix, sep="\t")
    if not set(samples).issubset(methylation.columns):
        missing = sorted(set(samples) - set(methylation.columns))
        raise ValueError(f"Methylation input lacks canonical samples: {missing}")
    n_features_input = len(methylation)
    values = methylation[samples].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    actual_n_observed = np.isfinite(values).sum(axis=1)
    minimum_observed = math.ceil(args.min_feature_fraction * len(samples))
    source_n_observed = pd.to_numeric(methylation.get("n_observed"), errors="coerce")
    stale_count_rows = int(
        np.count_nonzero(
            np.isfinite(source_n_observed.to_numpy(float))
            & (source_n_observed.to_numpy(float) != actual_n_observed)
        )
    )
    feature_keep = actual_n_observed >= minimum_observed
    methylation = methylation.loc[feature_keep].reset_index(drop=True)
    values = values[feature_keep]
    actual_n_observed = actual_n_observed[feature_keep]

    for path in (args.output, args.pvalues_output, args.null_output, args.summary):
        path.parent.mkdir(parents=True, exist_ok=True)
    if metadata.empty or methylation.empty:
        return empty_outputs(
            args,
            design,
            metadata,
            n_features_input,
            len(methylation),
            stale_count_rows,
            thresholds,
        )

    target_lookup: dict[str, set[str]] = {}
    if args.targets and args.targets.exists():
        targets = pd.read_csv(args.targets, sep="\t")
        for sv_id, group in targets.groupby("sv_id"):
            target_lookup[str(sv_id)] = set(group["feature_id"].astype(str))
    feature_ids = methylation["feature_id"].astype(str).to_numpy()

    n_svs = len(metadata)
    primary_maxima = np.zeros((n_replicates, n_svs), dtype=np.float32)
    sensitivity_maxima = np.zeros((n_replicates, n_svs), dtype=np.float32)
    primary_exceed = np.zeros(
        (n_replicates, n_svs, len(thresholds)), dtype=np.int32
    )
    sensitivity_exceed = np.zeros_like(primary_exceed)
    compact: dict[str, list[np.ndarray]] = defaultdict(list)
    retained: dict[tuple[int, int], dict[str, object]] = {}
    top_heaps: dict[int, list[tuple[float, int, dict[str, object]]]] = defaultdict(list)
    target_keys: list[str] = []
    target_statistics: list[np.ndarray] = []
    test_offset = 0
    heap_counter = 0
    observed_score_counts_primary = np.zeros(len(thresholds), dtype=np.int64)
    observed_score_counts_sensitivity = np.zeros(len(thresholds), dtype=np.int64)

    def row_for(
        sv_position: int,
        feature_index: int,
        compact_index: int,
        n_analyzed: int,
        n_carriers: int,
        dosage_mac: float,
        dosage_maf: float,
        residual_df: int,
        primary: bool,
        raw_beta: float,
        raw_se: float,
        raw_t: float,
        rint_beta: float,
        rint_se: float,
        rint_t: float,
        p_t: float,
        p_asymptotic: float,
        score: float,
    ) -> dict[str, object]:
        sv = metadata.iloc[sv_position]
        feature = methylation.iloc[feature_index]
        return {
            "_feature_index": feature_index,
            "compact_test_index": compact_index,
            "sv_index": int(sv["sv_index"]),
            "sv_id": str(sv["id"]),
            "sv_chrom": sv["chrom"],
            "sv_pos": sv["pos"],
            "sv_end": sv["end"],
            "svtype": sv["svtype"],
            "svlen": sv["svlen"],
            "upstream_n_carriers": sv["upstream_n_carriers"],
            "upstream_call_rate": sv["upstream_call_rate"],
            "upstream_dosage_r": sv["upstream_dosage_r"],
            "n_samples_canonical": len(samples),
            "feature_n_observed_actual": int(actual_n_observed[feature_index]),
            "n_samples_analyzed": n_analyzed,
            "n_carriers_analyzed": n_carriers,
            "dosage_mac_analyzed": dosage_mac,
            "dosage_maf_analyzed": dosage_maf,
            "canonical_dosage_mac": float(sv["canonical_dosage_mac"]),
            "canonical_dosage_maf": float(sv["canonical_dosage_maf"]),
            "primary_family": primary,
            "sensitivity_family": True,
            "residual_df": residual_df,
            "feature_id": str(feature["feature_id"]),
            "methyl_chrom": feature["chrom"],
            "methyl_start": feature["start"],
            "methyl_end": feature["end"],
            "raw_beta": raw_beta,
            "raw_hc3_se": raw_se,
            "raw_hc3_t": raw_t,
            "raw_hc3_p_t": 2.0 * stats.t.sf(abs(raw_t), df=residual_df),
            "rint_beta": rint_beta,
            "rint_hc3_se": rint_se,
            "rint_hc3_t": rint_t,
            "rint_hc3_p_t": p_t,
            "rint_hc3_p_asymptotic": p_asymptotic,
            "rint_hc3_null_score": score,
            "retained_reason": set(),
        }

    genotype_groups: dict[bytes, list[int]] = defaultdict(list)
    for sv_position in range(n_svs):
        genotype_groups[np.packbits(np.isfinite(x_values[sv_position])).tobytes()].append(
            sv_position
        )

    for feature_start in range(0, len(methylation), args.feature_batch_size):
        feature_stop = min(feature_start + args.feature_batch_size, len(methylation))
        raw_block = values[feature_start:feature_stop]
        outcome_masks = np.isfinite(raw_block)
        for sv_positions in genotype_groups.values():
            genotype_mask = np.isfinite(x_values[sv_positions[0]])
            masks = outcome_masks & genotype_mask[None, :]
            rint_block = row_rank_inverse_normal(np.where(masks, raw_block, np.nan))
            rint_residual, leverage, bread, rank, reduced_design = batched_reduced_residuals(
                rint_block, covariates, masks
            )
            raw_residual = residualize_outcome_batched(
                raw_block, masks, bread, reduced_design
            )
            n_analyzed = masks.sum(axis=1).astype(np.int16)
            residual_df = n_analyzed - rank - 1

            for sv_position in sv_positions:
                x = x_values[sv_position]
                x_residual = residualize_predictor_batched(
                    x, masks, bread, reduced_design
                )
                carrier = x > 0.0
                n_carriers = np.sum(masks & carrier[None, :], axis=1).astype(np.int16)
                dosage_mac, dosage_maf = complete_case_dosage_frequency(x, masks)
                x_min = np.min(np.where(masks, x[None, :], np.inf), axis=1)
                x_max = np.max(np.where(masks, x[None, :], -np.inf), axis=1)
                nonmonomorphic = (x_max - x_min) > 1e-12
                raw_beta, raw_se, raw_t = fwl_hc3_statistics(
                    x_residual, raw_residual, leverage
                )
                rint_beta, rint_se, rint_t = fwl_hc3_statistics(
                    x_residual, rint_residual, leverage
                )
                influence, score = hc3_null_score_influence(
                    x_residual, rint_residual, leverage
                )
                with np.errstate(invalid="ignore"):
                    p_t = 2.0 * stats.t.sf(np.abs(rint_t), df=residual_df)
                    p_asymptotic = 2.0 * stats.norm.sf(np.abs(rint_t))
                sensitivity = (
                    (n_carriers >= args.min_carriers)
                    & (dosage_mac >= args.min_mac)
                    & nonmonomorphic
                    & (residual_df > 0)
                    & np.isfinite(p_t)
                    & np.isfinite(score)
                )
                primary = sensitivity & (
                    (n_carriers >= args.primary_min_carriers)
                    & (dosage_mac >= args.primary_min_mac)
                )
                valid_local = np.flatnonzero(sensitivity)
                if not len(valid_local):
                    continue

                n_valid = len(valid_local)
                compact_indices = np.arange(
                    test_offset, test_offset + n_valid, dtype=np.int64
                )
                compact["sv_position"].append(
                    np.full(n_valid, sv_position, dtype=np.int16)
                )
                compact["feature_index"].append(
                    (feature_start + valid_local).astype(np.int32)
                )
                compact["primary_family"].append(primary[valid_local])
                compact["rint_hc3_t"].append(rint_t[valid_local].astype(np.float32))
                compact["rint_hc3_p_t"].append(p_t[valid_local].astype(np.float64))
                compact["rint_hc3_p_asymptotic"].append(
                    p_asymptotic[valid_local].astype(np.float64)
                )
                compact["rint_score"].append(score[valid_local].astype(np.float32))
                test_offset += n_valid

                abs_score = np.abs(score[valid_local])
                for threshold_index, threshold in enumerate(thresholds):
                    observed_score_counts_sensitivity[threshold_index] += np.count_nonzero(
                        abs_score >= threshold
                    )
                    observed_score_counts_primary[threshold_index] += np.count_nonzero(
                        np.abs(score[valid_local[primary[valid_local]]]) >= threshold
                    )

                influence_valid = influence[valid_local]
                primary_valid = primary[valid_local]
                target_set = target_lookup.get(str(metadata.iloc[sv_position]["id"]), set())
                target_local = np.asarray(
                    [
                        index
                        for index, local in enumerate(valid_local)
                        if feature_ids[feature_start + local] in target_set
                    ],
                    dtype=int,
                )
                target_buffer = (
                    np.empty((n_replicates, len(target_local)), dtype=np.float32)
                    if len(target_local)
                    else None
                )
                for replicate_start in range(
                    0, n_replicates, args.multiplier_batch_size
                ):
                    replicate_stop = min(
                        replicate_start + args.multiplier_batch_size, n_replicates
                    )
                    sensitivity_stats = np.abs(
                        multiplier_score_statistics(
                            rademacher[replicate_start:replicate_stop],
                            influence_valid,
                        )
                    )
                    sensitivity_maxima[replicate_start:replicate_stop, sv_position] = np.maximum(
                        sensitivity_maxima[replicate_start:replicate_stop, sv_position],
                        np.max(sensitivity_stats, axis=1),
                    )
                    for threshold_index, threshold in enumerate(thresholds):
                        sensitivity_exceed[
                            replicate_start:replicate_stop, sv_position, threshold_index
                        ] += np.sum(sensitivity_stats >= threshold, axis=1)
                    if target_buffer is not None:
                        target_buffer[replicate_start:replicate_stop] = sensitivity_stats[
                            :, target_local
                        ]
                    if np.any(primary_valid):
                        primary_stats = np.abs(
                            multiplier_score_statistics(
                                mammen[replicate_start:replicate_stop],
                                influence_valid[primary_valid],
                            )
                        )
                        primary_maxima[
                            replicate_start:replicate_stop, sv_position
                        ] = np.maximum(
                            primary_maxima[
                                replicate_start:replicate_stop, sv_position
                            ],
                            np.max(primary_stats, axis=1),
                        )
                        for threshold_index, threshold in enumerate(thresholds):
                            primary_exceed[
                                replicate_start:replicate_stop,
                                sv_position,
                                threshold_index,
                            ] += np.sum(primary_stats >= threshold, axis=1)

                if target_buffer is not None:
                    for target_position, local_position in enumerate(target_local):
                        local = valid_local[local_position]
                        target_keys.append(
                            f"{metadata.iloc[sv_position]['id']}|{feature_ids[feature_start + local]}"
                        )
                        target_statistics.append(target_buffer[:, target_position].copy())

                for valid_position, local in enumerate(valid_local):
                    feature_index = feature_start + local
                    row = row_for(
                        sv_position,
                        feature_index,
                        int(compact_indices[valid_position]),
                        int(n_analyzed[local]),
                        int(n_carriers[local]),
                        float(dosage_mac[local]),
                        float(dosage_maf[local]),
                        int(residual_df[local]),
                        bool(primary[local]),
                        float(raw_beta[local]),
                        float(raw_se[local]),
                        float(raw_t[local]),
                        float(rint_beta[local]),
                        float(rint_se[local]),
                        float(rint_t[local]),
                        float(p_t[local]),
                        float(p_asymptotic[local]),
                        float(score[local]),
                    )
                    key = (sv_position, feature_index)
                    if p_t[local] <= args.pvalue_threshold:
                        row["retained_reason"].add("threshold")
                        retained[key] = row
                    if feature_ids[feature_index] in target_lookup.get(
                        str(metadata.iloc[sv_position]["id"]), set()
                    ):
                        if key in retained:
                            retained[key]["retained_reason"].add("target")
                        else:
                            row["retained_reason"].add("target")
                            retained[key] = row
                    heap_counter += 1
                    heap = top_heaps[sv_position]
                    heapq.heappush(heap, (-float(p_t[local]), heap_counter, row))
                    if len(heap) > args.top_per_sv:
                        heapq.heappop(heap)

    for sv_position, heap in top_heaps.items():
        for _, _, row in heap:
            key = (sv_position, int(row["_feature_index"]))
            if key in retained:
                retained[key]["retained_reason"].add("top")
            else:
                row["retained_reason"].add("top")
                retained[key] = row

    rows = list(retained.values())
    for row in rows:
        row.pop("_feature_index")
        row["retained_reason"] = ",".join(sorted(row["retained_reason"]))
    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("rint_hc3_p_t")
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")

    compact_arrays = {
        name: np.concatenate(parts) if parts else np.asarray([])
        for name, parts in compact.items()
    }
    np.savez_compressed(
        args.pvalues_output,
        sample_ids=design["sample_ids"],
        sample_order_fingerprint=design["sample_order_fingerprint"],
        design_fingerprint=design["design_fingerprint"],
        mammen_fingerprint=design["mammen_fingerprint"],
        rademacher_fingerprint=design["rademacher_fingerprint"],
        sv_ids=metadata["id"].to_numpy(dtype=str),
        **compact_arrays,
    )
    target_matrix = (
        np.column_stack(target_statistics)
        if target_statistics
        else np.empty((n_replicates, 0), dtype=np.float32)
    )
    np.savez_compressed(
        args.null_output,
        sample_ids=design["sample_ids"],
        sample_order_fingerprint=design["sample_order_fingerprint"],
        design_fingerprint=design["design_fingerprint"],
        mammen_fingerprint=design["mammen_fingerprint"],
        rademacher_fingerprint=design["rademacher_fingerprint"],
        replicate_seeds=design["replicate_seeds"],
        sv_ids=metadata["id"].to_numpy(dtype=str),
        score_thresholds=thresholds,
        primary_mammen_maxima=primary_maxima,
        sensitivity_rademacher_maxima=sensitivity_maxima,
        primary_mammen_exceedances=primary_exceed,
        sensitivity_rademacher_exceedances=sensitivity_exceed,
        target_keys=np.asarray(target_keys),
        target_rademacher_statistics=target_matrix,
    )
    summary = {
        "method_version": METHOD_VERSION,
        "n_samples": len(samples),
        "sample_ids": samples,
        "sample_order_fingerprint": scalar_string(design["sample_order_fingerprint"]),
        "design_fingerprint": scalar_string(design["design_fingerprint"]),
        "mammen_fingerprint": scalar_string(design["mammen_fingerprint"]),
        "rademacher_fingerprint": scalar_string(design["rademacher_fingerprint"]),
        "n_replicates": n_replicates,
        "n_svs_upstream_qc": n_svs,
        "n_features_input": n_features_input,
        "n_features_actual_ge_min_fraction": len(methylation),
        "minimum_actual_observed": minimum_observed,
        "n_rows_stale_source_n_observed": stale_count_rows,
        "n_sensitivity_tests": int(len(compact_arrays["rint_hc3_p_t"])),
        "n_primary_tests": int(np.sum(compact_arrays["primary_family"])),
        "n_retained": len(result),
        "n_targets": len(target_keys),
        "min_carriers": args.min_carriers,
        "primary_min_carriers": args.primary_min_carriers,
        "min_mac": args.min_mac,
        "primary_min_mac": args.primary_min_mac,
        "min_feature_fraction": args.min_feature_fraction,
        "pvalue_threshold": args.pvalue_threshold,
        "score_thresholds": thresholds.tolist(),
        "observed_primary_score_exceedances": observed_score_counts_primary.tolist(),
        "observed_sensitivity_score_exceedances": observed_score_counts_sensitivity.tolist(),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
