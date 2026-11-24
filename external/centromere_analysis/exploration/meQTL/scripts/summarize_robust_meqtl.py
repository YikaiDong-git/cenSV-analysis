#!/usr/bin/env python3
"""Aggregate robust meQTL chunks with exact BH and multiplier calibration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from robust_meqtl import METHOD_VERSION, bh_adjust


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-manifest", type=Path, required=True)
    parser.add_argument("--jobs-dir", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-carriers", type=int, required=True)
    parser.add_argument("--primary-min-carriers", type=int, required=True)
    parser.add_argument("--min-mac", type=float, required=True)
    parser.add_argument("--primary-min-mac", type=float, required=True)
    parser.add_argument("--min-feature-fraction", type=float, required=True)
    parser.add_argument("--pvalue-threshold", type=float, required=True)
    return parser.parse_args()


def scalar_string(value: np.ndarray) -> str:
    return str(np.asarray(value).item())


def empirical_p(null: np.ndarray, observed: float) -> float:
    if not len(null) or not np.isfinite(observed):
        return float("nan")
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def read_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except (EmptyDataError, FileNotFoundError):
        return pd.DataFrame()


def assert_provenance(
    loaded: object,
    path: Path,
    expected: dict[str, object],
) -> None:
    samples = loaded["sample_ids"].astype(str)
    if not np.array_equal(samples, expected["sample_ids"]):
        raise ValueError(f"Canonical sample IDs/order differ in {path}")
    for field in (
        "sample_order_fingerprint",
        "design_fingerprint",
        "mammen_fingerprint",
        "rademacher_fingerprint",
    ):
        if scalar_string(loaded[field]) != expected[field]:
            raise ValueError(f"{field} differs in {path}")


def main() -> int:
    args = parse_args()
    manifest = pd.read_csv(args.job_manifest, sep="\t")
    with np.load(args.design) as design:
        expected = {
            "sample_ids": design["sample_ids"].astype(str).copy(),
            "sample_order_fingerprint": scalar_string(design["sample_order_fingerprint"]),
            "design_fingerprint": scalar_string(design["design_fingerprint"]),
            "mammen_fingerprint": scalar_string(design["mammen_fingerprint"]),
            "rademacher_fingerprint": scalar_string(
                design["rademacher_fingerprint"]
            ),
            "replicate_seeds": design["replicate_seeds"].copy(),
        }
    expected_summary_parameters = {
        "min_carriers": args.min_carriers,
        "primary_min_carriers": args.primary_min_carriers,
        "min_mac": args.min_mac,
        "primary_min_mac": args.primary_min_mac,
        "min_feature_fraction": args.min_feature_fraction,
        "pvalue_threshold": args.pvalue_threshold,
        "n_replicates": len(expected["replicate_seeds"]),
    }

    frames: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []
    p_t_parts: list[np.ndarray] = []
    p_asymptotic_parts: list[np.ndarray] = []
    score_parts: list[np.ndarray] = []
    primary_parts: list[np.ndarray] = []
    job_offsets = [0]
    pvalue_files: list[tuple[str, int, int]] = []
    primary_maxima_by_sv: dict[str, np.ndarray] = {}
    sensitivity_maxima_by_sv: dict[str, np.ndarray] = {}
    primary_exceed_by_sv: dict[str, np.ndarray] = {}
    sensitivity_exceed_by_sv: dict[str, np.ndarray] = {}
    target_null: dict[str, np.ndarray] = {}
    thresholds: np.ndarray | None = None

    for row in manifest.itertuples(index=False):
        job_dir = args.jobs_dir / row.job_id
        frame = read_table(job_dir / "observed.tsv.gz")
        if not frame.empty:
            frame.insert(0, "job_id", row.job_id)
            frames.append(frame)
        summary_path = job_dir / "summary.json"
        if not summary_path.exists():
            raise FileNotFoundError(summary_path)
        summary = json.loads(summary_path.read_text())
        if summary.get("method_version") != METHOD_VERSION:
            raise ValueError(f"Method version differs in {summary_path}")
        if summary.get("sample_ids") != expected["sample_ids"].tolist():
            raise ValueError(f"Canonical sample IDs/order differ in {summary_path}")
        for field in (
            "sample_order_fingerprint",
            "design_fingerprint",
            "mammen_fingerprint",
            "rademacher_fingerprint",
        ):
            if summary.get(field) != expected[field]:
                raise ValueError(f"{field} differs in {summary_path}")
        for field, value in expected_summary_parameters.items():
            if summary.get(field) != value:
                raise ValueError(f"{field} differs in {summary_path}")
        summary["job_id"] = row.job_id
        summaries.append(summary)

        pvalue_path = job_dir / "all_pvalues.npz"
        if not pvalue_path.exists():
            raise FileNotFoundError(pvalue_path)
        with np.load(pvalue_path) as values:
            assert_provenance(values, pvalue_path, expected)
            p_t = values["rint_hc3_p_t"].astype(np.float64)
            p_asymptotic = values["rint_hc3_p_asymptotic"].astype(np.float64)
            score = values["rint_score"].astype(np.float64)
            primary = values["primary_family"].astype(bool)
        if not (
            len(p_t) == len(p_asymptotic) == len(score) == len(primary)
        ):
            raise ValueError(f"Compact p-value arrays differ in {pvalue_path}")
        start = job_offsets[-1]
        stop = start + len(p_t)
        pvalue_files.append((row.job_id, start, stop))
        job_offsets.append(stop)
        p_t_parts.append(p_t)
        p_asymptotic_parts.append(p_asymptotic)
        score_parts.append(score)
        primary_parts.append(primary)

        null_path = job_dir / "null_maxima.npz"
        if not null_path.exists():
            raise FileNotFoundError(null_path)
        with np.load(null_path) as null:
            assert_provenance(null, null_path, expected)
            if not np.array_equal(
                null["replicate_seeds"], expected["replicate_seeds"]
            ):
                raise ValueError(f"Replicate seeds differ in {null_path}")
            chunk_thresholds = null["score_thresholds"].astype(float)
            if thresholds is None:
                thresholds = chunk_thresholds.copy()
            elif not np.array_equal(thresholds, chunk_thresholds):
                raise ValueError(f"Score thresholds differ in {null_path}")
            sv_ids = null["sv_ids"].astype(str)
            primary_maxima = null["primary_mammen_maxima"]
            sensitivity_maxima = null["sensitivity_rademacher_maxima"]
            primary_exceed = null["primary_mammen_exceedances"]
            sensitivity_exceed = null["sensitivity_rademacher_exceedances"]
            for index, sv_id in enumerate(sv_ids):
                if sv_id in primary_maxima_by_sv:
                    primary_maxima_by_sv[sv_id] = np.maximum(
                        primary_maxima_by_sv[sv_id], primary_maxima[:, index]
                    )
                    sensitivity_maxima_by_sv[sv_id] = np.maximum(
                        sensitivity_maxima_by_sv[sv_id],
                        sensitivity_maxima[:, index],
                    )
                    primary_exceed_by_sv[sv_id] += primary_exceed[:, index]
                    sensitivity_exceed_by_sv[sv_id] += sensitivity_exceed[:, index]
                else:
                    primary_maxima_by_sv[sv_id] = primary_maxima[:, index].copy()
                    sensitivity_maxima_by_sv[sv_id] = sensitivity_maxima[
                        :, index
                    ].copy()
                    primary_exceed_by_sv[sv_id] = primary_exceed[:, index].copy()
                    sensitivity_exceed_by_sv[sv_id] = sensitivity_exceed[
                        :, index
                    ].copy()
            for index, key in enumerate(null["target_keys"].astype(str)):
                if key in target_null:
                    raise ValueError(
                        f"Target {key} appeared in more than one methylation chunk"
                    )
                target_null[key] = null["target_rademacher_statistics"][
                    :, index
                ].copy()

    if thresholds is None:
        raise ValueError("No multiplier-bootstrap outputs were found")
    p_t_all = np.concatenate(p_t_parts)
    p_asymptotic_all = np.concatenate(p_asymptotic_parts)
    score_all = np.concatenate(score_parts)
    primary_all = np.concatenate(primary_parts)
    q_sensitivity_t = bh_adjust(p_t_all)
    q_sensitivity_asymptotic = bh_adjust(p_asymptotic_all)
    q_primary_t = np.full(len(p_t_all), np.nan)
    q_primary_asymptotic = np.full(len(p_t_all), np.nan)
    q_primary_t[primary_all] = bh_adjust(p_t_all[primary_all])
    q_primary_asymptotic[primary_all] = bh_adjust(
        p_asymptotic_all[primary_all]
    )

    associations = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    by_job = {job_id: (start, stop) for job_id, start, stop in pvalue_files}
    if not associations.empty:
        global_indices = np.asarray(
            [
                by_job[str(job_id)][0] + int(local_index)
                for job_id, local_index in zip(
                    associations["job_id"],
                    associations["compact_test_index"],
                    strict=True,
                )
            ],
            dtype=np.int64,
        )
        associations["rint_hc3_bh_q_primary_t"] = q_primary_t[global_indices]
        associations[
            "rint_hc3_bh_q_primary_asymptotic"
        ] = q_primary_asymptotic[global_indices]
        associations[
            "rint_hc3_bh_q_sensitivity_t"
        ] = q_sensitivity_t[global_indices]
        associations[
            "rint_hc3_bh_q_sensitivity_asymptotic"
        ] = q_sensitivity_asymptotic[global_indices]

    sv_ids = sorted(primary_maxima_by_sv)
    primary_global = np.max(
        np.column_stack([primary_maxima_by_sv[sv_id] for sv_id in sv_ids]),
        axis=1,
    )
    sensitivity_global = np.max(
        np.column_stack([sensitivity_maxima_by_sv[sv_id] for sv_id in sv_ids]),
        axis=1,
    )
    if not associations.empty:
        associations["per_sv_primary_mammen_maxT_p"] = [
            empirical_p(primary_maxima_by_sv[sv_id], abs(score))
            if bool(primary)
            else np.nan
            for sv_id, score, primary in zip(
                associations["sv_id"],
                associations["rint_hc3_null_score"],
                associations["primary_family"],
                strict=True,
            )
        ]
        associations["global_primary_mammen_maxT_p"] = [
            empirical_p(primary_global, abs(score)) if bool(primary) else np.nan
            for score, primary in zip(
                associations["rint_hc3_null_score"],
                associations["primary_family"],
                strict=True,
            )
        ]
        associations["per_sv_sensitivity_rademacher_maxT_p"] = [
            empirical_p(sensitivity_maxima_by_sv[sv_id], abs(score))
            for sv_id, score in zip(
                associations["sv_id"],
                associations["rint_hc3_null_score"],
                strict=True,
            )
        ]
        associations["global_sensitivity_rademacher_maxT_p"] = [
            empirical_p(sensitivity_global, abs(score))
            for score in associations["rint_hc3_null_score"]
        ]
        associations["target_rademacher_pointwise_p"] = np.nan
        associations["target_rademacher_family_maxT_p"] = np.nan
        target_global = (
            np.max(np.column_stack(list(target_null.values())), axis=1)
            if target_null
            else np.asarray([])
        )
        target_rows = associations[
            associations["retained_reason"].astype(str).str.contains("target")
        ]
        for index, row in target_rows.iterrows():
            key = f"{row['sv_id']}|{row['feature_id']}"
            if key in target_null:
                associations.loc[
                    index, "target_rademacher_pointwise_p"
                ] = empirical_p(target_null[key], abs(row["rint_hc3_null_score"]))
                associations.loc[
                    index, "target_rademacher_family_maxT_p"
                ] = empirical_p(target_global, abs(row["rint_hc3_null_score"]))

    primary_null_counts = np.sum(
        [primary_exceed_by_sv[sv_id] for sv_id in sv_ids], axis=0
    )
    sensitivity_null_counts = np.sum(
        [sensitivity_exceed_by_sv[sv_id] for sv_id in sv_ids], axis=0
    )
    fdr_rows: list[dict[str, object]] = []
    primary_family_label = (
        f"primary_mac_ge_{args.primary_min_mac:g}_carriers_ge_{args.primary_min_carriers}"
    )
    sensitivity_family_label = (
        f"sensitivity_mac_ge_{args.min_mac:g}_carriers_ge_{args.min_carriers}_inclusive"
    )
    for family, observed_mask, null_counts, operator in (
        (primary_family_label, primary_all, primary_null_counts, "mammen"),
        (
            sensitivity_family_label,
            np.ones(len(primary_all), dtype=bool),
            sensitivity_null_counts,
            "rademacher",
        ),
    ):
        observed_scores = np.abs(score_all[observed_mask])
        for threshold_index, threshold in enumerate(thresholds):
            observed_count = int(np.count_nonzero(observed_scores >= threshold))
            mean_null = float(np.mean(null_counts[:, threshold_index]))
            fdr_rows.append(
                {
                    "family": family,
                    "operator": operator,
                    "abs_score_threshold": threshold,
                    "observed_discoveries": observed_count,
                    "mean_null_exceedances": mean_null,
                    "empirical_fdr": (
                        min(1.0, mean_null / observed_count)
                        if observed_count
                        else np.nan
                    ),
                }
            )
    empirical_fdr = pd.DataFrame(fdr_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    associations.to_csv(
        args.output_dir / "association_results.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    np.savez_compressed(
        args.output_dir / "full_family_qvalues.npz",
        job_ids=np.asarray([item[0] for item in pvalue_files]),
        job_offsets=np.asarray(job_offsets, dtype=np.int64),
        primary_family=primary_all,
        rint_hc3_bh_q_primary_t=q_primary_t,
        rint_hc3_bh_q_primary_asymptotic=q_primary_asymptotic,
        rint_hc3_bh_q_sensitivity_t=q_sensitivity_t,
        rint_hc3_bh_q_sensitivity_asymptotic=q_sensitivity_asymptotic,
        design_fingerprint=np.asarray(expected["design_fingerprint"]),
    )
    null_rows = []
    for sv_id in sv_ids:
        for replicate_index, seed in enumerate(expected["replicate_seeds"]):
            null_rows.append(
                {
                    "replicate_index": replicate_index,
                    "replicate_seed": seed,
                    "sv_id": sv_id,
                    "primary_mammen_max_abs_score": primary_maxima_by_sv[sv_id][
                        replicate_index
                    ],
                    "sensitivity_rademacher_max_abs_score": sensitivity_maxima_by_sv[
                        sv_id
                    ][replicate_index],
                }
            )
    pd.DataFrame(null_rows).to_csv(
        args.output_dir / "null_maxima_by_sv.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    empirical_fdr.to_csv(
        args.output_dir / "empirical_fdr_curve.tsv", sep="\t", index=False
    )
    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(
        args.output_dir / "job_summary.tsv", sep="\t", index=False
    )

    lead = pd.DataFrame()
    if not associations.empty:
        lead = associations.sort_values("rint_hc3_p_t").drop_duplicates("sv_id")
    lead.to_csv(args.output_dir / "lead_sv_results.tsv", sep="\t", index=False)
    targets = (
        associations[
            associations["retained_reason"].astype(str).str.contains("target")
        ].copy()
        if not associations.empty
        else pd.DataFrame()
    )
    targets.to_csv(args.output_dir / "target_validation.tsv", sep="\t", index=False)

    calibration = {
        "method_version": METHOD_VERSION,
        "analysis": (
            "RINT full-model HC3 t with mask-specific complete-case covariate "
            "adjustment; global-null HC3 score multiplier bootstrap"
        ),
        "maxT_scope": (
            "Asymptotic global-null maxT calibration; not an exact permutation "
            "test and not guaranteed strong FWER control under alternatives"
        ),
        "primary_family": (
            "actual methylation completeness >= configured fraction and "
            f">={args.primary_min_carriers} aligned analyzed dosage carriers and "
            f"complete-case dosage MAC >= {args.primary_min_mac:g}"
        ),
        "sensitivity_family": (
            "inclusive actual methylation completeness >= configured fraction, "
            f">={args.min_carriers} aligned analyzed dosage carriers, and "
            f"complete-case dosage MAC >= {args.min_mac:g}"
        ),
        "primary_multiplier": "Mammen",
        "sensitivity_multiplier": "Rademacher",
        "n_replicates": len(expected["replicate_seeds"]),
        "n_tests_primary": int(primary_all.sum()),
        "n_tests_sensitivity_inclusive": int(len(primary_all)),
        "pvalue_retention_threshold": args.pvalue_threshold,
        "min_carriers": args.min_carriers,
        "primary_min_carriers": args.primary_min_carriers,
        "min_mac": args.min_mac,
        "primary_min_mac": args.primary_min_mac,
        "min_feature_fraction": args.min_feature_fraction,
        "n_primary_bh_q_t_le_0_05": int(
            np.count_nonzero(q_primary_t[primary_all] <= 0.05)
        ),
        "n_primary_bh_q_asymptotic_le_0_05": int(
            np.count_nonzero(q_primary_asymptotic[primary_all] <= 0.05)
        ),
        "n_sensitivity_bh_q_t_le_0_05": int(
            np.count_nonzero(q_sensitivity_t <= 0.05)
        ),
        "n_sensitivity_bh_q_asymptotic_le_0_05": int(
            np.count_nonzero(q_sensitivity_asymptotic <= 0.05)
        ),
        "minimum_empirical_p": 1.0 / (len(expected["replicate_seeds"]) + 1.0),
        "sample_order_fingerprint": expected["sample_order_fingerprint"],
        "design_fingerprint": expected["design_fingerprint"],
        "mammen_fingerprint": expected["mammen_fingerprint"],
        "rademacher_fingerprint": expected["rademacher_fingerprint"],
    }
    (args.output_dir / "calibration.json").write_text(
        json.dumps(calibration, indent=2) + "\n"
    )
    pd.DataFrame([calibration]).to_csv(
        args.output_dir / "calibration.tsv", sep="\t", index=False
    )
    print(json.dumps(calibration, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
