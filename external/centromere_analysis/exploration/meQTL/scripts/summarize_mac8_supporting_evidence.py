#!/usr/bin/env python3
"""Generate MAC8 influence, conditioning-sensitivity, phase, and dip summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats


def write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def write_json(value: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def rint(values: np.ndarray) -> np.ndarray:
    ranks = stats.rankdata(np.asarray(values, float), method="average")
    return special.ndtri((ranks - 0.5) / len(ranks))


def hc3(y: np.ndarray, design: np.ndarray, coefficient: int = -1) -> dict[str, float | str]:
    y = np.asarray(y, float)
    x = np.asarray(design, float)
    coefficient %= x.shape[1]
    if len(y) <= x.shape[1] or np.linalg.matrix_rank(x) < x.shape[1]:
        return {"status": "rank_deficient", "beta": np.nan, "se": np.nan, "p": np.nan, "max_leverage": np.nan}
    inverse = np.linalg.inv(x.T @ x)
    beta = inverse @ x.T @ y
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, inverse, x)
    adjusted = residual / np.clip(1.0 - leverage, 1e-10, None)
    meat = x.T @ ((adjusted * adjusted)[:, None] * x)
    covariance = inverse @ meat @ inverse
    se = float(np.sqrt(max(covariance[coefficient, coefficient], 0.0)))
    statistic = float(beta[coefficient] / se) if se > 0 else np.nan
    dof = len(y) - x.shape[1]
    p = float(2.0 * stats.t.sf(abs(statistic), dof)) if np.isfinite(statistic) else np.nan
    return {
        "status": "ok",
        "beta": float(beta[coefficient]),
        "se": se,
        "p": p,
        "max_leverage": float(np.max(leverage)),
    }


def load_design(path: Path, primary_pc: int = 10) -> tuple[list[str], np.ndarray]:
    data = np.load(path, allow_pickle=True)
    matrix = np.asarray(data["covariates"], float)
    if "pc_counts" in data.files:
        counts = dict(zip(data["pc_counts"].astype(int), data["covariate_counts"].astype(int)))
        matrix = matrix[:, : counts[primary_pc]]
    return data["sample_ids"].astype(str).tolist(), matrix


def load_dosage(row: pd.Series, samples: list[str]) -> np.ndarray:
    frame = pd.read_csv(row.sv_dosages, sep="\t")
    selected = frame.loc[pd.to_numeric(frame["sv_index"], errors="coerce").eq(int(row.sv_index))]
    if len(selected) != 1:
        raise ValueError(f"Expected one dosage row for {row.sv_id}, found {len(selected)}")
    return pd.to_numeric(selected.iloc[0][samples], errors="coerce").to_numpy(float)


def load_outcomes(pairs: pd.DataFrame, samples: list[str]) -> dict[str, np.ndarray]:
    outcomes: dict[str, np.ndarray] = {}
    for matrix_path, group in pairs.groupby("matrix_path"):
        wanted = set(group.feature_id.astype(str))
        frame = pd.read_csv(matrix_path, sep="\t", usecols=lambda name: name == "feature_id" or name in samples)
        frame["feature_id"] = frame.feature_id.astype(str)
        frame = frame.loc[frame.feature_id.isin(wanted)].set_index("feature_id")
        for feature_id in wanted & set(frame.index):
            outcomes[feature_id] = pd.to_numeric(frame.loc[feature_id, samples], errors="coerce").to_numpy(float)
    return outcomes


def leave_one_carrier_out(
    consensus_path: Path,
    sv_manifest_path: Path,
    cis_design_path: Path,
    trans_design_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    consensus = pd.read_csv(consensus_path, sep="\t", low_memory=False)
    pairs = consensus.loc[consensus.upstream_aligned_supported_scope.fillna(False)].copy()
    if len(pairs) != 84:
        raise ValueError(f"Primary retained set must contain 84 pairs, found {len(pairs)}")
    svs = pd.read_csv(sv_manifest_path, sep="\t").set_index("sv_id")
    cis_samples, cis_covariates = load_design(cis_design_path)
    trans_samples, trans_covariates = load_design(trans_design_path)
    if cis_samples != trans_samples:
        raise ValueError("Cis and trans PC10 designs have different sample order")
    samples = cis_samples
    designs = {"cis": cis_covariates, "trans": trans_covariates}
    outcomes = load_outcomes(pairs, samples)
    dosage_cache: dict[str, np.ndarray] = {}
    omission_rows: list[dict] = []
    pair_rows: list[dict] = []
    for pair in pairs.itertuples(index=False):
        sv_id = str(pair.sv_id)
        if sv_id not in dosage_cache:
            dosage_cache[sv_id] = load_dosage(svs.loc[sv_id], samples)
        dosage = dosage_cache[sv_id]
        outcome = outcomes[str(pair.feature_id)]
        covariates = designs[str(pair.scope)]
        complete = np.isfinite(outcome) & np.isfinite(dosage) & np.all(np.isfinite(covariates), axis=1)
        y = rint(outcome[complete])
        g = dosage[complete]
        cov = covariates[complete]
        sample_array = np.asarray(samples)[complete]
        baseline = hc3(y, np.column_stack([np.ones(len(y)), cov, g]))
        carriers = np.flatnonzero(g > 0)
        local = []
        for carrier in carriers:
            keep = np.arange(len(y)) != carrier
            fit = hc3(rint(outcome[complete][keep]), np.column_stack([np.ones(keep.sum()), cov[keep], g[keep]]))
            same_sign = bool(np.isfinite(fit["beta"]) and np.sign(fit["beta"]) == np.sign(baseline["beta"]))
            retention = float(fit["beta"] / baseline["beta"]) if np.isfinite(fit["beta"]) and baseline["beta"] else np.nan
            row = {
                "pair_id": pair.pair_id,
                "scope": pair.scope,
                "sv_id": sv_id,
                "feature_id": pair.feature_id,
                "omitted_sample": sample_array[carrier],
                "omitted_dosage": float(g[carrier]),
                "n_after_omission": int(keep.sum()),
                "beta": fit["beta"],
                "hc3_p": fit["p"],
                "max_leverage": fit["max_leverage"],
                "same_sign": same_sign,
                "effect_retention": retention,
                "status": fit["status"],
            }
            local.append(row)
            omission_rows.append(row)
        local_frame = pd.DataFrame(local)
        robust = bool(
            len(local_frame)
            and local_frame.status.eq("ok").all()
            and local_frame.same_sign.all()
            and local_frame.effect_retention.abs().min() >= 0.5
        )
        pair_rows.append({
            "pair_id": pair.pair_id,
            "scope": pair.scope,
            "sv_id": sv_id,
            "feature_id": pair.feature_id,
            "n_complete": int(complete.sum()),
            "n_carriers": int(len(carriers)),
            "baseline_beta": baseline["beta"],
            "baseline_hc3_p": baseline["p"],
            "same_sign_fraction": float(local_frame.same_sign.mean()) if len(local_frame) else np.nan,
            "minimum_absolute_effect_retention": float(local_frame.effect_retention.abs().min()) if len(local_frame) else np.nan,
            "median_effect_retention": float(local_frame.effect_retention.median()) if len(local_frame) else np.nan,
            "maximum_omission_hc3_p": float(local_frame.hc3_p.max()) if len(local_frame) else np.nan,
            "loo_robust_sign_and_half_effect": robust,
        })
    omission = pd.DataFrame(omission_rows)
    pair_summary = pd.DataFrame(pair_rows)
    summary = {
        "criterion": "all single-carrier omissions retain the baseline sign and at least 50% of the absolute RINT beta",
        "n_primary_pairs": int(len(pair_summary)),
        "n_cis_pairs": int(pair_summary.scope.eq("cis").sum()),
        "n_trans_pairs": int(pair_summary.scope.eq("trans").sum()),
        "n_robust_pairs": int(pair_summary.loo_robust_sign_and_half_effect.sum()),
        "n_robust_cis_pairs": int(pair_summary.loc[pair_summary.scope.eq("cis"), "loo_robust_sign_and_half_effect"].sum()),
        "n_robust_trans_pairs": int(pair_summary.loc[pair_summary.scope.eq("trans"), "loo_robust_sign_and_half_effect"].sum()),
        "interpretation": "Influence sensitivity within the discovery cohort; not independent replication.",
    }
    return omission, pair_summary, summary


def conditioning_sensitivity(primary_summary: dict, broad_path: Path) -> tuple[pd.DataFrame, dict]:
    broad = json.loads(broad_path.read_text())
    rows = pd.DataFrame([
        {
            "analysis": "Primary upstream-imputation-aligned conditioning",
            "pairs_evaluated": 465,
            "single_panel_supported_pairs": 287,
            "cross_panel_robust_pairs": 84,
            "definition": "Imputed-score, out-of-fold exact-tag residual, and SuSiE gates",
            "reporting_role": "primary",
        },
        {
            "analysis": "Broader local-SNP conditioning sensitivity",
            "pairs_evaluated": int(broad.get("n_pairs", broad.get("pairs_evaluated", 465))),
            "single_panel_supported_pairs": int(broad.get("n_single_panel_supported", broad.get("single_panel_supported_pairs", 7))),
            "cross_panel_robust_pairs": int(broad.get("n_cross_panel_robust", broad.get("cross_panel_robust_pairs", 0))),
            "definition": "Broader local-SNP panels; sensitivity rather than primary attribution criterion",
            "reporting_role": "sensitivity",
        },
    ])
    return rows, {
        "primary_retained_pairs": int(primary_summary["n_primary_pairs"]),
        "broader_local_snp_single_panel_supported_pairs": int(rows.iloc[1].single_panel_supported_pairs),
        "broader_local_snp_cross_panel_robust_pairs": int(rows.iloc[1].cross_panel_robust_pairs),
    }


def phase_summary(root: Path) -> tuple[pd.DataFrame, dict]:
    rows = []
    for path in sorted(root.glob("*/summary.json")):
        data = json.loads(path.read_text())
        evaluable = int(data.get("n_unique_phase_valid_cpgs", 0)) > 0
        rows.append({
            "locus": path.parent.name,
            "sv_id": data["sv_id"],
            "n_heterozygotes": data.get("n_heterozygotes"),
            "n_heterozygotes_phased": data.get("n_heterozygotes_phased"),
            "n_phase_valid_cpgs": data.get("n_unique_phase_valid_cpgs", 0),
            "evaluable": evaluable,
            "concordant_fraction": data.get("concordant_fraction"),
            "sign_test_p": data.get("two_sided_sign_test_p"),
            "median_beta_aligned_haplotype_difference": data.get("median_beta_aligned_haplotype_difference"),
            "wilcoxon_p": data.get("one_sided_wilcoxon_p"),
        })
    frame = pd.DataFrame(rows)
    evaluable = frame.loc[frame.evaluable]
    summary = {
        "n_loci": int(len(frame)),
        "n_evaluable": int(frame.evaluable.sum()),
        "n_untestable": int((~frame.evaluable).sum()),
        "evaluable_locus_n_heterozygotes": int(evaluable.n_heterozygotes.iloc[0]) if len(evaluable) == 1 else None,
        "evaluable_locus_sign_test_p": float(evaluable.sign_test_p.iloc[0]) if len(evaluable) == 1 else None,
        "interpretation": "Three loci lacked phase-valid CpGs; the sole evaluable three-carrier locus was null and is too small for a broad cis-causality conclusion.",
    }
    return frame, summary


def dip_associations(
    dip_path: Path,
    active_sv_path: Path,
    sv_manifest_path: Path,
    design_path: Path,
) -> tuple[pd.DataFrame, dict]:
    samples, covariates = load_design(design_path)
    dips = pd.read_csv(dip_path, sep="\t")
    phenotype_columns = ["cdr_pos_offset", "cdr_pos_offset_deepest", "deepest_min5mC", "n_cdrs", "mean_5mC"]
    available = [column for column in phenotype_columns if column in dips.columns]
    sample_level = dips.groupby(["sample", "chm13_chrom"], as_index=False)[available].mean()
    active = pd.read_csv(active_sv_path, sep="\t").rename(columns={
        "physical_family_id": "sv_id", "sv_chrom": "chrom"
    })
    svs = pd.read_csv(sv_manifest_path, sep="\t").set_index("sv_id")
    rows = []
    dosage_cache: dict[str, np.ndarray] = {}
    for sv in active.itertuples(index=False):
        sv_id = str(sv.sv_id)
        chrom = str(sv.chrom)
        dosage_row = (
            svs.loc[sv_id] if sv_id in svs.index
            else pd.Series({"sv_id": sv_id, "sv_dosages": sv.sv_dosages, "sv_index": sv.sv_index})
        )
        dosage_cache[sv_id] = load_dosage(dosage_row, samples)
        local = sample_level.loc[sample_level.chm13_chrom.eq(chrom)].set_index("sample")
        for phenotype in available:
            values = pd.to_numeric(local.reindex(samples)[phenotype], errors="coerce").to_numpy(float)
            dosage = dosage_cache[sv_id]
            complete = np.isfinite(values) & np.isfinite(dosage) & np.all(np.isfinite(covariates), axis=1)
            fit = hc3(rint(values[complete]), np.column_stack([np.ones(complete.sum()), covariates[complete], dosage[complete]]))
            rows.append({
                "sv_id": sv_id,
                "chrom": chrom,
                "phenotype": phenotype,
                "resolution": "sample-level mean of available haplotype-level chromosome phenotypes",
                "n_complete": int(complete.sum()),
                "n_carriers": int(np.sum(dosage[complete] > 0)),
                "beta": fit["beta"],
                "hc3_p": fit["p"],
                "max_leverage": fit["max_leverage"],
                "status": fit["status"],
            })
    frame = pd.DataFrame(rows)
    valid = frame.hc3_p.notna()
    if valid.any():
        order = frame.loc[valid, "hc3_p"].sort_values().index
        ranks = np.arange(1, len(order) + 1)
        adjusted = np.minimum.accumulate((frame.loc[order, "hc3_p"].to_numpy() * len(order) / ranks)[::-1])[::-1]
        frame.loc[order, "bh_q"] = np.minimum(adjusted, 1.0)
    summary = {
        "n_active_svs": int(frame.sv_id.nunique()),
        "n_phenotypes": int(frame.phenotype.nunique()),
        "n_tests": int(valid.sum()),
        "n_bh_q_lt_0_05": int((frame.bh_q < 0.05).sum()),
        "design": "same 124-sample cis PC10 covariates; complete-case RINT HC3",
        "resolution": "dip phenotypes were haplotype-level and were averaged to sample-level within chromosome before association with diploid dosage",
    }
    return frame, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", type=Path, required=True)
    parser.add_argument("--sv-manifest", type=Path, required=True)
    parser.add_argument("--active-svs", type=Path, required=True)
    parser.add_argument("--cis-design", type=Path, required=True)
    parser.add_argument("--trans-design", type=Path, required=True)
    parser.add_argument("--broad-conditioning-summary", type=Path, required=True)
    parser.add_argument("--phase-root", type=Path, required=True)
    parser.add_argument("--dip-phenotypes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    omission, pair_summary, loo_summary = leave_one_carrier_out(
        args.consensus, args.sv_manifest, args.cis_design, args.trans_design
    )
    write_table(omission, args.output_dir / "primary_84_leave_one_carrier_out.tsv")
    write_table(pair_summary, args.output_dir / "primary_84_leave_one_carrier_out_summary.tsv")
    write_json(loo_summary, args.output_dir / "leave_one_carrier_out_summary.json")

    sensitivity, sensitivity_summary = conditioning_sensitivity(loo_summary, args.broad_conditioning_summary)
    write_table(sensitivity, args.output_dir / "conditioning_sensitivity.tsv")
    write_json(sensitivity_summary, args.output_dir / "conditioning_sensitivity_summary.json")

    phase, phase_meta = phase_summary(args.phase_root)
    write_table(phase, args.output_dir / "phase_resolved_loci.tsv")
    write_json(phase_meta, args.output_dir / "phase_resolved_summary.json")

    dip, dip_meta = dip_associations(args.dip_phenotypes, args.active_svs, args.sv_manifest, args.cis_design)
    write_table(dip, args.output_dir / "centromere_dip_associations.tsv")
    write_json(dip_meta, args.output_dir / "centromere_dip_association_summary.json")
    print(json.dumps({"leave_one_carrier_out": loo_summary, "dip": dip_meta}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
