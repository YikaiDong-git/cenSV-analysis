#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import yaml
from scipy import special, stats


BED_COVERAGE_COL = 9
BED_METHYLATION_COL = 10
SATELLITE_ORDER = ["live_hor", "divergent_hor", "monomeric_alpha", "other_satellite", "centromere_transition", "non_satellite"]
SCOPE_COLORS = {"local_satellite": "#007A78", "distal_cis": "#6B7280", "trans": "#C2410C"}


def read_table(path: str | Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, **kwargs)


def write_table(frame: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def parse_feature_id(value: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"(chr(?:[0-9]+|X|Y)):(\d+)-(\d+)", str(value))
    if match is None:
        raise ValueError(f"Unparseable CpG feature ID: {value}")
    return match.group(1), int(match.group(2)), int(match.group(3))


def parse_methylation_fraction(value: str | float) -> float:
    result = float(value)
    return result / 100.0 if result > 1.0 else result


def classify_censat_name(name: str) -> tuple[str, str]:
    stem, _, detail = str(name).partition("(")
    detail = detail.rstrip(")")
    if stem.startswith("hor_"):
        category = "live_hor"
    elif stem.startswith("dhor_"):
        category = "divergent_hor"
    elif stem.startswith("mon_"):
        category = "monomeric_alpha"
    elif stem.startswith(("hsat_", "bsat_", "gsat_", "censat_")):
        category = "other_satellite"
    elif stem.startswith("ct_"):
        category = "centromere_transition"
    else:
        category = "other_satellite"
    return category, detail or stem


def load_censat(path: str | Path) -> dict[str, list[tuple[int, int, str, str, str]]]:
    intervals: dict[str, list[tuple[int, int, str, str, str]]] = defaultdict(list)
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            category, array_name = classify_censat_name(fields[3])
            intervals[fields[0]].append((int(fields[1]), int(fields[2]), fields[3], category, array_name))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def annotate_censat_interval(
    chrom: str,
    start: int,
    end: int,
    intervals: dict[str, list[tuple[int, int, str, str, str]]],
) -> dict[str, object]:
    start, end = int(start), max(int(start) + 1, int(end))
    candidates = intervals.get(str(chrom), [])
    overlaps = []
    nearest = None
    nearest_distance = math.inf
    for left, right, name, category, array_name in candidates:
        overlap = max(0, min(end, right) - max(start, left))
        if overlap:
            overlaps.append((overlap, -(right - left), left, right, name, category, array_name))
        distance = max(left - end, start - right, 0)
        if distance < nearest_distance:
            nearest_distance = distance
            nearest = (left, right, name, category, array_name)
    if overlaps:
        _, _, left, right, name, category, array_name = max(overlaps)
        boundary_distance = min(abs(start - left), abs(start - right), abs(end - left), abs(end - right))
        overlap_bp = max(0, min(end, right) - max(start, left))
        relation = "overlap"
    elif nearest is not None:
        left, right, name, _, array_name = nearest
        category = "non_satellite"
        boundary_distance = int(nearest_distance)
        overlap_bp = 0
        relation = "nearest"
    else:
        left = right = -1
        name = array_name = ""
        category = "non_satellite"
        boundary_distance = math.nan
        overlap_bp = 0
        relation = "unannotated"
    return {
        "satellite_class": category,
        "censat_relation": relation,
        "censat_feature": name,
        "array_name": array_name,
        "array_start": left,
        "array_end": right,
        "overlap_bp": overlap_bp,
        "distance_to_array_boundary_bp": boundary_distance,
    }


def monomer_register_distance(length: float, period: int = 171) -> float:
    if not np.isfinite(length) or length <= 0:
        return math.nan
    remainder = float(length) % period
    return min(remainder, period - remainder)


def periodicity_permutation(
    svs: pd.DataFrame,
    period: int = 171,
    tolerance: int = 10,
    permutations: int = 10_000,
    seed: int = 20260829,
) -> dict[str, object]:
    table = svs.copy()
    table["abs_svlen"] = pd.to_numeric(table["svlen"], errors="coerce").abs()
    table["register_distance_bp"] = table["abs_svlen"].map(lambda x: monomer_register_distance(x, period))
    table["periodic"] = table["register_distance_bp"].le(tolerance)
    live = table[(table["satellite_class"] == "live_hor") & table["abs_svlen"].ge(period)].copy()
    controls = table[(table["satellite_class"] != "live_hor") & table["abs_svlen"].ge(period)].copy()
    matched: list[np.ndarray] = []
    kept_live = []
    for row in live.itertuples(index=False):
        local = controls[controls["svtype"].astype(str).eq(str(row.svtype))]
        ratio = np.abs(np.log2(local["abs_svlen"].to_numpy(float) / float(row.abs_svlen)))
        indices = np.flatnonzero(np.isfinite(ratio) & (ratio <= 1.0))
        if len(indices):
            matched.append(local.iloc[indices].index.to_numpy())
            kept_live.append(row)
    if not kept_live:
        return {"status": "not_estimable", "n_live_hor": 0, "period_bp": period, "tolerance_bp": tolerance}
    observed_live = np.mean([bool(row.periodic) for row in kept_live])
    rng = np.random.default_rng(seed)
    null = np.empty(permutations)
    for index in range(permutations):
        selected = [rng.choice(candidates) for candidates in matched]
        null[index] = observed_live - float(table.loc[selected, "periodic"].mean())
    observed_controls = np.mean([table.loc[candidates, "periodic"].mean() for candidates in matched])
    observed = observed_live - observed_controls
    return {
        "status": "estimated",
        "period_bp": period,
        "tolerance_bp": tolerance,
        "n_live_hor": len(kept_live),
        "live_hor_periodic_fraction": float(observed_live),
        "matched_control_periodic_fraction": float(observed_controls),
        "fraction_difference": float(observed),
        "empirical_p_enrichment": float((1 + np.sum(null >= observed)) / (permutations + 1)),
        "permutations": permutations,
        "matching": "SV type and within one two-fold length band",
    }


def direct_rint(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    ranks = stats.rankdata(values, method="average")
    return special.ndtri((ranks - 0.5) / len(values))


def hc3_fit(y: np.ndarray, design: np.ndarray, coefficient: int = -1, weights: np.ndarray | None = None) -> dict[str, float | str]:
    y = np.asarray(y, float)
    x = np.asarray(design, float)
    if weights is not None:
        root = np.sqrt(np.asarray(weights, float))
        y = y * root
        x = x * root[:, None]
    rank = int(np.linalg.matrix_rank(x))
    df = len(y) - rank
    if rank < x.shape[1] or df <= 0:
        return {"status": "rank_deficient", "beta": math.nan, "se": math.nan, "p": math.nan, "max_leverage": math.nan}
    bread = np.linalg.inv(x.T @ x)
    coefficients = bread @ x.T @ y
    residual = y - x @ coefficients
    leverage = np.einsum("ij,jk,ik->i", x, bread, x)
    adjusted = residual / np.clip(1.0 - leverage, 1e-8, None)
    meat = (x * adjusted[:, None]).T @ (x * adjusted[:, None])
    covariance = bread @ meat @ bread
    index = coefficient if coefficient >= 0 else x.shape[1] + coefficient
    se = float(np.sqrt(max(0.0, covariance[index, index])))
    t_value = float(coefficients[index] / se) if se else math.nan
    return {
        "status": "ok" if np.isfinite(t_value) else "nonfinite",
        "beta": float(coefficients[index]),
        "se": se,
        "p": float(2 * stats.t.sf(abs(t_value), df)) if np.isfinite(t_value) else math.nan,
        "max_leverage": float(leverage.max()),
    }


def binomial_hc3_fit(successes: np.ndarray, trials: np.ndarray, design: np.ndarray, coefficient: int = -1) -> dict[str, float | str]:
    successes, trials, x = np.asarray(successes, float), np.asarray(trials, float), np.asarray(design, float)
    beta = np.zeros(x.shape[1])
    for _ in range(80):
        eta = np.clip(x @ beta, -20, 20)
        mu = special.expit(eta)
        variance = np.clip(trials * mu * (1 - mu), 1e-8, None)
        working = eta + (successes - trials * mu) / variance
        bread = np.linalg.pinv((x * np.sqrt(variance)[:, None]).T @ (x * np.sqrt(variance)[:, None]), rcond=1e-12)
        updated = bread @ (x.T @ (variance * working))
        if np.max(np.abs(updated - beta)) < 1e-9:
            beta = updated
            break
        beta = updated
    mu = special.expit(np.clip(x @ beta, -20, 20))
    variance = np.clip(trials * mu * (1 - mu), 1e-8, None)
    weighted_x = x * np.sqrt(variance)[:, None]
    rank = int(np.linalg.matrix_rank(weighted_x))
    if rank < x.shape[1]:
        return {"status": "rank_deficient", "beta": math.nan, "se": math.nan, "p": math.nan, "max_leverage": math.nan}
    bread = np.linalg.inv(weighted_x.T @ weighted_x)
    leverage = np.einsum("ij,jk,ik->i", weighted_x, bread, weighted_x)
    score = x * ((successes - trials * mu) / np.clip(1 - leverage, 1e-8, None))[:, None]
    covariance = bread @ (score.T @ score) @ bread
    index = coefficient if coefficient >= 0 else x.shape[1] + coefficient
    se = float(np.sqrt(max(0.0, covariance[index, index])))
    z_value = float(beta[index] / se) if se else math.nan
    return {"status": "ok", "beta": float(beta[index]), "se": se, "p": float(2 * stats.norm.sf(abs(z_value))), "max_leverage": float(leverage.max())}


def model_pair(
    methylation: np.ndarray,
    depth: np.ndarray,
    dosage: np.ndarray,
    covariates: np.ndarray,
    coverage_min: float,
) -> dict[str, object]:
    finite_base = np.isfinite(dosage) & np.all(np.isfinite(covariates), axis=1)
    callable_mask = finite_base & np.isfinite(methylation) & np.isfinite(depth) & (depth >= coverage_min)
    depth_mask = finite_base & np.isfinite(depth) & (depth > 0)
    result: dict[str, object] = {"n_callable": int(callable_mask.sum()), "n_depth_observed": int(depth_mask.sum())}
    base_design = np.column_stack([np.ones(callable_mask.sum()), covariates[callable_mask], dosage[callable_mask]])
    y = methylation[callable_mask]
    if len(y) <= base_design.shape[1] or np.nanstd(dosage[callable_mask]) == 0:
        result["status"] = "not_estimable"
        return result
    y_rint = direct_rint(y)
    baseline = hc3_fit(y_rint, base_design)
    depth_design = np.column_stack([np.ones(depth_mask.sum()), covariates[depth_mask], dosage[depth_mask]])
    depth_fit = hc3_fit(np.log1p(depth[depth_mask]), depth_design)
    passage = (np.nan_to_num(depth, nan=0.0) >= coverage_min).astype(float)
    passage_design = np.column_stack([np.ones(finite_base.sum()), covariates[finite_base], dosage[finite_base]])
    passage_fit = hc3_fit(passage[finite_base], passage_design)
    log_depth = np.log1p(depth[callable_mask])
    adjusted_design = np.column_stack([np.ones(len(y)), covariates[callable_mask], log_depth, dosage[callable_mask]])
    adjusted = hc3_fit(y_rint, adjusted_design)
    weights = np.clip(depth[callable_mask], 1, np.nanquantile(depth[callable_mask], 0.95))
    weighted = hc3_fit(y_rint, base_design, weights=weights)
    successes = np.rint(np.clip(y, 0, 1) * depth[callable_mask])
    count_aware = binomial_hc3_fit(successes, depth[callable_mask], base_design)
    result.update({
        "status": "ok",
        "baseline_beta": baseline["beta"], "baseline_p": baseline["p"],
        "depth_beta": depth_fit["beta"], "depth_p": depth_fit["p"],
        "threshold_pass_beta": passage_fit["beta"], "threshold_pass_p": passage_fit["p"],
        "depth_adjusted_beta": adjusted["beta"], "depth_adjusted_p": adjusted["p"],
        "coverage_weighted_beta": weighted["beta"], "coverage_weighted_p": weighted["p"],
        "count_aware_log_odds_beta": count_aware["beta"], "count_aware_p": count_aware["p"],
        "depth_adjusted_effect_retention": adjusted["beta"] / baseline["beta"] if baseline["beta"] else math.nan,
        "coverage_weighted_effect_retention": weighted["beta"] / baseline["beta"] if baseline["beta"] else math.nan,
    })
    return result


def deterministic_control_set(path: str | Path) -> pd.DataFrame:
    controls = read_table(path)
    controls = controls.sort_values(["matched_feature_id", "matching_level", "match_distance", "feature_id"])
    controls = controls.drop_duplicates("matched_feature_id", keep="first").copy()
    parsed = controls["feature_id"].map(parse_feature_id)
    controls["chrom"] = [x[0] for x in parsed]
    controls["start"] = [x[1] for x in parsed]
    controls["end"] = [x[2] for x in parsed]
    return controls


def prepare(args: argparse.Namespace) -> None:
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    pairs = read_table(args.pairs)
    cis = read_table(args.cis_associations)
    cis["rint_hc3_p_t"] = pd.to_numeric(cis["rint_hc3_p_t"], errors="coerce")
    cis["rint_hc3_bh_q_primary_t"] = pd.to_numeric(cis["rint_hc3_bh_q_primary_t"], errors="coerce")
    tail = cis[cis["rint_hc3_p_t"] <= args.tail_p].copy()
    jobs = read_table(args.source_jobs)[["job_id", "sv_dosages"]].drop_duplicates()
    tail = tail.merge(jobs, on="job_id", how="left", validate="many_to_one")
    audit_pairs = pd.DataFrame({
        "pair_id": [f"tail{index:06d}" for index in range(len(tail))],
        "sv_id": tail["sv_id"].astype(str), "feature_id": tail["feature_id"].astype(str),
        "scope": "cis", "source_chrom": tail["sv_chrom"].astype(str),
        "source_pos": tail["sv_pos"].astype(int), "source_end": tail["sv_end"].astype(int),
        "target_chrom": tail["methyl_chrom"].astype(str), "target_start": tail["methyl_start"].astype(int),
        "target_end": tail["methyl_end"].astype(int), "sv_index": tail["sv_index"].astype(int),
        "sv_dosages": tail["sv_dosages"].astype(str), "expected_rint_p": tail["rint_hc3_p_t"],
        "expected_rint_beta": tail["rint_beta"], "expected_raw_beta": tail["raw_beta"],
        "discovery_q": tail["rint_hc3_bh_q_primary_t"],
        "is_discovery": tail["rint_hc3_bh_q_primary_t"] < args.discovery_q,
    })
    write_table(audit_pairs, output / "audit_pairs.tsv.gz")

    features = audit_pairs[["feature_id", "target_chrom", "target_start", "target_end"]].drop_duplicates().rename(
        columns={"target_chrom": "chrom", "target_start": "start", "target_end": "end"}
    )
    features["feature_role"] = "cis_tail"
    controls = deterministic_control_set(args.matched_controls)
    control_features = controls[["feature_id", "chrom", "start", "end", "matched_feature_id"]].copy()
    control_features["feature_role"] = "matched_control"
    features = pd.concat([features, control_features], ignore_index=True, sort=False).drop_duplicates("feature_id")
    write_table(features, output / "depth_features.tsv.gz")

    censat = load_censat(args.censat)
    cpg_rows = []
    for row in features.itertuples(index=False):
        cpg_rows.append({**row._asdict(), **annotate_censat_interval(row.chrom, row.start, row.end, censat)})
    write_table(pd.DataFrame(cpg_rows), output / "cpg_censat_annotations.tsv.gz")

    tested_svs = read_table(args.tested_svs).copy()
    sv_rows = []
    for row in tested_svs.itertuples(index=False):
        annotation = annotate_censat_interval(row.sv_chrom, row.sv_pos, row.sv_end, censat)
        sv_rows.append({**row._asdict(), **annotation, "register_distance_171bp": monomer_register_distance(abs(float(row.svlen)))})
    tested_svs = pd.DataFrame(sv_rows)
    write_table(tested_svs, output / "tested_sv_censat_annotations.tsv")
    periodicity = periodicity_permutation(tested_svs, args.period_bp, args.period_tolerance, args.permutations)
    summary = {
        "n_tail_pairs": int(len(audit_pairs)), "n_tail_cpgs": int(audit_pairs.feature_id.nunique()),
        "n_discovery_pairs": int(audit_pairs.is_discovery.sum()), "n_depth_features": int(len(features)),
        "n_matched_controls": int((features.feature_role == "matched_control").sum()),
        "n_tested_svs": int(len(tested_svs)), "periodicity": periodicity,
        "resolution": "SV dosage is diploid sample-level; hp1/hp2 bedMethyl depth is summed before association",
    }
    (output / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def bed_paths_for_sample(row: pd.Series) -> list[tuple[Path, Path | None]]:
    beds = [Path(value) for value in str(row["methylation_beds"]).split(";") if value]
    indexes = [Path(value) for value in str(row.get("methylation_indexes", "")).split(";") if value and value != "nan"]
    indexes.extend([None] * (len(beds) - len(indexes)))
    return list(zip(beds, indexes, strict=True))


def extract_site(handles: Iterable[pysam.TabixFile], chrom: str, start: int, end: int) -> tuple[float, float]:
    methylated = 0.0
    depth = 0.0
    for handle in handles:
        try:
            lines = handle.fetch(chrom, int(start), int(end))
        except (ValueError, OSError):
            continue
        for line in lines:
            fields = line.split("\t")
            if len(fields) <= BED_METHYLATION_COL or int(fields[1]) != int(start):
                continue
            row_depth = float(fields[BED_COVERAGE_COL])
            methylated += parse_methylation_fraction(fields[BED_METHYLATION_COL]) * row_depth
            depth += row_depth
    return (methylated / depth if depth > 0 else math.nan, depth if depth > 0 else math.nan)


def extract_depth(args: argparse.Namespace) -> None:
    features = read_table(args.features)
    features = features[features["chrom"].astype(str).eq(args.chrom)].copy()
    manifest = read_table(args.methylation_manifest).set_index("sample_id")
    design = np.load(args.design, allow_pickle=True)
    samples = design["sample_ids"].astype(str).tolist()
    missing = sorted(set(samples) - set(manifest.index.astype(str)))
    if missing:
        raise ValueError(f"{len(missing)} canonical samples are absent from the methylation manifest")
    rows = []
    for sample in samples:
        opened = []
        try:
            for bed, index in bed_paths_for_sample(manifest.loc[sample]):
                opened.append(pysam.TabixFile(str(bed), index=str(index) if index else None))
            for feature in features.itertuples(index=False):
                methylation, depth = extract_site(opened, args.chrom, int(feature.start), int(feature.end))
                rows.append({"feature_id": feature.feature_id, "sample_id": sample, "methylation": methylation, "depth": depth})
        finally:
            for handle in opened:
                handle.close()
    write_table(pd.DataFrame(rows), args.output)
    payload = {"chrom": args.chrom, "n_features": int(len(features)), "n_samples": len(samples), "n_rows": len(rows)}
    Path(args.summary).write_text(json.dumps(payload, indent=2) + "\n")


def load_dosage_rows(pairs: pd.DataFrame, samples: list[str]) -> dict[tuple[str, int], np.ndarray]:
    cache: dict[Path, pd.DataFrame] = {}
    result = {}
    for path, group in pairs.groupby("sv_dosages"):
        path = Path(path)
        frame = read_table(path)
        if "sv_index" in frame.columns:
            frame = frame.set_index("sv_index")
        for row in group[["sv_id", "sv_index"]].drop_duplicates().itertuples(index=False):
            result[(str(row.sv_id), int(row.sv_index))] = pd.to_numeric(frame.loc[int(row.sv_index), samples], errors="coerce").to_numpy(float)
    return result


def simulate_null_negative_fraction(
    methylation: np.ndarray,
    depth: np.ndarray,
    dosage: np.ndarray,
    covariates: np.ndarray,
    coverage_min: float,
    simulations: int,
    seed: int,
) -> np.ndarray:
    mask = np.isfinite(methylation) & np.isfinite(depth) & (depth >= coverage_min) & np.isfinite(dosage)
    x0 = np.column_stack([np.ones(mask.sum()), covariates[mask]])
    if mask.sum() <= x0.shape[1] + 1:
        return np.full(simulations, np.nan)
    fitted = np.clip(x0 @ np.linalg.lstsq(x0, methylation[mask], rcond=None)[0], 0.005, 0.995)
    trials = np.rint(depth[mask]).astype(int)
    design = np.column_stack([x0, dosage[mask]])
    rng = np.random.default_rng(seed)
    betas = np.empty(simulations)
    for index in range(simulations):
        simulated = rng.binomial(trials, fitted) / np.maximum(trials, 1)
        betas[index] = float(hc3_fit(direct_rint(simulated), design)["beta"])
    return betas


def analyze_chrom(args: argparse.Namespace) -> None:
    pairs = read_table(args.pairs)
    pairs = pairs[pairs.target_chrom.astype(str).eq(args.chrom)].copy()
    depth = read_table(args.depth).pivot(index="feature_id", columns="sample_id")
    design_data = np.load(args.design, allow_pickle=True)
    samples = design_data["sample_ids"].astype(str).tolist()
    covariates = np.asarray(design_data["covariates"], float)
    dosages = load_dosage_rows(pairs, samples)
    rows = []
    simulations = []
    for pair in pairs.itertuples(index=False):
        methylation = depth["methylation"].loc[str(pair.feature_id), samples].to_numpy(float)
        coverage = depth["depth"].loc[str(pair.feature_id), samples].to_numpy(float)
        dosage = dosages[(str(pair.sv_id), int(pair.sv_index))]
        fit = model_pair(methylation, coverage, dosage, covariates, args.coverage_min)
        rows.append({**pair._asdict(), **fit})
        if bool(pair.is_discovery) and fit.get("status") == "ok" and args.simulations > 0:
            betas = simulate_null_negative_fraction(
                methylation, coverage, dosage, covariates, args.coverage_min,
                args.simulations, args.seed + int(str(pair.pair_id).replace("tail", "")),
            )
            simulations.append({
                "pair_id": pair.pair_id, "n_simulations": args.simulations,
                "null_negative_fraction": float(np.mean(betas < 0)),
                "null_beta_median": float(np.nanmedian(betas)),
                "null_beta_q025": float(np.nanquantile(betas, 0.025)),
                "null_beta_q975": float(np.nanquantile(betas, 0.975)),
            })
    write_table(pd.DataFrame(rows), args.output)
    simulation_columns = ["pair_id", "n_simulations", "null_negative_fraction", "null_beta_median", "null_beta_q025", "null_beta_q975"]
    write_table(pd.DataFrame(simulations).reindex(columns=simulation_columns), args.simulation_output)


def bh_qvalues(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    q = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    if len(finite):
        order = finite[np.argsort(p[finite], kind="mergesort")]
        adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    return pd.Series(q, index=values.index)


def compartment_matched_permutation(
    annotations: pd.DataFrame,
    controls: pd.DataFrame,
    group_col: str,
    permutations: int = 10_000,
    seed: int = 20260831,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for group, hits in annotations.groupby(group_col, dropna=False):
        observed = []
        pools = []
        for hit in hits.itertuples(index=False):
            value = float(getattr(hit, "cpg_compartment"))
            local = controls[controls["matched_feature_id"].astype(str).eq(str(hit.feature_id))]
            pool = pd.to_numeric(local["cpg_compartment"], errors="coerce").dropna().to_numpy(float)
            pool = pool[pool != 0]
            if np.isfinite(value) and value != 0 and len(pool):
                observed.append(value > 0)
                pools.append(pool)
        if not pools:
            rows.append({group_col: group, "status": "not_estimable", "n_evaluable": 0})
            continue
        null = np.empty(permutations)
        for index in range(permutations):
            null[index] = np.mean([rng.choice(pool) > 0 for pool in pools])
        obs = float(np.mean(observed))
        rows.append({
            group_col: group, "status": "estimated", "n_evaluable": len(pools),
            "observed_A_fraction": obs, "matched_expected_A_fraction": float(null.mean()),
            "empirical_p_enrichment": float((1 + np.sum(null >= obs)) / (permutations + 1)),
            "n_permutations": permutations,
        })
    return pd.DataFrame(rows)


def summarize(args: argparse.Namespace) -> None:
    output = args.output_dir
    result_paths = sorted(output.glob("chromosomes/*.audit.tsv.gz"))
    results = pd.concat([read_table(path) for path in result_paths if path.stat().st_size > 50], ignore_index=True)
    simulation_paths = sorted(output.glob("chromosomes/*.simulation.tsv.gz"))
    simulations = pd.concat([read_table(path) for path in simulation_paths if path.stat().st_size > 80], ignore_index=True)
    for column in ["baseline_p", "depth_p", "threshold_pass_p", "depth_adjusted_p", "coverage_weighted_p", "count_aware_p"]:
        results[f"{column}_bh_q_tail"] = bh_qvalues(results[column])
    discoveries = results[results.is_discovery.astype(bool)].copy()
    material = (
        discoveries["depth_adjusted_effect_retention"].abs().median() < args.material_retention
        or (discoveries["depth_adjusted_p"] >= 0.05).mean() >= args.material_loss_fraction
    )
    write_table(results, output / "depth_audit_pair_results.tsv.gz")

    cpg_annotations = read_table(output / "cpg_censat_annotations.tsv.gz")
    discovery_features = discoveries[["feature_id", "source_pos", "source_end", "expected_raw_beta"]].drop_duplicates("feature_id")
    direction = discovery_features.merge(cpg_annotations, on="feature_id", how="left")
    direction["distance_to_sv_bp"] = np.maximum(
        direction["start"].to_numpy() - direction["source_end"].to_numpy(),
        direction["source_pos"].to_numpy() - direction["end"].to_numpy(),
    ).clip(min=0)
    direction["distance_class"] = pd.cut(
        direction["distance_to_sv_bp"], [-1, 100_000, 1_000_000, 5_000_000, np.inf],
        labels=["<=100 kb", "0.1-1 Mb", "1-5 Mb", ">5 Mb"],
    )
    direction_summary = (
        direction.groupby(["satellite_class", "distance_class"], observed=True)
        .agg(n_cpgs=("feature_id", "size"), negative_beta_fraction=("expected_raw_beta", lambda x: float((x < 0).mean())), median_beta=("expected_raw_beta", "median"))
        .reset_index()
    )
    write_table(direction_summary, output / "effect_direction_by_satellite_distance.tsv")

    summary = {
        "n_tail_pairs": int(len(results)), "n_discovery_pairs": int(len(discoveries)),
        "baseline_validation": {
            "median_abs_log10_p_error": float(np.nanmedian(np.abs(np.log10(results.baseline_p) - np.log10(results.expected_rint_p)))),
            "median_abs_beta_error": float(np.nanmedian(np.abs(results.baseline_beta - results.expected_rint_beta))),
        },
        "coverage_association": {
            "depth_bh_q_lt_0_05": int((results.depth_p_bh_q_tail < 0.05).sum()),
            "threshold_pass_bh_q_lt_0_05": int((results.threshold_pass_p_bh_q_tail < 0.05).sum()),
        },
        "discovery_robustness": {
            "median_depth_adjusted_effect_retention": float(discoveries.depth_adjusted_effect_retention.median()),
            "fraction_depth_adjusted_p_ge_0_05": float((discoveries.depth_adjusted_p >= 0.05).mean()),
            "fraction_coverage_weighted_p_lt_0_05": float((discoveries.coverage_weighted_p < 0.05).mean()),
            "fraction_count_aware_p_lt_0_05": float((discoveries.count_aware_p < 0.05).mean()),
        },
        "null_simulation": {
            "n_pairs": int(len(simulations)),
            "mean_null_negative_fraction": float(simulations.null_negative_fraction.mean()) if len(simulations) else None,
            "observed_negative_fraction": float((discoveries.expected_raw_beta < 0).mean()),
        },
        "depth_materially_explains_signal": bool(material),
        "decision": "rerun_complete_cis_family_depth_aware" if material else "retain_completed_PC10_discovery_layer",
        "materiality_rule": f"median |depth-adjusted/baseline beta| < {args.material_retention} or >= {args.material_loss_fraction:.0%} of discoveries have depth-adjusted p >= 0.05",
    }
    (output / "coverage_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_review_figures(output, discoveries, direction_summary)
    print(json.dumps(summary, indent=2))


def plot_review_figures(output: Path, discoveries: pd.DataFrame, direction_summary: pd.DataFrame) -> None:
    svs = read_table(output / "tested_sv_censat_annotations.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)
    counts = svs.satellite_class.value_counts().reindex(SATELLITE_ORDER, fill_value=0)
    axes[0].bar(np.arange(len(counts)), counts, color="#007A78", edgecolor="black", linewidth=0.5)
    axes[0].set_xticks(np.arange(len(counts)), [x.replace("_", "\n") for x in counts.index], fontsize=8)
    axes[0].set_ylabel("Tested physical SVs")
    axes[0].set_title("cenSat context of all tested MAC8 SVs")
    distances = pd.to_numeric(svs.register_distance_171bp, errors="coerce").dropna()
    axes[1].hist(distances, bins=np.arange(0, 87, 5), color="#D55E00", edgecolor="white")
    axes[1].axvline(10, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Distance to nearest 171-bp multiple (bp)")
    axes[1].set_ylabel("Tested physical SVs")
    axes[1].set_title("Alpha-satellite register periodicity audit")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"fig1_censat_periodicity.{suffix}", dpi=240 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)
    plot = direction_summary.copy()
    plot["label"] = plot.satellite_class.astype(str) + "\n" + plot.distance_class.astype(str)
    axes[0].bar(np.arange(len(plot)), plot.negative_beta_fraction, color="#2563A6")
    axes[0].axhline(0.5, color="black", linestyle="--", linewidth=1)
    axes[0].set_xticks(np.arange(len(plot)), plot.label, rotation=70, ha="right", fontsize=7)
    axes[0].set_ylabel("Fraction with negative raw beta")
    axes[0].set_title("Effect direction by satellite context and distance")
    axes[1].scatter(discoveries.baseline_beta, discoveries.depth_adjusted_beta, s=15, alpha=0.45, color="#007A78")
    lim = np.nanmax(np.abs(np.r_[discoveries.baseline_beta, discoveries.depth_adjusted_beta]))
    axes[1].plot([-lim, lim], [-lim, lim], color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Baseline PC10 RINT beta")
    axes[1].set_ylabel("Log-depth-adjusted PC10 RINT beta")
    axes[1].set_title("Coverage robustness of cis discoveries")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"fig2_satellite_direction_coverage.{suffix}", dpi=240 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("--pairs", type=Path, required=True)
    prepare_parser.add_argument("--cis-associations", type=Path, required=True)
    prepare_parser.add_argument("--source-jobs", type=Path, required=True)
    prepare_parser.add_argument("--matched-controls", type=Path, required=True)
    prepare_parser.add_argument("--tested-svs", type=Path, required=True)
    prepare_parser.add_argument("--censat", type=Path, required=True)
    prepare_parser.add_argument("--tail-p", type=float, default=1e-4)
    prepare_parser.add_argument("--discovery-q", type=float, default=0.01)
    prepare_parser.add_argument("--period-bp", type=int, default=171)
    prepare_parser.add_argument("--period-tolerance", type=int, default=10)
    prepare_parser.add_argument("--permutations", type=int, default=10_000)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)

    extract_parser = sub.add_parser("extract-depth")
    extract_parser.add_argument("--features", type=Path, required=True)
    extract_parser.add_argument("--methylation-manifest", type=Path, required=True)
    extract_parser.add_argument("--design", type=Path, required=True)
    extract_parser.add_argument("--chrom", required=True)
    extract_parser.add_argument("--output", type=Path, required=True)
    extract_parser.add_argument("--summary", type=Path, required=True)

    analyze_parser = sub.add_parser("analyze-chrom")
    analyze_parser.add_argument("--pairs", type=Path, required=True)
    analyze_parser.add_argument("--depth", type=Path, required=True)
    analyze_parser.add_argument("--design", type=Path, required=True)
    analyze_parser.add_argument("--chrom", required=True)
    analyze_parser.add_argument("--coverage-min", type=float, default=10)
    analyze_parser.add_argument("--simulations", type=int, default=200)
    analyze_parser.add_argument("--seed", type=int, default=20260829)
    analyze_parser.add_argument("--output", type=Path, required=True)
    analyze_parser.add_argument("--simulation-output", type=Path, required=True)

    summary_parser = sub.add_parser("summarize")
    summary_parser.add_argument("--output-dir", type=Path, required=True)
    summary_parser.add_argument("--material-retention", type=float, default=0.75)
    summary_parser.add_argument("--material-loss-fraction", type=float, default=0.20)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "prepare":
        prepare(args)
    elif args.command == "extract-depth":
        extract_depth(args)
    elif args.command == "analyze-chrom":
        analyze_chrom(args)
    else:
        summarize(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
