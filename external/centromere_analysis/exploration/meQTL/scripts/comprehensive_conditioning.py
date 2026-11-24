#!/usr/bin/env python3
"""Shared deterministic helpers for comprehensive cis SNP conditioning."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import special, stats
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

METHOD_VERSION = "comprehensive-cis-conditioning-v1"
PAIR_KEYS = ["sv_id", "feature_id"]
PC_COUNTS = (5, 10, 20, 30)


def sha256_file(path: str | Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_frame_hash(frame: pd.DataFrame, columns: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for row in frame[list(columns)].itertuples(index=False, name=None):
        payload = json.dumps([str(value) for value in row], separators=(",", ":"))
        digest.update(len(payload).to_bytes(8, "little"))
        digest.update(payload.encode("utf-8"))
    return digest.hexdigest()


def build_pair_union(hierarchical: pd.DataFrame, global_pairs: pd.DataFrame) -> pd.DataFrame:
    for name, frame in (("hierarchical", hierarchical), ("global", global_pairs)):
        missing = set(PAIR_KEYS) - set(frame)
        if missing:
            raise ValueError(f"{name} discoveries lack keys: {sorted(missing)}")
        if frame.duplicated(PAIR_KEYS).any():
            raise ValueError(f"{name} discoveries contain duplicate pair identities")
    h = hierarchical.copy()
    h["discovery_scope"] = "hierarchical_primary"
    h["hierarchical_primary"] = True
    h["global_pair_discovery"] = h.get("global_pair_discovery", False).fillna(False).astype(bool)
    h_keys = pd.MultiIndex.from_frame(h[PAIR_KEYS])
    g = global_pairs.loc[~pd.MultiIndex.from_frame(global_pairs[PAIR_KEYS]).isin(h_keys)].copy()
    g["discovery_scope"] = "global_only_secondary"
    g["hierarchical_primary"] = False
    g["global_pair_discovery"] = True
    union = pd.concat([h, g], ignore_index=True, sort=False)
    union = union.sort_values(["canonical_pair_index", *PAIR_KEYS], kind="stable").reset_index(drop=True)
    union.insert(0, "pair_id", [f"ccpair{i:06d}" for i in range(len(union))])
    if union.duplicated(PAIR_KEYS).any():
        raise AssertionError("Pair union is not unique")
    return union


def make_batches(pairs: pd.DataFrame, maximum: int = 10) -> pd.DataFrame:
    if maximum < 1 or maximum > 10:
        raise ValueError("Outcome batch size must be between 1 and 10")
    rows: list[dict[str, object]] = []
    batch_index = 0
    ordered = pairs.sort_values(["sv_id", "canonical_pair_index", "feature_id"], kind="stable")
    for sv_id, group in ordered.groupby("sv_id", sort=False):
        ids = group["pair_id"].astype(str).tolist()
        for offset in range(0, len(ids), maximum):
            selected = ids[offset : offset + maximum]
            rows.append({
                "outcome_batch_id": f"obatch{batch_index:05d}",
                "sv_id": str(sv_id),
                "pair_ids": ",".join(selected),
                "n_outcomes": len(selected),
            })
            batch_index += 1
    return pd.DataFrame(rows)


def expand_analysis_tasks(
    batches: pd.DataFrame, panels: Iterable[str], windows: Iterable[int]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    svs = batches["sv_id"].drop_duplicates().tolist()
    snp_rows: list[dict[str, object]] = []
    snp_lookup: dict[tuple[str, str, int], str] = {}
    for sv_id in svs:
        for panel in panels:
            for window in windows:
                task_id = f"snp{len(snp_rows):05d}"
                snp_lookup[(str(sv_id), str(panel), int(window))] = task_id
                snp_rows.append({"snp_task_id": task_id, "sv_id": sv_id, "panel": panel, "window_bp": int(window)})
    analysis_rows: list[dict[str, object]] = []
    for batch in batches.itertuples(index=False):
        for panel in panels:
            for window in windows:
                analysis_rows.append({
                    "analysis_task_id": f"analysis{len(analysis_rows):05d}",
                    "outcome_batch_id": batch.outcome_batch_id,
                    "snp_task_id": snp_lookup[(str(batch.sv_id), str(panel), int(window))],
                    "sv_id": batch.sv_id,
                    "pair_ids": batch.pair_ids,
                    "n_outcomes": int(batch.n_outcomes),
                    "panel": panel,
                    "window_bp": int(window),
                })
    return pd.DataFrame(snp_rows), pd.DataFrame(analysis_rows)


def direct_rint(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("direct_rint requires a finite one-dimensional complete-case vector")
    ranks = stats.rankdata(values, method="average")
    return special.ndtri((ranks - 0.5) / len(values))


def hc3_fit(y: np.ndarray, design: np.ndarray, coefficient: int = -1) -> dict[str, float | int]:
    y = np.asarray(y, dtype=float)
    x = np.asarray(design, dtype=float)
    if y.ndim != 1 or x.ndim != 2 or len(y) != len(x) or not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError("HC3 inputs must be finite and aligned")
    rank = int(np.linalg.matrix_rank(x))
    df = int(len(y) - rank)
    if rank < x.shape[1] or df <= 0:
        return {"status": "rank_deficient", "beta": np.nan, "hc3_se": np.nan, "hc3_t": np.nan,
                "hc3_p": np.nan, "rank": rank, "df": df, "max_leverage": np.nan,
                "mean_leverage": np.nan, "condition_number": np.inf}
    bread = np.linalg.inv(x.T @ x)
    beta = bread @ x.T @ y
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, bread, x)
    adjusted = residual / np.clip(1.0 - leverage, 1e-8, None)
    meat = (x * adjusted[:, None]).T @ (x * adjusted[:, None])
    covariance = bread @ meat @ bread
    index = coefficient if coefficient >= 0 else x.shape[1] + coefficient
    se = float(np.sqrt(max(0.0, covariance[index, index])))
    t_value = float(beta[index] / se) if se > 0 else np.nan
    return {
        "status": "ok" if np.isfinite(t_value) else "nonfinite",
        "beta": float(beta[index]), "hc3_se": se, "hc3_t": t_value,
        "hc3_p": float(2 * stats.t.sf(abs(t_value), df)) if np.isfinite(t_value) else np.nan,
        "rank": rank, "df": df, "max_leverage": float(leverage.max()),
        "mean_leverage": float(leverage.mean()), "condition_number": float(np.linalg.cond(x)),
    }


def residualize_train_test(train: np.ndarray, test: np.ndarray, cov_train: np.ndarray, cov_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_train = np.column_stack([np.ones(len(train)), cov_train])
    x_test = np.column_stack([np.ones(len(test)), cov_test])
    coefficient = np.linalg.lstsq(x_train, train, rcond=None)[0]
    return train - x_train @ coefficient, test - x_test @ coefficient


def residualize_matrix(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(covariates)), covariates])
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def abs_correlations(genotypes: np.ndarray, sv: np.ndarray) -> np.ndarray:
    numerator = np.abs(genotypes.T @ sv)
    denominator = np.sqrt(np.sum(genotypes * genotypes, axis=0) * np.sum(sv * sv))
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def phenotype_independent_snp_selection(
    genotypes: np.ndarray,
    sv: np.ndarray,
    covariates: np.ndarray,
    maximum: int,
    variant_ids: Iterable[str],
) -> tuple[np.ndarray, np.ndarray]:
    if maximum < 1:
        return np.array([], dtype=int), np.array([], dtype=float)
    g_resid = residualize_matrix(np.asarray(genotypes, float), covariates)
    sv_resid = residualize_matrix(np.asarray(sv, float)[:, None], covariates)[:, 0]
    scores = abs_correlations(g_resid, sv_resid)
    ids = np.asarray([str(value) for value in variant_ids])
    order = np.lexsort((ids, -scores))
    return order[: min(maximum, len(order))], scores


def ld_prune(indices: np.ndarray, genotypes: np.ndarray, maximum: int, r2_threshold: float) -> np.ndarray:
    chosen: list[int] = []
    standardized = np.asarray(genotypes, float)
    standardized = (standardized - standardized.mean(axis=0)) / np.where(standardized.std(axis=0) > 0, standardized.std(axis=0), 1)
    for index in indices:
        if len(chosen) >= maximum:
            break
        if not chosen:
            chosen.append(int(index))
            continue
        correlations = standardized[:, chosen].T @ standardized[:, int(index)] / len(standardized)
        if np.all(correlations * correlations < r2_threshold):
            chosen.append(int(index))
    return np.asarray(chosen, dtype=int)


def crossfit_residual_sv(
    sv: np.ndarray, genotypes: np.ndarray, covariates: np.ndarray, folds: int = 5, seed: int = 104729
) -> tuple[np.ndarray, float]:
    if genotypes.shape[1] == 0:
        return np.asarray(sv, float).copy(), np.nan
    residual = np.empty(len(sv), dtype=float)
    observed = np.empty(len(sv), dtype=float)
    predicted = np.empty(len(sv), dtype=float)
    splitter = KFold(n_splits=min(folds, len(sv)), shuffle=True, random_state=seed)
    for train, test in splitter.split(genotypes):
        g_train, g_test = residualize_train_test(genotypes[train], genotypes[test], covariates[train], covariates[test])
        sv_train, sv_test = residualize_train_test(sv[train, None], sv[test, None], covariates[train], covariates[test])
        scale = np.where(g_train.std(axis=0) > 0, g_train.std(axis=0), 1)
        model = RidgeCV(alphas=np.logspace(-3, 4, 16)).fit(g_train / scale, sv_train[:, 0])
        fold_pred = model.predict(g_test / scale)
        residual[test] = sv_test[:, 0] - fold_pred
        observed[test] = sv_test[:, 0]
        predicted[test] = fold_pred
    denominator = np.sum((observed - observed.mean()) ** 2)
    r2 = 1 - np.sum((observed - predicted) ** 2) / denominator if denominator > 0 else np.nan
    return residual, float(r2)


def interval_covered(path: str | Path | None, chrom: str, start: int, end: int) -> bool:
    if path is None:
        return True
    bed = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    intervals = bed.loc[bed.chrom.astype(str).eq(chrom), ["start", "end"]].sort_values("start")
    cursor = int(start)
    for row in intervals.itertuples(index=False):
        if int(row.end) <= cursor:
            continue
        if int(row.start) > cursor:
            return False
        cursor = max(cursor, int(row.end))
        if cursor >= end:
            return True
    return False
