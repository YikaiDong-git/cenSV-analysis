#!/usr/bin/env python3
"""Statistical helpers for the robust meQTL multiplier-bootstrap scan."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Sequence

import numpy as np
from scipy import special, stats


METHOD_VERSION = "rint-hc3-score-multiplier-mac-v4"


def residual_basis(covariates: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    covariates = np.asarray(covariates, dtype=float)
    intercept = np.ones((covariates.shape[0], 1), dtype=float)
    design = np.column_stack([intercept, covariates]) if covariates.size else intercept
    u, singular_values, _ = np.linalg.svd(design, full_matrices=True)
    tol = np.finfo(float).eps * max(design.shape) * float(singular_values[0])
    rank = int(np.sum(singular_values > tol))
    return u[:, :rank], u[:, rank:], rank


def complete_case_dosage_frequency(
    predictor: np.ndarray, masks: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    predictor = np.asarray(predictor, dtype=float)
    masks = np.asarray(masks, dtype=bool)
    if predictor.ndim != 1 or masks.ndim != 2 or masks.shape[1] != len(predictor):
        raise ValueError("Predictor/mask dimensions are incompatible")
    finite = np.isfinite(predictor)
    if np.any(masks & ~finite[None, :]):
        raise ValueError("Complete-case masks include non-finite dosages")
    observed = predictor[finite]
    if np.any((observed < -1e-8) | (observed > 2.0 + 1e-8)):
        raise ValueError("Dosages must lie in [0, 2]")
    dosage = np.clip(np.where(finite, predictor, 0.0), 0.0, 2.0)
    n = masks.sum(axis=1).astype(float)
    alternate_count = masks @ dosage
    mac = np.minimum(alternate_count, 2.0 * n - alternate_count)
    with np.errstate(divide="ignore", invalid="ignore"):
        maf = mac / (2.0 * n)
    maf[n == 0] = np.nan
    return mac, maf


def row_rank_inverse_normal(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    ranks = stats.rankdata(values, axis=1, method="average", nan_policy="omit")
    n_observed = np.isfinite(values).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        transformed = special.ndtri((ranks - 0.5) / n_observed[:, None])
    transformed[~np.isfinite(values)] = np.nan
    return transformed


def canonical_sample_ids(
    covariate_ids: Sequence[str], sample_collections: Iterable[Sequence[str]]
) -> list[str]:
    covariate_ids = [str(sample) for sample in covariate_ids]
    if len(covariate_ids) != len(set(covariate_ids)):
        raise ValueError("Covariate sample IDs are not unique")
    common = set(covariate_ids)
    for samples in sample_collections:
        samples = [str(sample) for sample in samples]
        if len(samples) != len(set(samples)):
            raise ValueError("An input contains duplicate sample IDs")
        common.intersection_update(samples)
    if not common:
        raise ValueError("No samples are shared by all robust meQTL inputs")
    return sorted(common)


def _hash_strings(hasher: object, values: Sequence[str]) -> None:
    for value in values:
        encoded = str(value).encode("utf-8")
        hasher.update(len(encoded).to_bytes(8, "little"))
        hasher.update(encoded)


def sample_order_fingerprint(sample_ids: Sequence[str]) -> str:
    hasher = hashlib.sha256()
    _hash_strings(hasher, sample_ids)
    return hasher.hexdigest()


def replicate_seeds(base_seed: int, count: int) -> np.ndarray:
    rng = np.random.default_rng(base_seed)
    return rng.integers(1, np.iinfo(np.int64).max, size=count, dtype=np.int64)


def multiplier_operators(
    sample_ids: Sequence[str], base_seed: int, count: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    seeds = replicate_seeds(base_seed, count)
    n_samples = len(sample_ids)
    mammen = np.empty((count, n_samples), dtype=np.float32)
    rademacher = np.empty((count, n_samples), dtype=np.int8)
    sqrt_five = np.sqrt(5.0)
    low = np.float32((1.0 - sqrt_five) / 2.0)
    high = np.float32((1.0 + sqrt_five) / 2.0)
    probability_low = (sqrt_five + 1.0) / (2.0 * sqrt_five)
    for index, seed in enumerate(seeds):
        rng = np.random.default_rng(int(seed))
        mammen[index] = np.where(rng.random(n_samples) < probability_low, low, high)
        rademacher[index] = np.where(rng.random(n_samples) < 0.5, -1, 1)
    return seeds, mammen, rademacher


def array_fingerprint(*arrays: np.ndarray, labels: Sequence[str] = ()) -> str:
    hasher = hashlib.sha256()
    _hash_strings(hasher, [METHOD_VERSION, *labels])
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        _hash_strings(hasher, [str(contiguous.dtype), json.dumps(contiguous.shape)])
        hasher.update(contiguous.tobytes())
    return hasher.hexdigest()


def design_fingerprint(
    sample_ids: Sequence[str],
    covariate_names: Sequence[str],
    covariates: np.ndarray,
    seeds: np.ndarray,
    mammen: np.ndarray,
    rademacher: np.ndarray,
) -> str:
    labels = [*sample_ids, "COVARIATES", *covariate_names]
    return array_fingerprint(covariates, seeds, mammen, rademacher, labels=labels)


def batched_reduced_residuals(
    values: np.ndarray,
    covariates: np.ndarray,
    masks: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    masks = np.asarray(masks, dtype=bool)
    covariates = np.asarray(covariates, dtype=float)
    if values.shape != masks.shape:
        raise ValueError("values and masks must have identical feature-by-sample shape")
    if values.shape[1] != covariates.shape[0]:
        raise ValueError("covariates are not aligned to the value columns")
    design = np.column_stack([np.ones(len(covariates), dtype=float), covariates])
    weights = masks.astype(float)
    gram = np.einsum("bn,ni,nj->bij", weights, design, design, optimize=True)
    bread = np.linalg.pinv(gram, rcond=1e-12, hermitian=True)
    rank = np.linalg.matrix_rank(gram, hermitian=True)
    filled = np.where(masks, values, 0.0)
    rhs = np.einsum("bn,ni->bi", filled, design, optimize=True)
    coefficients = np.einsum("bij,bj->bi", bread, rhs, optimize=True)
    residuals = np.where(masks, filled - coefficients @ design.T, 0.0)
    leverage = np.einsum("ni,bij,nj->bn", design, bread, design, optimize=True)
    leverage = np.where(masks, leverage, 0.0)
    return residuals, leverage, bread, rank.astype(np.int16), design


def residualize_predictor_batched(
    predictor: np.ndarray,
    masks: np.ndarray,
    bread: np.ndarray,
    design: np.ndarray,
) -> np.ndarray:
    predictor = np.asarray(predictor, dtype=float)
    masks = np.asarray(masks, dtype=bool)
    rhs = np.einsum("bn,n,ni->bi", masks, predictor, design, optimize=True)
    coefficients = np.einsum("bij,bj->bi", bread, rhs, optimize=True)
    return np.where(masks, predictor[None, :] - coefficients @ design.T, 0.0)


def fwl_hc3_statistics(
    predictor_residuals: np.ndarray,
    outcome_residuals: np.ndarray,
    reduced_leverage: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(predictor_residuals, dtype=float)
    y = np.asarray(outcome_residuals, dtype=float)
    x_ss = np.sum(x * x, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        beta = np.sum(x * y, axis=1) / x_ss
        residual = y - beta[:, None] * x
        leverage = reduced_leverage + (x * x) / x_ss[:, None]
        influence = x * residual / np.clip(1.0 - leverage, 1e-8, None)
        se = np.sqrt(np.sum(influence * influence, axis=1)) / x_ss
        t_stat = beta / se
    beta[~np.isfinite(beta)] = np.nan
    se[~np.isfinite(se)] = np.nan
    t_stat[~np.isfinite(t_stat)] = np.nan
    return beta, se, t_stat


def hc3_null_score_influence(
    predictor_residuals: np.ndarray,
    reduced_outcome_residuals: np.ndarray,
    reduced_leverage: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    raw_score = predictor_residuals * reduced_outcome_residuals
    influence = raw_score / np.clip(1.0 - reduced_leverage, 1e-8, None)
    scale = np.sqrt(np.sum(influence * influence, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = np.sum(raw_score, axis=1) / scale
    statistic[~np.isfinite(statistic)] = np.nan
    return influence, statistic


def multiplier_score_statistics(
    multipliers: np.ndarray, influence: np.ndarray
) -> np.ndarray:
    influence = np.asarray(influence, dtype=float)
    scale = np.sqrt(np.sum(influence * influence, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        statistics = np.asarray(multipliers, dtype=float) @ influence.T
        statistics /= scale[None, :]
    statistics[~np.isfinite(statistics)] = 0.0
    return statistics


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not len(values):
        return values.copy()
    if not np.isfinite(values).all():
        raise ValueError("BH input contains non-finite p-values")
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1]
    output = np.empty(len(values), dtype=float)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def residualize_outcome_batched(
    values: np.ndarray,
    masks: np.ndarray,
    bread: np.ndarray,
    design: np.ndarray,
) -> np.ndarray:
    filled = np.where(masks, values, 0.0)
    rhs = np.einsum("bn,ni->bi", filled, design, optimize=True)
    coefficients = np.einsum("bij,bj->bi", bread, rhs, optimize=True)
    return np.where(masks, filled - coefficients @ design.T, 0.0)
