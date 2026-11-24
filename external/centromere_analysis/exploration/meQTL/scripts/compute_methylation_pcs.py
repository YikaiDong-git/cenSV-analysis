#!/usr/bin/env python3
"""Compute methylation sample PCs from chunk matrices with bounded memory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


MATRIX_META_COLUMNS = {"chrom", "start", "end", "feature_id", "n_observed"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices", nargs="+", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--n-pcs", type=int, default=5)
    parser.add_argument("--max-sites", type=int, default=50000)
    parser.add_argument("--chunksize", type=int, default=20000)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def matrix_summary_path(matrix_path: Path) -> Path:
    name = matrix_path.name
    if name.endswith(".methylation.tsv.gz"):
        return matrix_path.with_name(name.replace(".methylation.tsv.gz", ".summary.json"))
    return matrix_path.with_suffix(".summary.json")


def count_matrix_rows(matrix_path: Path) -> int:
    summary_path = matrix_summary_path(matrix_path)
    if summary_path.exists():
        summary = json.loads(summary_path.read_text())
        if "n_output_sites" in summary:
            return int(summary["n_output_sites"])
    with pd.read_csv(matrix_path, sep="\t", chunksize=100000) as reader:
        return sum(len(chunk) for chunk in reader)


def select_global_indices(total_sites: int, max_sites: int) -> np.ndarray:
    if total_sites <= 0:
        return np.array([], dtype=np.int64)
    if max_sites <= 0 or total_sites <= max_sites:
        return np.arange(total_sites, dtype=np.int64)
    return np.unique(np.linspace(0, total_sites - 1, num=max_sites, dtype=np.int64))


def standardize_site_rows(values: np.ndarray) -> np.ndarray:
    means = np.nanmean(values, axis=1)
    keep = ~np.isnan(means)
    if not np.any(keep):
        return np.empty((0, values.shape[1]), dtype=float)
    values = values[keep].copy()
    means = means[keep]
    missing = np.isnan(values)
    if np.any(missing):
        values[missing] = means[np.where(missing)[0]]
    values -= values.mean(axis=1, keepdims=True)
    std = values.std(axis=1, ddof=0)
    keep = std > 0
    if not np.any(keep):
        return np.empty((0, values.shape[1]), dtype=float)
    return values[keep] / std[keep, None]


def accumulate_sample_gram(
    matrices: list[Path],
    samples: list[str],
    selected_indices: np.ndarray,
    matrix_counts: list[int],
    chunksize: int,
) -> tuple[np.ndarray, int]:
    sample_count = len(samples)
    gram = np.zeros((sample_count, sample_count), dtype=np.float64)
    n_sites_used = 0
    offset = 0

    for matrix_path, matrix_rows in zip(matrices, matrix_counts, strict=True):
        local_keep = selected_indices[
            (selected_indices >= offset) & (selected_indices < offset + matrix_rows)
        ] - offset
        offset += matrix_rows
        if local_keep.size == 0:
            continue
        local_keep_set = set(local_keep.astype(int).tolist())
        row_start = 0
        for chunk in pd.read_csv(matrix_path, sep="\t", chunksize=chunksize):
            row_end = row_start + len(chunk)
            wanted = [idx - row_start for idx in local_keep_set if row_start <= idx < row_end]
            if wanted:
                missing_samples = [sample for sample in samples if sample not in chunk.columns]
                if missing_samples:
                    raise ValueError(f"{matrix_path} is missing {len(missing_samples)} sample columns")
                values = chunk.iloc[sorted(wanted)][samples].to_numpy(dtype=float)
                standardized = standardize_site_rows(values)
                if standardized.size:
                    gram += standardized.T @ standardized
                    n_sites_used += standardized.shape[0]
            row_start = row_end
    return gram, n_sites_used


def gram_to_pc_scores(
    gram: np.ndarray,
    samples: list[str],
    n_pcs: int,
) -> tuple[pd.DataFrame, list[float], list[float]]:
    output = pd.DataFrame({"sample_id": samples})
    if n_pcs <= 0 or gram.size == 0 or not np.any(gram):
        return output, [], []

    gram = (gram + gram.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > max(float(eigenvalues[0]) * 1e-12, 0.0)
    n_components = min(n_pcs, len(samples) - 1, int(np.sum(positive)))
    if n_components <= 0:
        return output, [], []

    scores = eigenvectors[:, :n_components] * np.sqrt(eigenvalues[:n_components])
    for pc_index in range(n_components):
        anchor = int(np.argmax(np.abs(scores[:, pc_index])))
        if scores[anchor, pc_index] < 0:
            scores[:, pc_index] *= -1

    for pc_index in range(n_components):
        output[f"methylPC{pc_index + 1}"] = scores[:, pc_index]

    total_variance = float(np.sum(eigenvalues))
    ratios = (
        (eigenvalues[:n_components] / total_variance).astype(float).tolist()
        if total_variance > 0
        else []
    )
    explained = eigenvalues[:n_components].astype(float).tolist()
    return output, explained, ratios


def main() -> int:
    args = parse_args()
    samples = read_samples(args.samples)
    matrix_counts = [count_matrix_rows(path) for path in args.matrices]
    total_sites = int(sum(matrix_counts))
    selected_indices = select_global_indices(total_sites, args.max_sites)
    gram, n_sites_used = accumulate_sample_gram(
        args.matrices,
        samples,
        selected_indices,
        matrix_counts,
        args.chunksize,
    )
    pcs, explained, ratios = gram_to_pc_scores(gram, samples, args.n_pcs)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pcs.to_csv(args.output, sep="\t", index=False)
    summary = {
        "method": "standardized_site_gram_eigendecomposition",
        "n_samples": len(samples),
        "n_matrices": len(args.matrices),
        "matrix_rows": dict(zip([str(path) for path in args.matrices], matrix_counts, strict=True)),
        "n_sites_available": total_sites,
        "n_sites_selected": int(len(selected_indices)),
        "n_sites_used": int(n_sites_used),
        "n_pcs_requested": int(args.n_pcs),
        "n_pcs": int(len([col for col in pcs.columns if col != "sample_id"])),
        "max_sites": int(args.max_sites),
        "chunksize": int(args.chunksize),
        "explained_variance": explained,
        "explained_variance_ratio": ratios,
        "output": str(args.output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
