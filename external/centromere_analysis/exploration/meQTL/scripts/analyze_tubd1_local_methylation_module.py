#!/usr/bin/env python3
"""Analyze local methylation modules across the VMP1/MIR21/TUBD1 locus."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats
from sklearn.utils.extmath import randomized_svd

from run_sv_methylation_assoc import MATRIX_META_COLUMNS, residualize


GENES = {
    "VMP1": (60_575_940, 60_711_009, "+"),
    "MIR21": (60_710_019, 60_710_091, "+"),
    "TUBD1": (60_728_230, 60_761_693, "-"),
}
PLOT_GENES = {
    "YPEL2": (60_199_794, 60_269_870, "+"),
    **GENES,
}
LOCUS_START = min(value[0] for value in GENES.values())
LOCUS_END = max(value[1] for value in GENES.values())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methylation-matrix", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--sv-dosages", type=Path, required=True)
    parser.add_argument("--sv-index", type=int, default=3)
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--expression-covariates", type=Path, required=True)
    parser.add_argument("--screen-ccres", type=Path, required=True)
    parser.add_argument("--gm12878-ccres", type=Path, required=True)
    parser.add_argument("--monocyte-ccres", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--components", type=int, default=20)
    parser.add_argument("--minimum-observed-fraction", type=float, default=0.8)
    parser.add_argument("--minimum-sd", type=float, default=0.02)
    parser.add_argument("--component-permutations", type=int, default=10_000)
    parser.add_argument("--cpg-permutations", type=int, default=1_000)
    parser.add_argument("--seed", type=int, default=1_047_291)
    return parser.parse_args()


def bh_qvalues(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(values), dtype=float)
    result[order] = np.minimum(ranked, 1.0)
    return result


def association_statistics(
    x_resid: np.ndarray, y_resid: np.ndarray, covariate_rank: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if x_resid.ndim == 1:
        x_resid = x_resid[:, None]
    if y_resid.ndim == 1:
        y_resid = y_resid[:, None]
    x_ss = np.sum(x_resid * x_resid, axis=0)
    y_ss = np.sum(y_resid * y_resid, axis=0)
    xy = x_resid.T @ y_resid
    beta = xy / x_ss[:, None]
    df = len(x_resid) - covariate_rank - 1
    rss = np.maximum(y_ss[None, :] - beta * beta * x_ss[:, None], 0.0)
    se = np.sqrt((rss / df) / x_ss[:, None])
    t_value = beta / se
    pvalue = 2.0 * stats.t.sf(np.abs(t_value), df=df)
    return beta, se, t_value, pvalue


def permuted_max_statistics(
    dosage: np.ndarray,
    outcomes_residual: np.ndarray,
    covariates: np.ndarray,
    covariate_rank: int,
    permutations: int,
    rng: np.random.Generator,
    batch_size: int = 100,
) -> np.ndarray:
    maxima: list[np.ndarray] = []
    for batch_start in range(0, permutations, batch_size):
        count = min(batch_size, permutations - batch_start)
        shuffled = np.column_stack(
            [dosage[rng.permutation(len(dosage))] for _ in range(count)]
        )
        shuffled_residual, _ = residualize(shuffled, covariates)
        _, _, t_value, _ = association_statistics(
            shuffled_residual, outcomes_residual, covariate_rank
        )
        maxima.append(np.max(np.abs(t_value), axis=1))
    return np.concatenate(maxima) if maxima else np.empty(0, dtype=float)


def load_locus_matrix(path: Path, start: int, end: int) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=25_000):
        mask = chunk["chrom"].eq("chr17") & chunk["start"].between(start, end)
        if mask.any():
            selected.append(chunk.loc[mask].copy())
    if not selected:
        raise ValueError(f"No CpGs found in chr17:{start}-{end}")
    return pd.concat(selected, ignore_index=True).sort_values("start").reset_index(drop=True)


def prepare_methylation(
    table: pd.DataFrame,
    samples: list[str],
    minimum_observed_fraction: float,
    minimum_sd: float,
) -> tuple[pd.DataFrame, np.ndarray]:
    values = table[samples].to_numpy(dtype=float)
    observed = np.sum(np.isfinite(values), axis=1)
    keep = observed >= int(np.ceil(minimum_observed_fraction * len(samples)))
    table = table.loc[keep].reset_index(drop=True)
    values = values[keep]
    observed = observed[keep]
    means = np.nansum(values, axis=1) / observed
    rows, columns = np.where(~np.isfinite(values))
    values[rows, columns] = means[rows]
    keep = np.std(values, axis=1, ddof=1) >= minimum_sd
    return table.loc[keep].reset_index(drop=True), values[keep]


def load_intervals(path: Path, start: int, end: int) -> pd.DataFrame:
    compression = "gzip" if str(path).endswith(".gz") else None
    table = pd.read_csv(path, sep="\t", header=None, compression=compression, comment="#")
    table = table.iloc[:, : max(3, min(table.shape[1], 6))]
    names = ["chrom", "start", "end", "id", "accession", "class"][: table.shape[1]]
    table.columns = names
    return table[
        table["chrom"].eq("chr17")
        & (table["start"] < end)
        & (table["end"] > start)
    ].sort_values("start").reset_index(drop=True)


def annotate_interval_overlap(positions: np.ndarray, intervals: pd.DataFrame) -> np.ndarray:
    if intervals.empty:
        return np.zeros(len(positions), dtype=bool)
    starts = intervals["start"].to_numpy(dtype=int)
    ends = intervals["end"].to_numpy(dtype=int)
    indices = np.searchsorted(starts, positions, side="right") - 1
    valid = indices >= 0
    result = np.zeros(len(positions), dtype=bool)
    result[valid] = positions[valid] < ends[indices[valid]]
    return result


def orient_components(scores: np.ndarray, loadings: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    signs = np.sign(np.sum(loadings, axis=1))
    signs[signs == 0] = 1
    return scores * signs[None, :], loadings * signs[:, None]


def make_predefined_modules(
    methylation: pd.DataFrame,
    residuals: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    definitions = {
        "VMP1_gene_body": methylation["start"].between(*GENES["VMP1"][:2]),
        "VMP1_promoter_2kb": methylation["start"].between(
            GENES["VMP1"][0] - 2_000, GENES["VMP1"][0] + 2_000
        ),
        "MIR21_plus_20kb": methylation["start"].between(
            GENES["MIR21"][0] - 20_000, GENES["MIR21"][1] + 20_000
        ),
        "TUBD1_gene_body": methylation["start"].between(*GENES["TUBD1"][:2]),
        "TUBD1_promoter_2kb": methylation["start"].between(
            GENES["TUBD1"][1] - 2_000, GENES["TUBD1"][1] + 2_000
        ),
        "SCREEN_cCRE_primary": methylation["screen_ccre"],
        "GM12878_active_cCRE_primary": methylation["gm12878_active"],
        "monocyte_active_cCRE_primary": methylation["monocyte_active"],
    }
    score_columns: dict[str, np.ndarray] = {}
    metadata_rows: list[dict] = []
    for name, mask_series in definitions.items():
        mask = mask_series.to_numpy(dtype=bool, copy=True)
        if name.endswith("_primary"):
            mask &= methylation["start"].between(
                LOCUS_START - 500_000, LOCUS_END + 500_000
            ).to_numpy()
        if not mask.any():
            continue
        matrix = residuals[:, mask]
        scale = np.std(matrix, axis=0, ddof=1)
        usable = scale > np.finfo(float).eps
        if not usable.any():
            continue
        standardized = matrix[:, usable] / scale[usable]
        score_columns[name] = np.mean(standardized, axis=1)
        selected = methylation.loc[mask].iloc[np.where(usable)[0]]
        metadata_rows.append(
            {
                "module": name,
                "n_cpgs": int(usable.sum()),
                "start": int(selected["start"].min()),
                "end": int(selected["end"].max()),
            }
        )
    return pd.DataFrame(score_columns), pd.DataFrame(metadata_rows)


def fit_ols(y: np.ndarray, predictors: np.ndarray) -> dict[str, np.ndarray | int]:
    design = np.column_stack([np.ones(len(y)), predictors])
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ coefficients
    rank = np.linalg.matrix_rank(design)
    df = len(y) - rank
    inverse = np.linalg.pinv(design.T @ design)
    variance = residual @ residual / df
    se = np.sqrt(np.maximum(np.diag(inverse) * variance, 0.0))
    t_value = coefficients / se
    pvalue = 2.0 * stats.t.sf(np.abs(t_value), df=df)
    return {"beta": coefficients, "se": se, "pvalue": pvalue, "df": int(df)}


def load_expression_rows(path: Path, targets: dict[str, str]) -> dict[str, pd.Series]:
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        result: dict[str, pd.Series] = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            base_id = fields[3].split(".")[0]
            if base_id in targets:
                result[targets[base_id]] = pd.Series(
                    np.asarray(fields[4:], dtype=float), index=header[4:]
                )
    return result


def load_expression_covariates(path: Path) -> tuple[pd.DataFrame, list[str], list[str]]:
    table = pd.read_csv(path, sep="\t", compression="gzip", index_col=0).astype(object)
    table.loc["sex"] = (
        table.loc["sex"].map({"XY": 1.0, "XX": 0.0, "male": 1.0, "female": 0.0})
        .fillna(0.0)
        .to_numpy()
    )
    genotype_pcs = [value for value in table.index if re.match(r"^PC\d+$", str(value))]
    peer_factors = [value for value in table.index if re.match(r"^PEER\d+$", str(value))]
    table = table.apply(pd.to_numeric, errors="coerce").fillna(0.0).T
    return table, genotype_pcs, peer_factors


def plot_results(
    cpg: pd.DataFrame,
    components: pd.DataFrame,
    significant_scores: pd.DataFrame,
    dosage: pd.Series,
    expression_plot: pd.DataFrame,
    screen_ccres: pd.DataFrame,
    gm12878_ccres: pd.DataFrame,
    expression_pvalue: float,
    output_dir: Path,
) -> None:
    fig = plt.figure(figsize=(14, 13))
    grid = fig.add_gridspec(
        4, 2, height_ratios=[1.65, 2.8, 2.2, 2.2], hspace=0.62, wspace=0.30
    )

    ax = fig.add_subplot(grid[0, :])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(
        "A  Evidence map: centromeric sequence to a centrosomal candidate",
        loc="left",
        fontweight="bold",
    )
    boxes = [
        (0.09, 0.70, "chr3 alpha-satellite HOR\n336-bp deletion; CENP-B box lost", "#dce9f2"),
        (0.45, 0.70, "Two corrected chr17 CpGs\nVMP1 intron; YPEL2 upstream", "#e1f0e8"),
        (0.70, 0.70, "TUBD1 expression\nnominal and PEER-sensitive", "#fff1d6"),
        (0.91, 0.70, "TUBD1 / delta-tubulin\ncentriole triplets and stability", "#dce9f2"),
        (0.17, 0.20, "Predicted loss of\nCENP-B occupancy", "#fff1d6"),
        (0.41, 0.20, "Candidate centromere state\nchromatin / kinetochore /\nnuclear organization", "#fff1d6"),
    ]
    for x_position, y_position, label, color in boxes:
        ax.text(
            x_position,
            y_position,
            label,
            ha="center",
            va="center",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": color, "edgecolor": "#54616b"},
        )
    solid_arrow = {"arrowstyle": "-|>", "color": "#315f9b", "linewidth": 1.5}
    dashed_arrow = {
        "arrowstyle": "-|>",
        "color": "#a86614",
        "linewidth": 1.3,
        "linestyle": "--",
    }
    ax.annotate("", xy=(0.36, 0.70), xytext=(0.19, 0.70), arrowprops=solid_arrow)
    ax.text(0.275, 0.77, "measured trans-meQTL", ha="center", fontsize=8, color="#315f9b")
    ax.annotate("", xy=(0.61, 0.70), xytext=(0.54, 0.70), arrowprops=dashed_arrow)
    ax.text(0.575, 0.77, "provisional", ha="center", fontsize=8, color="#a86614")
    ax.annotate("", xy=(0.83, 0.70), xytext=(0.78, 0.70), arrowprops=solid_arrow)
    ax.text(0.805, 0.77, "established biology", ha="center", fontsize=8, color="#315f9b")
    ax.annotate("", xy=(0.15, 0.33), xytext=(0.10, 0.56), arrowprops=dashed_arrow)
    ax.annotate("", xy=(0.31, 0.20), xytext=(0.24, 0.20), arrowprops=dashed_arrow)
    ax.annotate("", xy=(0.43, 0.55), xytext=(0.42, 0.33), arrowprops=dashed_arrow)
    ax.text(
        0.61,
        0.20,
        "Solid blue = measured association or established biology\n"
        "Dashed amber = unproven mechanistic step",
        ha="left",
        va="center",
        fontsize=8,
    )

    ax = fig.add_subplot(grid[1, :])
    colors = np.where(cpg["gm12878_active"], "#c43c39", np.where(cpg["screen_ccre"], "#d6a21d", "#637381"))
    ax.scatter(cpg["start"] / 1e6, -np.log10(cpg["pvalue"].clip(lower=1e-300)), c=colors, s=8, alpha=0.7, linewidths=0)
    ax.axhline(-np.log10(0.01 / len(cpg)), color="#333333", linestyle="--", linewidth=0.8)
    ax.set_ylabel(r"SV-meQTL $-\log_{10}(p)$")
    ax.set_title("B  chr3 centromeric deletion and local chr17 methylation", loc="left", fontweight="bold")
    ax.set_xlim(cpg["start"].min() / 1e6, cpg["start"].max() / 1e6)

    track = ax.inset_axes([0, -0.35, 1, 0.25], transform=ax.transAxes)
    track.set_xlim(ax.get_xlim())
    track.set_ylim(0, 1)
    for row in screen_ccres.itertuples(index=False):
        track.add_patch(Rectangle((row.start / 1e6, 0.76), (row.end - row.start) / 1e6, 0.07, color="#d6a21d", alpha=0.55))
    for row in gm12878_ccres.itertuples(index=False):
        track.add_patch(Rectangle((row.start / 1e6, 0.88), (row.end - row.start) / 1e6, 0.07, color="#c43c39", alpha=0.75))
    for index, (gene, (start, end, strand)) in enumerate(PLOT_GENES.items()):
        y = 0.08 + 0.14 * (index % 2)
        track.add_patch(Rectangle((start / 1e6, y), max((end - start) / 1e6, 0.002), 0.09, color="#2b6f77"))
        track.text((start + end) / 2e6, y - 0.015, gene, ha="center", va="top", fontsize=8, fontweight="bold")
    track.axis("off")

    annotation_legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#637381", markeredgecolor="none", markersize=6, label="Other CpG"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#d6a21d", markeredgecolor="none", markersize=6, label="SCREEN cCRE"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#c43c39", markeredgecolor="none", markersize=6, label="GM12878-active cCRE"),
    ]
    ax.legend(handles=annotation_legend, frameon=False, fontsize=8, loc="upper right", ncol=3)

    significant = cpg[cpg["qvalue_bh"] < 0.01].sort_values("pvalue")
    labels = ["VMP1 intron", "YPEL2 upstream"]
    for row, label in zip(significant.itertuples(index=False), labels):
        ax.annotate(
            label,
            xy=(row.start / 1e6, -np.log10(row.pvalue)),
            xytext=(0, 12),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            arrowprops={"arrowstyle": "-", "color": "#333333", "linewidth": 0.7},
        )

    ax = fig.add_subplot(grid[2, 0])
    palette = ["#2b6f77", "#c43c39"]
    for cpg_index, column in enumerate(significant_scores.columns):
        for dose in sorted(dosage.unique()):
            mask = dosage.eq(dose)
            jitter = np.linspace(-0.055, 0.055, int(mask.sum()))
            position = dose + (cpg_index - 0.5) * 0.18
            ax.scatter(
                position + jitter,
                significant_scores.loc[mask, column],
                s=16,
                alpha=0.55,
                color=palette[cpg_index],
                label=labels[cpg_index] if dose == 0 else None,
            )
            ax.plot(
                [position - 0.07, position + 0.07],
                [significant_scores.loc[mask, column].mean()] * 2,
                color="#111111",
                linewidth=1.4,
            )
    ax.axhline(0, color="#999999", linewidth=0.7)
    ax.set_xticks([0, 1, 2])
    ax.set_xlabel("Observed chr3 deletion dosage")
    ax.set_ylabel("PC30-adjusted methylation")
    ax.set_title("C  Two locus-wide corrected CpG associations", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(grid[2, 1])
    predefined = components[components["window"].eq("predefined_annotation")].copy()
    predefined = predefined.sort_values("beta")
    y = np.arange(len(predefined))
    ax.errorbar(predefined["beta"], y, xerr=1.96 * predefined["se"], fmt="o", color="#2b6f77", ecolor="#637381", capsize=2)
    ax.axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    ax.set_yticks(y, predefined["component"].str.replace("_", " "), fontsize=8)
    ax.set_xlabel("Mean methylation change per SV allele")
    ax.set_title("D  Predefined regional modules", loc="left", fontweight="bold")
    ax.text(0.98, 0.04, "No module: BH q < 0.05", transform=ax.transAxes, ha="right", fontsize=8)

    ax = fig.add_subplot(grid[3, :])
    for dose, group in expression_plot.groupby("dosage"):
        ax.scatter(group["cpg_residual"], group["expression_residual"], s=28, alpha=0.65, label=f"dosage {dose:g}")
    if len(expression_plot) > 1:
        slope, intercept = np.polyfit(expression_plot["cpg_residual"], expression_plot["expression_residual"], 1)
        x_line = np.linspace(expression_plot["cpg_residual"].min(), expression_plot["cpg_residual"].max(), 100)
        ax.plot(x_line, intercept + slope * x_line, color="#333333", linewidth=1.2)
    ax.set_xlabel("YPEL2-upstream CpG methylation residual")
    ax.set_ylabel("TUBD1 expression residual")
    ax.set_title("E  CpG-expression sensitivity after covariates and SV dosage", loc="left", fontweight="bold")
    ax.text(0.98, 0.94, f"PEER10 partial p = {expression_pvalue:.3g}", transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(output_dir / "tubd1_local_methylation_module.png", dpi=220, bbox_inches="tight")
    fig.savefig(output_dir / "tubd1_local_methylation_module.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    windows = {
        "gene_span_plus_100kb": (LOCUS_START - 100_000, LOCUS_END + 100_000),
        "gene_span_plus_500kb": (LOCUS_START - 500_000, LOCUS_END + 500_000),
        "gene_span_plus_1mb": (LOCUS_START - 1_000_000, LOCUS_END + 1_000_000),
    }
    maximum_start = min(value[0] for value in windows.values())
    maximum_end = max(value[1] for value in windows.values())

    covariates = pd.read_csv(args.covariates, sep="\t").set_index("sample_id")
    dosage_table = pd.read_csv(args.sv_dosages, sep="\t").set_index("sv_index")
    methylation = load_locus_matrix(args.methylation_matrix, maximum_start, maximum_end)
    samples = [sample for sample in covariates.index if sample in dosage_table.columns and sample in methylation.columns]
    methylation, values = prepare_methylation(
        methylation, samples, args.minimum_observed_fraction, args.minimum_sd
    )
    cov = covariates.loc[samples].to_numpy(dtype=float)
    x = dosage_table.loc[args.sv_index, samples].to_numpy(dtype=float)
    y_residual, covariate_rank = residualize(values.T, cov)
    x_residual, _ = residualize(x[:, None], cov)
    x_residual = x_residual[:, 0]

    screen = load_intervals(args.screen_ccres, maximum_start, maximum_end)
    gm12878 = load_intervals(args.gm12878_ccres, maximum_start, maximum_end)
    monocyte = load_intervals(args.monocyte_ccres, maximum_start, maximum_end)
    positions = methylation["start"].to_numpy(dtype=int)
    methylation["screen_ccre"] = annotate_interval_overlap(positions, screen)
    methylation["gm12878_active"] = annotate_interval_overlap(positions, gm12878)
    methylation["monocyte_active"] = annotate_interval_overlap(positions, monocyte)

    primary_start, primary_end = windows["gene_span_plus_500kb"]
    primary_mask = methylation["start"].between(primary_start, primary_end).to_numpy()
    primary = methylation.loc[primary_mask].copy().reset_index(drop=True)
    primary_residual = y_residual[:, primary_mask]
    beta, se, t_value, pvalue = association_statistics(
        x_residual, primary_residual, covariate_rank
    )
    primary["beta"] = beta[0]
    primary["se"] = se[0]
    primary["t"] = t_value[0]
    primary["pvalue"] = pvalue[0]
    primary["qvalue_bh"] = bh_qvalues(primary["pvalue"].to_numpy())

    rng = np.random.default_rng(args.seed)
    observed_max = float(np.max(np.abs(primary["t"])))
    null_maxima = permuted_max_statistics(
        x,
        primary_residual,
        cov,
        covariate_rank,
        args.cpg_permutations,
        rng,
    )
    primary["global_permutation_fwer"] = (
        1 + np.sum(null_maxima[:, None] >= np.abs(primary["t"].to_numpy())[None, :], axis=0)
    ) / (args.cpg_permutations + 1)
    primary_output_columns = list(MATRIX_META_COLUMNS) + [
        "screen_ccre",
        "gm12878_active",
        "monocyte_active",
        "beta",
        "se",
        "t",
        "pvalue",
        "qvalue_bh",
        "global_permutation_fwer",
    ]
    primary[primary_output_columns].to_csv(
        args.output_dir / "local_cpg_associations.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    component_rows: list[dict] = []
    loading_rows: list[pd.DataFrame] = []
    score_tables: dict[str, pd.DataFrame] = {}
    standardized_tables: dict[str, np.ndarray] = {}
    for window_name, (start, end) in windows.items():
        mask = methylation["start"].between(start, end).to_numpy()
        matrix = y_residual[:, mask]
        scale = np.std(matrix, axis=0, ddof=1)
        usable = scale > np.finfo(float).eps
        matrix = matrix[:, usable] / scale[usable]
        metadata = methylation.loc[mask].loc[usable].reset_index(drop=True)
        component_count = min(args.components, matrix.shape[0] - 1, matrix.shape[1])
        left, singular, right = randomized_svd(
            matrix, n_components=component_count, random_state=args.seed
        )
        scores, right = orient_components(left * singular[None, :], right)
        score_table = pd.DataFrame(
            scores,
            index=samples,
            columns=[f"PC{index + 1}" for index in range(component_count)],
        )
        score_tables[window_name] = score_table
        standardized_tables[window_name] = matrix
        score_residual, _ = residualize(scores, cov)
        comp_beta, comp_se, comp_t, comp_p = association_statistics(
            x_residual, score_residual, covariate_rank
        )
        variance = singular * singular / np.sum(singular * singular)
        for index in range(component_count):
            component_rows.append(
                {
                    "window": window_name,
                    "window_start": start,
                    "window_end": end,
                    "n_cpgs": matrix.shape[1],
                    "component": f"PC{index + 1}",
                    "variance_fraction_among_retained": variance[index],
                    "beta": comp_beta[0, index],
                    "se": comp_se[0, index],
                    "t": comp_t[0, index],
                    "pvalue": comp_p[0, index],
                }
            )
        loading = pd.DataFrame(right.T, columns=[f"PC{index + 1}" for index in range(component_count)])
        loading = pd.concat([metadata.reset_index(drop=True), loading], axis=1)
        loading.insert(0, "window", window_name)
        loading_rows.append(loading)

    predefined_scores, predefined_metadata = make_predefined_modules(
        methylation, y_residual
    )
    predefined_scores.index = samples
    if not predefined_scores.empty:
        predefined_residual, _ = residualize(predefined_scores.to_numpy(dtype=float), cov)
        predefined_beta, predefined_se, predefined_t, predefined_p = association_statistics(
            x_residual, predefined_residual, covariate_rank
        )
        for index, module in enumerate(predefined_scores.columns):
            metadata_row = predefined_metadata.set_index("module").loc[module]
            component_rows.append(
                {
                    "window": "predefined_annotation",
                    "window_start": metadata_row["start"],
                    "window_end": metadata_row["end"],
                    "n_cpgs": metadata_row["n_cpgs"],
                    "component": module,
                    "variance_fraction_among_retained": np.nan,
                    "beta": predefined_beta[0, index],
                    "se": predefined_se[0, index],
                    "t": predefined_t[0, index],
                    "pvalue": predefined_p[0, index],
                }
            )
        score_tables["predefined_annotation"] = predefined_scores

    components = pd.DataFrame(component_rows)
    components["qvalue_bh"] = bh_qvalues(components["pvalue"].to_numpy())
    observed_component_max = float(np.max(np.abs(components["t"])))
    component_null_maxima: list[float] = []
    score_residuals = [
        residualize(score_tables[name].to_numpy(dtype=float), cov)[0]
        for name in score_tables
    ]
    combined_scores = np.column_stack(score_residuals)
    component_null_maxima = permuted_max_statistics(
        x,
        combined_scores,
        cov,
        covariate_rank,
        args.component_permutations,
        rng,
        batch_size=1_000,
    )
    components["global_permutation_fwer"] = (
        1 + np.sum(
            component_null_maxima[:, None] >= np.abs(components["t"].to_numpy())[None, :],
            axis=0,
        )
    ) / (args.component_permutations + 1)
    components.to_csv(args.output_dir / "component_associations.tsv", sep="\t", index=False)
    loadings = pd.concat(loading_rows, ignore_index=True)
    loadings.to_csv(
        args.output_dir / "component_loadings.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    top = components.sort_values("pvalue").iloc[0]
    top_scores = score_tables[top["window"]][top["component"]].copy()
    top_scores.name = f"{top['window']}:{top['component']}"
    if top["window"] == "predefined_annotation":
        module_name = top["component"]
        if module_name == "VMP1_gene_body":
            top_mask = methylation["start"].between(*GENES["VMP1"][:2])
        elif module_name == "VMP1_promoter_2kb":
            top_mask = methylation["start"].between(GENES["VMP1"][0] - 2_000, GENES["VMP1"][0] + 2_000)
        elif module_name == "MIR21_plus_20kb":
            top_mask = methylation["start"].between(GENES["MIR21"][0] - 20_000, GENES["MIR21"][1] + 20_000)
        elif module_name == "TUBD1_gene_body":
            top_mask = methylation["start"].between(*GENES["TUBD1"][:2])
        elif module_name == "TUBD1_promoter_2kb":
            top_mask = methylation["start"].between(GENES["TUBD1"][1] - 2_000, GENES["TUBD1"][1] + 2_000)
        elif module_name == "SCREEN_cCRE_primary":
            top_mask = methylation["screen_ccre"] & methylation["start"].between(primary_start, primary_end)
        elif module_name == "GM12878_active_cCRE_primary":
            top_mask = methylation["gm12878_active"] & methylation["start"].between(primary_start, primary_end)
        else:
            top_mask = methylation["monocyte_active"] & methylation["start"].between(primary_start, primary_end)
        top_loading = methylation.loc[top_mask, list(MATRIX_META_COLUMNS) + ["screen_ccre", "gm12878_active", "monocyte_active"]].copy()
        top_loading["loading"] = 1.0 / max(int(top_mask.sum()), 1)
    else:
        top_loading = loadings[
            loadings["window"].eq(top["window"])
        ][list(MATRIX_META_COLUMNS) + ["screen_ccre", "gm12878_active", "monocyte_active", top["component"]]].copy()
        top_loading = top_loading.rename(columns={top["component"]: "loading"})
    top_loading["absolute_loading"] = top_loading["loading"].abs()
    top_loading = top_loading.sort_values("absolute_loading", ascending=False)
    top_loading.to_csv(args.output_dir / "top_component_cpgs.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"sample_id": samples, "sv_dosage": x, "module_score": top_scores.to_numpy(dtype=float)}
    ).to_csv(args.output_dir / "top_component_scores.tsv", sep="\t", index=False)

    expression = load_expression_rows(
        args.expression,
        {
            "ENSG00000108423": "TUBD1",
            "ENSG00000062716": "VMP1",
            "ENSG00000175155": "YPEL2",
        },
    )
    expression_covariates, genotype_pcs, peer_factors = load_expression_covariates(
        args.expression_covariates
    )
    expression_rows: list[dict] = []
    for gene, expression_values in expression.items():
        for peer_count in [0, 5, 10, 20, 30, 40, 60]:
            columns = ["sex"] + genotype_pcs[:5] + peer_factors[:peer_count]
            common = [
                sample for sample in samples
                if sample in expression_values.index and sample in expression_covariates.index
            ]
            outcome = expression_values.loc[common].to_numpy(dtype=float)
            local_score = top_scores.loc[common].to_numpy(dtype=float)
            local_dosage = pd.Series(x, index=samples).loc[common].to_numpy(dtype=float)
            expression_cov = expression_covariates.loc[common, columns].to_numpy(dtype=float)
            complete = np.isfinite(outcome) & np.isfinite(local_score) & np.isfinite(local_dosage)
            outcome = outcome[complete]
            local_score = local_score[complete]
            local_dosage = local_dosage[complete]
            expression_cov = expression_cov[complete]
            base = fit_ols(outcome, np.column_stack([expression_cov, local_dosage]))
            mediated = fit_ols(
                outcome, np.column_stack([expression_cov, local_dosage, local_score])
            )
            expression_rows.append(
                {
                    "gene": gene,
                    "module": top_scores.name,
                    "peer_factors": peer_count,
                    "n": len(outcome),
                    "df": mediated["df"],
                    "sv_beta": base["beta"][-1],
                    "sv_pvalue": base["pvalue"][-1],
                    "sv_beta_adjusted_module": mediated["beta"][-2],
                    "sv_pvalue_adjusted_module": mediated["pvalue"][-2],
                    "module_beta": mediated["beta"][-1],
                    "module_pvalue": mediated["pvalue"][-1],
                    "sv_attenuation_ratio": mediated["beta"][-2] / base["beta"][-1]
                    if base["beta"][-1] != 0
                    else np.nan,
                }
            )
    expression_results = pd.DataFrame(expression_rows)
    expression_results.to_csv(
        args.output_dir / "expression_mediation_sensitivity.tsv", sep="\t", index=False
    )

    significant_cpg_rows: list[dict] = []
    significant_cpgs = primary[primary["qvalue_bh"] < 0.01].sort_values("pvalue")
    primary_feature_to_index = {
        feature: index for index, feature in enumerate(primary["feature_id"])
    }
    for cpg in significant_cpgs.itertuples(index=False):
        cpg_score = pd.Series(
            primary_residual[:, primary_feature_to_index[cpg.feature_id]], index=samples
        )
        for gene, expression_values in expression.items():
            for peer_count in [0, 5, 10, 20, 30, 40, 60]:
                columns = ["sex"] + genotype_pcs[:5] + peer_factors[:peer_count]
                common = [
                    sample
                    for sample in samples
                    if sample in expression_values.index and sample in expression_covariates.index
                ]
                outcome = expression_values.loc[common].to_numpy(dtype=float)
                local_cpg = cpg_score.loc[common].to_numpy(dtype=float)
                local_dosage = pd.Series(x, index=samples).loc[common].to_numpy(dtype=float)
                expression_cov = expression_covariates.loc[common, columns].to_numpy(dtype=float)
                complete = np.isfinite(outcome) & np.isfinite(local_cpg) & np.isfinite(local_dosage)
                outcome = outcome[complete]
                local_cpg = local_cpg[complete]
                local_dosage = local_dosage[complete]
                expression_cov = expression_cov[complete]
                base = fit_ols(outcome, np.column_stack([expression_cov, local_dosage]))
                mediated = fit_ols(
                    outcome, np.column_stack([expression_cov, local_dosage, local_cpg])
                )
                significant_cpg_rows.append(
                    {
                        "cpg": cpg.feature_id,
                        "cpg_beta_sv": cpg.beta,
                        "cpg_qvalue_bh": cpg.qvalue_bh,
                        "cpg_permutation_fwer": cpg.global_permutation_fwer,
                        "gene": gene,
                        "peer_factors": peer_count,
                        "n": len(outcome),
                        "sv_beta": base["beta"][-1],
                        "sv_pvalue": base["pvalue"][-1],
                        "sv_beta_adjusted_cpg": mediated["beta"][-2],
                        "sv_pvalue_adjusted_cpg": mediated["pvalue"][-2],
                        "cpg_beta_expression": mediated["beta"][-1],
                        "cpg_pvalue_expression": mediated["pvalue"][-1],
                    }
                )
    significant_cpg_results = pd.DataFrame(significant_cpg_rows)
    significant_cpg_results.to_csv(
        args.output_dir / "significant_cpg_expression_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    dosage_series = pd.Series(x, index=samples, name="SV dosage")
    significant_scores = pd.DataFrame(
        {
            row.feature_id: primary_residual[
                :, primary_feature_to_index[row.feature_id]
            ]
            for row in significant_cpgs.itertuples(index=False)
        },
        index=samples,
    )
    expression_cpg = significant_cpgs.iloc[min(1, len(significant_cpgs) - 1)]["feature_id"]
    expression_values = expression["TUBD1"]
    expression_columns = ["sex"] + genotype_pcs[:5] + peer_factors[:10]
    expression_common = [
        sample
        for sample in samples
        if sample in expression_values.index and sample in expression_covariates.index
    ]
    expression_outcome = expression_values.loc[expression_common].to_numpy(dtype=float)
    expression_methylation = significant_scores.loc[
        expression_common, expression_cpg
    ].to_numpy(dtype=float)
    expression_dosage = dosage_series.loc[expression_common].to_numpy(dtype=float)
    expression_cov = expression_covariates.loc[
        expression_common, expression_columns
    ].to_numpy(dtype=float)
    expression_design = np.column_stack(
        [np.ones(len(expression_common)), expression_cov, expression_dosage]
    )
    expression_residual = expression_outcome - expression_design @ np.linalg.lstsq(
        expression_design, expression_outcome, rcond=None
    )[0]
    methylation_residual = expression_methylation - expression_design @ np.linalg.lstsq(
        expression_design, expression_methylation, rcond=None
    )[0]
    expression_plot = pd.DataFrame(
        {
            "cpg_residual": methylation_residual,
            "expression_residual": expression_residual,
            "dosage": expression_dosage,
        },
        index=expression_common,
    )
    expression_pvalue = float(
        significant_cpg_results[
            significant_cpg_results["cpg"].eq(expression_cpg)
            & significant_cpg_results["gene"].eq("TUBD1")
            & significant_cpg_results["peer_factors"].eq(10)
        ]["cpg_pvalue_expression"].iloc[0]
    )
    plot_results(
        primary,
        components,
        significant_scores,
        dosage_series,
        expression_plot,
        screen,
        gm12878,
        expression_pvalue,
        args.output_dir,
    )

    summary = {
        "locus_chm13": f"chr17:{LOCUS_START}-{LOCUS_END}",
        "primary_window": f"chr17:{primary_start}-{primary_end}",
        "n_samples": len(samples),
        "n_primary_cpgs": len(primary),
        "n_primary_bh_q_lt_0_01": int((primary["qvalue_bh"] < 0.01).sum()),
        "n_primary_permutation_fwer_lt_0_05": int((primary["global_permutation_fwer"] < 0.05).sum()),
        "n_modules_bh_q_lt_0_05": int((components["qvalue_bh"] < 0.05).sum()),
        "n_modules_permutation_fwer_lt_0_05": int((components["global_permutation_fwer"] < 0.05).sum()),
        "top_cpg": primary.sort_values("pvalue").iloc[0]["feature_id"],
        "top_cpg_pvalue": float(primary["pvalue"].min()),
        "top_component": top_scores.name,
        "top_component_beta": float(top["beta"]),
        "top_component_pvalue": float(top["pvalue"]),
        "top_component_qvalue_bh": float(top["qvalue_bh"]),
        "top_component_global_permutation_fwer": float(top["global_permutation_fwer"]),
        "observed_cpg_max_abs_t": observed_max,
        "observed_component_max_abs_t": observed_component_max,
        "component_permutations": args.component_permutations,
        "cpg_permutations": args.cpg_permutations,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
