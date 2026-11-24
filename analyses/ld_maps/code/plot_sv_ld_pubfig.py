"""Publication-grade variant of plot_sv_ld.py."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import os as _censv_os
CENSV_COLLAB_ROOT = _censv_os.environ["CENSV_COLLAB_ROOT"]

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch, Rectangle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chr-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sv-category", default="original_decomposed")
    parser.add_argument("--min-maf", type=float, default=0.0)
    parser.add_argument("--max-variants-for-plot", type=int, default=250)
    parser.add_argument("--top-pairs", type=int, default=500)
    parser.add_argument("--snp-vcf", type=Path)
    parser.add_argument("--snp-region")
    parser.add_argument("--snp-min-maf", type=float, default=0.01)
    parser.add_argument("--cross-max-sv-plot", type=int, default=80)
    parser.add_argument("--cross-max-snp-plot", type=int, default=160)
    parser.add_argument("--upstream-extension-bp", type=int, default=100_000)
    parser.add_argument(
        "--allow-sv-only-fallback", action="store_true",
        help="If --snp-region returns zero biallelic SNPs at the requested MAF, "
        "fall back to an SV-only triangle for this chrom instead of erroring out.",
    )
    parser.add_argument(
        "--active-bed",
        type=Path,
        default=Path(CENSV_COLLAB_ROOT + "/resources/chm13v2.0_censat_v2.1.active.bed"),
    )
    parser.add_argument("--kinetochore-region")
    return parser.parse_args()


def run_command(cmd: list[str]) -> str:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return result.stdout


def infer_prefix(chr_dir: Path) -> str:
    npz_files = sorted(chr_dir.glob("*_sv_genotypes.npz"))
    if len(npz_files) != 1:
        raise FileNotFoundError(f"Expected exactly one *_sv_genotypes.npz in {chr_dir}")
    return npz_files[0].name.removesuffix("_sv_genotypes.npz")


def load_sv_inputs(chr_dir: Path) -> tuple[str, np.ndarray, pd.DataFrame, list[str]]:
    prefix = infer_prefix(chr_dir)
    genotypes = np.load(chr_dir / f"{prefix}_sv_genotypes.npz")["genotypes"]
    variants = pd.read_csv(chr_dir / f"{prefix}_sv_variants.tsv", sep="\t")
    samples = [
        line.strip()
        for line in (chr_dir / f"{prefix}_sv_samples.txt").read_text().splitlines()
        if line.strip()
    ]
    if genotypes.shape[0] != len(variants):
        raise ValueError("Variant table and genotype matrix row counts do not match")
    if genotypes.shape[1] != len(samples):
        raise ValueError("Sample list and genotype matrix column counts do not match")
    return prefix, genotypes, variants, samples


def normalize_variant_frame(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    frame = frame.copy()
    if "chrom" not in frame.columns:
        frame["chrom"] = ""
    if "end" not in frame.columns:
        frame["end"] = frame["pos"]
    if "svtype" not in frame.columns:
        frame["svtype"] = source
    if "af" not in frame.columns:
        frame["af"] = np.nan
    if "maf" not in frame.columns:
        frame["maf"] = np.nan
    frame["source"] = source
    return frame


def finalize_variant_matrix(genotypes: np.ndarray, variants: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    if genotypes.shape[0] != len(variants):
        raise ValueError("Variant frame and genotype matrix row counts do not match")
    row_means = np.nanmean(genotypes, axis=1)
    keep = ~np.isnan(row_means)
    genotypes = genotypes[keep, :]
    variants = variants.loc[keep].copy()
    row_means = np.nanmean(genotypes, axis=1)
    missing_mask = np.isnan(genotypes)
    if np.any(missing_mask):
        genotypes = genotypes.copy()
        genotypes[missing_mask] = np.take(row_means, np.where(missing_mask)[0])
    variable_mask = genotypes.var(axis=1) > 0
    genotypes = genotypes[variable_mask, :]
    variants = variants.loc[variable_mask].copy()
    order = np.lexsort(
        (
            variants["id"].astype(str).to_numpy(),
            variants["end"].to_numpy(),
            variants["pos"].to_numpy(),
        )
    )
    return genotypes[order, :], variants.iloc[order].reset_index(drop=True)


def filter_sv_variants(
    genotypes: np.ndarray,
    variants: pd.DataFrame,
    sv_category: str,
    min_maf: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    mask = variants["sv_category"].eq(sv_category)
    if "maf" in variants.columns:
        mask &= variants["maf"].fillna(0.0).ge(min_maf)
    filtered_variants = normalize_variant_frame(variants.loc[mask].copy(), source="centromere_sv")
    filtered_genotypes = genotypes[mask.to_numpy(), :]
    if filtered_genotypes.size == 0:
        raise ValueError(f"No SVs remain after filtering sv_category={sv_category!r} and min_maf={min_maf}")
    return finalize_variant_matrix(filtered_genotypes, filtered_variants)


def gt_to_dosage(gt: str) -> float:
    if gt in {".", "./.", ".|."}:
        return np.nan
    gt_core = gt.split(":", 1)[0]
    sep = "/" if "/" in gt_core else "|" if "|" in gt_core else None
    if sep is None:
        return np.nan
    alleles = gt_core.split(sep)
    if len(alleles) != 2 or "." in alleles:
        return np.nan
    try:
        vals = [int(a) for a in alleles]
    except ValueError:
        return np.nan
    if max(vals) > 1:
        return np.nan
    return float(sum(vals))


def query_vcf_samples(vcf_path: Path) -> list[str]:
    output = run_command(["bcftools", "query", "-l", str(vcf_path)])
    return [line.strip() for line in output.splitlines() if line.strip()]


def load_snp_panel_from_vcf(
    vcf_path: Path,
    region: str,
    target_samples: list[str],
    min_maf: float,
) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    vcf_samples = query_vcf_samples(vcf_path)
    vcf_set = set(vcf_samples)
    keep_samples = [s for s in target_samples if s in vcf_set]
    if not keep_samples:
        raise ValueError("No overlapping samples between the SV matrix and the SNP VCF")
    with tempfile.NamedTemporaryFile("w", suffix=".samples.txt", delete=False) as handle:
        handle.write("\n".join(keep_samples) + "\n")
        sample_file = Path(handle.name)
    try:
        output = run_command([
            "bcftools", "query", "-r", region, "-S", str(sample_file),
            "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n", str(vcf_path),
        ])
    finally:
        sample_file.unlink(missing_ok=True)

    variant_rows: list[dict[str, object]] = []
    genotype_rows: list[list[float]] = []
    for line in output.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5 + len(keep_samples):
            continue
        chrom, pos, variant_id, ref, alt = fields[:5]
        gt_fields = fields[5:]
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        dosage = np.array([gt_to_dosage(gt) for gt in gt_fields], dtype=float)
        if np.all(np.isnan(dosage)):
            continue
        af = float(np.nanmean(dosage) / 2.0)
        maf = min(af, 1.0 - af)
        if maf < min_maf:
            continue
        variant_rows.append({
            "chrom": chrom, "pos": int(pos), "end": int(pos),
            "id": variant_id if variant_id != "." else f"{chrom}:{pos}:{ref}:{alt}",
            "svtype": "SNP", "ref": ref, "alt": alt,
            "af": af, "maf": maf, "source": "peri_snp",
        })
        genotype_rows.append(dosage.tolist())
    if not genotype_rows:
        raise ValueError(f"No biallelic SNPs remain in {region} after filtering min_maf={min_maf}")
    genotypes = np.asarray(genotype_rows, dtype=float)
    variants = pd.DataFrame(variant_rows)
    genotypes, variants = finalize_variant_matrix(genotypes, variants)
    return genotypes, variants, keep_samples


def standardize_rows(genotypes: np.ndarray) -> np.ndarray:
    centered = genotypes - genotypes.mean(axis=1, keepdims=True)
    scale = genotypes.std(axis=1, keepdims=True)
    scale[scale == 0] = 1.0
    return centered / scale


def compute_r(genotypes: np.ndarray) -> np.ndarray:
    z = standardize_rows(genotypes.astype(float))
    r = (z @ z.T) / z.shape[1]
    return np.clip(r, -1.0, 1.0)


def compute_cross_r(left_genotypes: np.ndarray, right_genotypes: np.ndarray) -> np.ndarray:
    left = standardize_rows(left_genotypes.astype(float))
    right = standardize_rows(right_genotypes.astype(float))
    n_samples = left.shape[1]
    r = (left @ right.T) / n_samples
    return np.clip(r, -1.0, 1.0)


def select_plot_subset(variants: pd.DataFrame, max_variants: int) -> tuple[pd.DataFrame, np.ndarray]:
    if len(variants) <= max_variants:
        return variants.copy(), np.arange(len(variants))
    subset = (
        variants.sort_values(["maf", "af", "pos"], ascending=[False, False, True])
        .head(max_variants)
        .sort_values(["pos", "end", "id"])
    )
    return subset.reset_index(drop=True), subset.index.to_numpy()


def parse_region_string(region: str) -> tuple[str, int, int]:
    chrom, coords = region.split(":")
    start_s, end_s = coords.replace(",", "").split("-")
    return chrom, int(start_s), int(end_s)


def load_active_span(active_bed: Path, chrom: str) -> dict[str, object] | None:
    if not active_bed.exists():
        return None
    best = None
    with active_bed.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] != chrom:
                continue
            start, end, label = int(fields[1]), int(fields[2]), fields[3]
            length = end - start
            if best is None or length > best["length"]:
                best = {"label": f"Active HOR ({label})", "start": start, "end": end, "length": length}
    if best is None:
        return None
    best.pop("length")
    return best


def get_centromere_annotations(active_bed: Path, chrom: str, kinetochore: str | None) -> list[dict[str, object]]:
    spans: list[dict[str, object]] = []
    active = load_active_span(active_bed, chrom)
    if active:
        spans.append(active)
    if kinetochore:
        k_chrom, k_start, k_end = parse_region_string(kinetochore)
        if k_chrom == chrom:
            spans.append({"label": "Kinetochore Region", "start": k_start, "end": k_end})
    return spans


def pos_to_column(variants: pd.DataFrame, target_pos: int) -> float | None:
    n = len(variants)
    if n == 0:
        return None
    pos_arr = variants["pos"].to_numpy()
    if target_pos <= pos_arr[0]:
        return 0.0
    if target_pos >= pos_arr[-1]:
        return float(n)
    right_idx = int(np.searchsorted(pos_arr, target_pos))
    left_idx = right_idx - 1
    left_pos = pos_arr[left_idx]
    right_pos = pos_arr[right_idx]
    if right_pos == left_pos:
        return float(right_idx)
    frac = (target_pos - left_pos) / (right_pos - left_pos)
    return left_idx + frac + 0.5


def pos_range_to_columns(variants: pd.DataFrame, start: int, end: int) -> tuple[float, float] | None:
    left = pos_to_column(variants, start)
    right = pos_to_column(variants, end)
    if left is None or right is None:
        return None
    return left, right


def build_triangle_polygons(r: np.ndarray) -> tuple[list[list[tuple[float, float]]], np.ndarray]:
    n = r.shape[0]
    polygons: list[list[tuple[float, float]]] = []
    values: list[float] = []
    for i in range(n):
        for j in range(i, n):
            corners = []
            for x, y in ((j, i), (j + 1, i), (j + 1, i + 1), (j, i + 1)):
                corners.append(((x + y) / 2.0, (x - y) / 2.0))
            polygons.append(corners)
            values.append(r[i, j])
    return polygons, np.asarray(values, dtype=float)


def genomic_ticks(variants: pd.DataFrame, max_ticks: int = 7) -> tuple[np.ndarray, list[str]]:
    n = len(variants)
    if n == 0:
        return np.array([]), []
    ticks = np.unique(np.linspace(0, n - 1, num=min(max_ticks, n), dtype=int))
    labels = [f"{variants.iloc[i]['pos'] / 1e6:.2f}" for i in ticks]
    return ticks, labels


def format_window_label(bp: int) -> str:
    if bp >= 1_000_000 and bp % 1_000_000 == 0:
        return f"{bp // 1_000_000} Mb"
    if bp >= 1_000 and bp % 1_000 == 0:
        return f"{bp // 1_000} kb"
    return f"{bp} bp"


def plot_triangle_pubfig(
    r: np.ndarray,
    variants: pd.DataFrame,
    title: str,
    output_jpg: Path,
    centromere_span: dict[str, object] | None,
    upstream_extension_pos: int | None,
    upstream_extension_bp: int,
    boundary_idx: int | None,
    fallback_note: str | None = None,
) -> None:
    n = len(variants)
    fig_w = max(18, min(36, n * 0.135))
    fig_h = max(10, min(20, n * 0.085))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    polygons, values = build_triangle_polygons(r)
    collection = PolyCollection(
        polygons, array=values, cmap="coolwarm", edgecolors="none", clim=(-1.0, 1.0),
    )
    ax.add_collection(collection)

    triangle_height = n / 2.0
    extra_top = max(4.0, triangle_height * 0.36)
    extra_bot = max(4.0, triangle_height * 0.28)
    xmargin_left = max(3.0, 0.06 * n) if boundary_idx is not None else max(1.0, 0.015 * n)
    xmargin_right = max(1.0, 0.015 * n)
    ax.set_xlim(-xmargin_left, n + xmargin_right)
    ax.set_ylim(-extra_bot, triangle_height + extra_top)
    ax.set_aspect("equal")

    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"

    legend_handles_top: list[Patch] = []
    legend_labels_top: list[str] = []
    window_label = format_window_label(upstream_extension_bp)
    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_height + extra_top * 0.16
        strip_h = max(1.4, extra_top * 0.10)
        ax.add_patch(Rectangle(
            (0, strip_y), boundary_idx, strip_h,
            facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=1.4, zorder=2,
        ))
        ax.add_patch(Rectangle(
            (boundary_idx, strip_y), n - boundary_idx, strip_h,
            facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=1.4, zorder=2,
        ))
        legend_handles_top.append(Patch(
            facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=1.4,
        ))
        legend_labels_top.append(
            f"peri-centromeric SNPs  ·  {window_label} upstream  ·  n = {boundary_idx}"
        )
        cen_mb = None
        if centromere_span is not None:
            cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6
        legend_handles_top.append(Patch(
            facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=1.4,
        ))
        legend_labels_top.append(
            f"Active HOR (centromere)  ·  {cen_mb:.2f} Mb  ·  n = {n - boundary_idx} SVs"
            if cen_mb is not None else
            f"Active HOR (centromeric SVs)  ·  n = {n - boundary_idx}"
        )

    if boundary_idx is not None and 0 < boundary_idx < n:
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=2.6, zorder=4)
        ax.annotate(
            "Active HOR start (centromere edge)",
            xy=(boundary_idx, triangle_height + extra_top * 0.03),
            xytext=(boundary_idx, triangle_height + extra_top * 0.55),
            ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
            arrowprops=dict(arrowstyle="-", lw=1.8, color=boundary_colour),
        )

    if boundary_idx is not None and 0 < boundary_idx < n:
        wedge_text_colour = "#2a2a2a"
        wedge_edge_colour = "#404040"
        bx = float(boundary_idx)
        fn = float(n)
        snp_wedge = [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)]
        ax.plot(
            [p[0] for p in snp_wedge], [p[1] for p in snp_wedge],
            color=wedge_edge_colour, lw=2.6, linestyle=(0, (5, 3)), zorder=5,
        )
        ax.text(
            bx * 0.50, bx * 0.20,
            "SNP × SNP",
            ha="center", va="center", fontsize=40, color=wedge_text_colour,
            style="italic", weight="bold", alpha=0.95, zorder=6,
        )
        sv_wedge = [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)]
        ax.plot(
            [p[0] for p in sv_wedge], [p[1] for p in sv_wedge],
            color=wedge_edge_colour, lw=2.6, linestyle=(0, (5, 3)), zorder=5,
        )
        ax.text(
            (bx + fn) / 2.0, (fn - bx) * 0.20,
            "SV × SV",
            ha="center", va="center", fontsize=40, color=wedge_text_colour,
            style="italic", weight="bold", alpha=0.95, zorder=6,
        )

    if boundary_idx is not None and 0 < boundary_idx < n:
        cross_colour = "#6c3baa"
        bx = float(boundary_idx)
        fn = float(n)
        rhombus = [
            (bx / 2.0, bx / 2.0),
            (fn / 2.0, fn / 2.0),
            ((bx + fn) / 2.0, (fn - bx) / 2.0),
            (bx, 0.0),
            (bx / 2.0, bx / 2.0),
        ]
        xs = [p[0] for p in rhombus]
        ys = [p[1] for p in rhombus]
        ax.plot(xs, ys, color=cross_colour, lw=1.6, linestyle=(0, (5, 3)), zorder=5)
        cross_label_x = (bx + fn) / 2.0 - bx / 4.0
        cross_label_y = (fn - bx) / 2.0 + bx / 4.0
        ax.annotate(
            "SV × SNP cross-LD",
            xy=(cross_label_x, cross_label_y),
            xytext=(cross_label_x + max(6, n * 0.04), cross_label_y + max(4, triangle_height * 0.10)),
            ha="left", va="bottom", fontsize=42, color=cross_colour, weight="bold",
            arrowprops=dict(arrowstyle="-", lw=1.8, color=cross_colour),
        )

    upstream_col = None
    if boundary_idx is not None and upstream_extension_pos is not None:
        upstream_col = pos_to_column(variants, upstream_extension_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=2.6, zorder=4)
        upstream_lbl = f"−{format_window_label(upstream_extension_bp)} upstream"
        ax.annotate(
            upstream_lbl,
            xy=(upstream_col, -extra_bot * 0.05),
            xytext=(upstream_col, -extra_bot * 0.92),
            ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
            arrowprops=dict(arrowstyle="-", lw=1.8, color=upstream_colour),
        )

    if boundary_idx is None and centromere_span is not None:
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6
        strip_y = triangle_height + extra_top * 0.18
        strip_h = max(1.0, extra_top * 0.10)
        ax.add_patch(Rectangle(
            (0, strip_y), n, strip_h,
            facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=0.8, zorder=2,
        ))
        ax.text(
            n / 2.0, strip_y + strip_h / 2.0,
            f"Active HOR (centromere, {cen_mb:.2f} Mb) — n={n} SVs",
            ha="center", va="center", fontsize=16, color="#1a1a1a", weight="bold",
        )

    if fallback_note:
        ax.text(
            0.5, 1.005,
            fallback_note,
            transform=ax.transAxes,
            ha="center", va="bottom", fontsize=14,
            color="#555555", style="italic",
        )

    tick_positions: list[float] = []
    tick_labels: list[str] = []
    if upstream_col is not None and upstream_extension_pos is not None:
        tick_positions.append(upstream_col)
        tick_labels.append(f"{upstream_extension_pos / 1e6:.3f}")
    if boundary_idx is not None and 0 < boundary_idx < n and centromere_span is not None:
        tick_positions.append(float(boundary_idx))
        tick_labels.append(f"{int(centromere_span['start']) / 1e6:.3f}")
    elif boundary_idx is None and centromere_span is not None and n > 0:
        first_pos = int(variants.iloc[0]["pos"])
        tick_positions.append(0.0)
        tick_labels.append(f"{first_pos / 1e6:.3f}")
    if n > 0:
        last_pos = int(variants.iloc[-1]["pos"])
        tick_positions.append(float(n))
        tick_labels.append(f"{last_pos / 1e6:.3f}")
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=17)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=18, weight="semibold")
    ax.set_yticks([])
    ax.set_ylabel("")
    title_pad = 34 if fallback_note else 14
    ax.set_title(title, fontsize=19, pad=title_pad, weight="bold")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#4a4a4a")
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="x", which="major", length=5, width=1.0, color="#4a4a4a", pad=5)

    cbar = fig.colorbar(collection, ax=ax, fraction=0.024, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=17, weight="semibold")
    cbar.ax.tick_params(labelsize=15)

    fig.tight_layout()
    fig.savefig(output_jpg, dpi=300, format="jpg", bbox_inches="tight", pad_inches=0.25)
    fig.savefig(output_jpg.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prefix, sv_genotypes, sv_variants_raw, sv_samples = load_sv_inputs(args.chr_dir)
    sv_genotypes, sv_variants = filter_sv_variants(
        genotypes=sv_genotypes, variants=sv_variants_raw,
        sv_category=args.sv_category, min_maf=args.min_maf,
    )
    chrom = str(sv_variants.iloc[0]["chrom"])
    centromere_spans = get_centromere_annotations(args.active_bed, chrom, args.kinetochore_region)
    active_span = centromere_spans[0] if centromere_spans else None

    snp_panel_loaded = False
    snp_genotypes = None
    snp_variants = None
    overlap_samples: list[str] = []
    fallback_note = None
    upstream_pos = None
    if active_span is not None:
        upstream_pos = int(active_span["start"]) - args.upstream_extension_bp

    if args.snp_vcf is not None and args.snp_region:
        try:
            snp_genotypes, snp_variants, overlap_samples = load_snp_panel_from_vcf(
                vcf_path=args.snp_vcf, region=args.snp_region,
                target_samples=sv_samples, min_maf=args.snp_min_maf,
            )
            snp_panel_loaded = True
        except ValueError as exc:
            if not args.allow_sv_only_fallback:
                raise
            fallback_note = (
                f"1KG SNP panel empty in {args.snp_region}; "
                "falling back to SV-only triangle."
            )

    snp_label = args.snp_region.replace(":", "_").replace("-", "_") if args.snp_region else "svonly"

    if snp_panel_loaded:
        sv_idx_overlap = np.array([sv_samples.index(s) for s in overlap_samples], dtype=int)
        sv_overlap_genotypes = sv_genotypes[:, sv_idx_overlap]
        plot_snp, snp_idx = select_plot_subset(snp_variants, args.cross_max_snp_plot)
        plot_sv, sv_plot_idx = select_plot_subset(sv_variants, args.cross_max_sv_plot)
        combined_variants = pd.concat([plot_snp, plot_sv], ignore_index=True)
        combined_genotypes = np.vstack(
            [snp_genotypes[snp_idx, :], sv_overlap_genotypes[sv_plot_idx, :]]
        )
        combined_r = compute_r(combined_genotypes)
        boundary_idx: int | None = len(plot_snp)
        out_prefix = args.output_dir / f"{prefix}_{args.sv_category}_plus_peri_snp_{snp_label}_pubfig"
        title = f"{prefix} peri-SNP + centromeric SV LD (Pearson r)"
    else:
        plot_sv, sv_plot_idx = select_plot_subset(sv_variants, args.cross_max_sv_plot)
        combined_variants = plot_sv.reset_index(drop=True)
        combined_r = compute_r(sv_genotypes[sv_plot_idx, :])
        boundary_idx = None
        plot_snp = combined_variants.iloc[:0].copy()
        out_prefix = args.output_dir / f"{prefix}_{args.sv_category}_svonly_pubfig"
        title = f"{prefix} centromeric SV LD (Pearson r) — SV-only (no 1KG SNPs upstream)"

    pd.DataFrame(combined_r, index=combined_variants["id"], columns=combined_variants["id"]).to_csv(
        f"{out_prefix}_r.tsv.gz", sep="\t", compression="gzip"
    )
    combined_variants.to_csv(f"{out_prefix}_variants.tsv", sep="\t", index=False)

    plot_triangle_pubfig(
        combined_r, combined_variants,
        title=title,
        output_jpg=Path(f"{out_prefix}_r_triangle.jpg"),
        centromere_span=active_span,
        upstream_extension_pos=upstream_pos,
        upstream_extension_bp=args.upstream_extension_bp,
        boundary_idx=boundary_idx,
        fallback_note=fallback_note,
    )

    summary = {
        "prefix": prefix,
        "chrom": chrom,
        "sv_category": args.sv_category,
        "n_samples_sv_panel": len(sv_samples),
        "n_samples_overlap": len(overlap_samples),
        "snp_vcf": str(args.snp_vcf) if args.snp_vcf else None,
        "snp_region": args.snp_region,
        "upstream_extension_bp": args.upstream_extension_bp,
        "upstream_extension_pos": upstream_pos,
        "active_span": active_span,
        "mode": "combined" if snp_panel_loaded else "sv_only_fallback",
        "fallback_note": fallback_note,
        "combined_panel": {
            "n_snp_variants_total": int(len(snp_variants)) if snp_variants is not None else 0,
            "n_sv_variants_total": int(len(sv_variants)),
            "n_snp_variants_plotted": int(len(plot_snp)),
            "n_sv_variants_plotted": int(len(combined_variants) - len(plot_snp)),
        },
        "output_jpg": str(Path(f"{out_prefix}_r_triangle.jpg")),
    }
    Path(f"{out_prefix}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
