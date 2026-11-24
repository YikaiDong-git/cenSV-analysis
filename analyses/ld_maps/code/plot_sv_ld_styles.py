"""Ten distinct figure-design variants of the chr peri-SNP + centromeric SV LD plot."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch, Polygon, Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from plot_sv_ld_pubfig import (  # noqa: E402
    compute_cross_r,
    compute_r,
    filter_sv_variants,
    format_window_label,
    get_centromere_annotations,
    load_snp_panel_from_vcf,
    load_sv_inputs,
    pos_to_column,
    run_command,
    select_plot_subset,
)


def find_first_snp_pos(vcf_path: Path, region: str) -> int | None:
    import subprocess as _sp
    try:
        proc = _sp.run(
            [
                "bash", "-c",
                f"set -o pipefail; bcftools query -r {region} -f '%POS\\n' "
                f"{vcf_path} 2>/dev/null | head -1",
            ],
            check=False, capture_output=True, text=True,
        )
    except Exception:
        return None
    line = (proc.stdout or "").strip()
    if not line:
        return None
    try:
        return int(line)
    except ValueError:
        return None




def build_rotated_polygons(r: np.ndarray) -> tuple[list[list[tuple[float, float]]], np.ndarray]:
    n = r.shape[0]
    polygons: list[list[tuple[float, float]]] = []
    values: list[float] = []
    for i in range(n):
        for j in range(i, n):
            corners = [
                ((j + i) / 2.0, (j - i) / 2.0),
                ((j + 1 + i) / 2.0, (j + 1 - i) / 2.0),
                ((j + 1 + i + 1) / 2.0, (j + 1 - i - 1) / 2.0),
                ((j + i + 1) / 2.0, (j - i - 1) / 2.0),
            ]
            polygons.append(corners)
            values.append(r[i, j])
    return polygons, np.asarray(values, dtype=float)


def minimal_ticks(
    ax: plt.Axes,
    variants: pd.DataFrame,
    boundary_idx: int | None,
    centromere_start: int | None,
    upstream_pos: int | None,
    upstream_col: float | None,
    n: int,
    boundary_idx_right: int | None = None,
    centromere_end: int | None = None,
    downstream_pos: int | None = None,
    downstream_col: float | None = None,
) -> None:
    positions: list[float] = []
    labels: list[str] = []
    if upstream_col is not None and upstream_pos is not None:
        positions.append(upstream_col)
        labels.append(f"{upstream_pos / 1e6:.3f}")
    if boundary_idx is not None and centromere_start is not None and 0 < boundary_idx < n:
        positions.append(float(boundary_idx))
        labels.append(f"{centromere_start / 1e6:.3f}")
    if (boundary_idx_right is not None and centromere_end is not None
            and 0 < boundary_idx_right < n):
        positions.append(float(boundary_idx_right))
        labels.append(f"{centromere_end / 1e6:.3f}")
    if downstream_col is not None and downstream_pos is not None:
        positions.append(downstream_col)
        labels.append(f"{downstream_pos / 1e6:.3f}")
    if n > 0 and boundary_idx_right is None:
        last_pos = int(variants.iloc[-1]["pos"])
        positions.append(float(n))
        labels.append(f"{last_pos / 1e6:.3f}")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=42)




def render_v01(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]
    fallback_note = payload["fallback_note"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(4.0, triangle_h * 0.30)
    xmargin_left = max(4.0, 0.07 * n) if boundary_idx is not None else max(1.0, 0.02 * n)
    xmargin_right = max(1.0, 0.02 * n)
    ax.set_xlim(-xmargin_left, n + xmargin_right)
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    wedge_edge = "#404040"
    wedge_txt = "#2a2a2a"

    legend_handles: list[Patch] = []
    legend_labels: list[str] = []
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0, zorder=2))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0, zorder=2))
        legend_handles.append(Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        legend_labels.append(payload["snp_legend_label"])
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        if cen_mb:
            legend_handles.append(Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
            legend_labels.append(payload["sv_legend_label"])

    if boundary_idx is not None and 0 < boundary_idx < n:
        bx, fn = float(boundary_idx), float(n)
        for wedge in [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
            [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)],
        ]:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=3.0, linestyle=(0, (6, 4)), zorder=5)
        if n < 280:
            ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                    fontsize=32, color="#555", style="italic", weight="medium",
                    alpha=0.85, zorder=6)
            ax.text((bx + fn) / 2.0, (fn - bx) * 0.22, "SV × SV", ha="center", va="center",
                    fontsize=32, color="#555", style="italic", weight="medium",
                    alpha=0.85, zorder=6)
        rhombus = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                   ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        ax.plot([p[0] for p in rhombus], [p[1] for p in rhombus],
                color=cross_colour, lw=3.0, linestyle=(0, (6, 3)), zorder=5)
        clx, cly = (bx + fn) / 2.0 - bx / 4.0, (fn - bx) / 2.0 + bx / 4.0
        ax.annotate("SV × SNP cross-LD", xy=(clx, cly),
                    xytext=(clx + max(6, n * 0.04), cly + max(4, triangle_h * 0.10)),
                    ha="left", va="bottom", fontsize=68, color=cross_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=cross_colour))
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=len(legend_handles), frameon=True, framealpha=0.95,
                        edgecolor="0.4", fontsize=36,
                        handlelength=2.5, handleheight=1.8, borderpad=0.8,
                        columnspacing=3.0)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#3a3a3a")
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", which="major", length=10, width=2.0, color="#3a3a3a", pad=10)

    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold")

    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    if fallback_note:
        ax.text(0.5, 1.005, fallback_note, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=36, color="#555", style="italic")

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v02(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(11, min(20, n * 0.085))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(3.0, triangle_h * 0.22)
    extra_bot = max(3.0, triangle_h * 0.22)
    xmargin_left = max(4.0, 0.07 * n) if boundary_idx is not None else max(1.0, 0.02 * n)
    xmargin_right = max(1.0, 0.02 * n)
    ax.set_xlim(-xmargin_left, n + xmargin_right)
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    boundary_colour = "#1a1a1a"
    upstream_colour = "#2a9d8f"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))

    if boundary_idx is not None and 0 < boundary_idx < n:
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.0, zorder=4)
        ax.text(boundary_idx, triangle_h + extra_top * 0.45, "Active HOR start",
                ha="center", va="bottom", fontsize=46, color=boundary_colour, weight="bold")
        ax.annotate("", xy=(0, -extra_bot * 0.15), xytext=(boundary_idx, -extra_bot * 0.15),
                    arrowprops=dict(arrowstyle="-", lw=3.0, color="#666666"))
        ax.text(boundary_idx / 2.0, -extra_bot * 0.40, "SNPs", ha="center", va="top",
                fontsize=42, color="#444", weight="bold")
        ax.annotate("", xy=(boundary_idx, -extra_bot * 0.15), xytext=(n, -extra_bot * 0.15),
                    arrowprops=dict(arrowstyle="-", lw=3.0, color="#666666"))
        ax.text((boundary_idx + n) / 2.0, -extra_bot * 0.40, "SVs", ha="center", va="top",
                fontsize=42, color="#444", weight="bold")

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.0, zorder=4)
        ax.text(upstream_col, triangle_h + extra_top * 0.45, f"−{window_label}",
                ha="center", va="bottom", fontsize=46, color=upstream_colour, weight="bold")

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=44, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", which="major", length=10, width=2.0, color="#3a3a3a", pad=10)
    ax.set_title(title, fontsize=54, pad=80, weight="bold")
    if payload.get("snp_legend_label") and payload.get("sv_legend_label"):
        sub = f"{payload['snp_legend_label']}    ·    {payload['sv_legend_label']}"
        ax.text(0.5, 1.005, sub, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=32, color="#404040")

    cbar = fig.colorbar(coll, ax=ax, fraction=0.024, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v03(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(4.0, triangle_h * 0.30)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
        ]
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
        ]

        bx, fn = float(boundary_idx), float(n)
        corners = [
            (bx / 2.0, bx / 2.0),
            (fn / 2.0, fn / 2.0),
            ((bx + fn) / 2.0, (fn - bx) / 2.0),
            (bx, 0.0),
        ]
        tick_len = max(3.0, n * 0.025)
        for cx, cy in corners:
            for dx, dy in [(tick_len, tick_len), (-tick_len, tick_len),
                           (tick_len, -tick_len), (-tick_len, -tick_len)]:
                inside_x = cx + dx
                inside_y = cy + dy
                if inside_x < (bx / 2.0) - 0.5 or inside_x > ((bx + fn) / 2.0) + 0.5:
                    continue
                if inside_y < -0.5 or inside_y > fn / 2.0 + 0.5:
                    continue
                ax.plot([cx, inside_x], [cy, inside_y], color=cross_colour,
                        lw=4.0, zorder=5, solid_capstyle="round")
        rhombus_pts = corners + [corners[0]]
        ax.plot([p[0] for p in rhombus_pts], [p[1] for p in rhombus_pts],
                color=cross_colour, lw=1.2, linestyle=(0, (2, 5)), alpha=0.40, zorder=4)
        ax.text((bx + fn) / 2.0 - bx / 6.0, (fn - bx) / 2.0 + bx / 6.0,
                "SV × SNP cross-LD", ha="center", va="bottom",
                fontsize=68, color=cross_colour, weight="bold")

        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=2, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=36, handlelength=2.5, handleheight=1.8,
                        borderpad=0.8, columnspacing=3.0)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold")

    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v04(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    boundary_idx_right = payload.get("boundary_idx_right")
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    downstream_pos = payload.get("downstream_pos")
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(22, min(40, n * 0.165))
    fig_h = max(13, min(24, n * 0.105))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="RdBu_r", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    has_right_hint = boundary_idx_right is not None
    extra_top = max(6.0, triangle_h * (0.55 if has_right_hint else 0.42))
    extra_bot = max(4.0, triangle_h * 0.30)
    xmargin_left = max(5.0, 0.08 * n)
    xmargin_right = max(5.0, 0.08 * n) if has_right_hint else max(1.0, 0.02 * n)
    ax.set_xlim(-xmargin_left, n + xmargin_right)
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")
    n_scale = min(1.0, max(0.72, n / 500.0))

    snp_strip_fill = "#b6c4d4"
    snp_strip_edge = "#2e4960"
    snp_strip_fill_r = "#b6d3ad"
    snp_strip_edge_r = "#2c5435"
    cen_strip_fill = "#eccc8e"
    cen_strip_edge = "#5e3f17"
    boundary_colour = "#1a1a1a"
    upstream_colour = "#0f7b6c"
    downstream_colour = "#a05a0d"
    cross_colour = "#5e3aa6"
    wedge_edge = "#404040"
    wedge_txt = "#2a2a2a"
    halo_effects = [path_effects.Stroke(linewidth=4, foreground="white"),
                    path_effects.Normal()]
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    window_label_r = format_window_label(
        payload.get("actual_window_bp_right", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    has_right = boundary_idx_right is not None and boundary_idx is not None \
        and 0 < boundary_idx < boundary_idx_right < n

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        strip_lw = 3.0
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=strip_lw))
        sv_right_edge = float(boundary_idx_right) if has_right else float(n)
        ax.add_patch(Rectangle((boundary_idx, strip_y), sv_right_edge - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=strip_lw))
        if has_right:
            ax.add_patch(Rectangle(
                (boundary_idx_right, strip_y), n - boundary_idx_right, strip_h,
                facecolor=snp_strip_fill_r, edgecolor=snp_strip_edge_r, lw=strip_lw,
            ))

        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=strip_lw),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=strip_lw),
        ]
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
        ]
        if has_right:
            legend_handles.append(
                Patch(facecolor=snp_strip_fill_r, edgecolor=snp_strip_edge_r, lw=strip_lw)
            )
            legend_labels.append(payload.get("snp_legend_label_right") or "peri-centromeric SNPs (R)")

        bx, fn = float(boundary_idx), float(n)
        bxr = float(boundary_idx_right) if has_right else None

        wedges_to_draw = [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
        ]
        if has_right:
            wedges_to_draw.append(
                [(bx, 0.0), ((bx + bxr) / 2.0, (bxr - bx) / 2.0), (bxr, 0.0), (bx, 0.0)]
            )
            wedges_to_draw.append(
                [(bxr, 0.0), ((bxr + fn) / 2.0, (fn - bxr) / 2.0), (fn, 0.0), (bxr, 0.0)]
            )
        else:
            wedges_to_draw.append(
                [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)]
            )
        for wedge in wedges_to_draw:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=3.0, linestyle=(0, (6, 4)), zorder=5)

        wedge_fs = max(26.0, 38.0 * n_scale)
        if bx >= 18:
            ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                    fontsize=wedge_fs, color=wedge_txt, style="italic", weight="bold", zorder=6)
        sv_width = (bxr - bx) if has_right else (fn - bx)
        sv_cx = (bx + bxr) / 2.0 if has_right else (bx + fn) / 2.0
        if sv_width >= 18:
            ax.text(sv_cx, sv_width * 0.22, "SV × SV", ha="center", va="center",
                    fontsize=wedge_fs, color=wedge_txt, style="italic", weight="bold", zorder=6)
        if has_right and (fn - bxr) >= 18:
            ax.text((bxr + fn) / 2.0, (fn - bxr) * 0.22, "SNP × SNP",
                    ha="center", va="center", fontsize=wedge_fs, color=wedge_txt,
                    style="italic", weight="bold", zorder=6)

        end_L = bxr if has_right else fn
        rhombus_L = [
            (bx / 2.0, bx / 2.0),
            (end_L / 2.0, end_L / 2.0),
            ((bx + end_L) / 2.0, (end_L - bx) / 2.0),
            (bx, 0.0),
            (bx / 2.0, bx / 2.0),
        ]
        ax.plot([p[0] for p in rhombus_L], [p[1] for p in rhombus_L],
                color=cross_colour, lw=4.5, linestyle=(0, (6, 3)), zorder=5)
        cross_cx_L = (2.0 * bx + end_L) / 4.0
        cross_cy_L = end_L / 4.0
        sv_width_for_fs = (bxr - bx) if has_right else (fn - bx)
        cross_fs = min(60.0, max(36.0 * n_scale, min(bx, sv_width_for_fs) * 0.38))
        cross_fs = max(cross_fs, 30.0)
        ax.text(cross_cx_L, cross_cy_L, "SNP × SV", ha="center", va="center",
                fontsize=cross_fs, color=cross_colour, style="italic", weight="bold",
                zorder=7, path_effects=halo_effects)

        if has_right:
            rhombus_R = [
                ((bx + bxr) / 2.0, (bxr - bx) / 2.0),
                ((bx + fn) / 2.0, (fn - bx) / 2.0),
                ((bxr + fn) / 2.0, (fn - bxr) / 2.0),
                (bxr, 0.0),
                ((bx + bxr) / 2.0, (bxr - bx) / 2.0),
            ]
            ax.plot([p[0] for p in rhombus_R], [p[1] for p in rhombus_R],
                    color=cross_colour, lw=4.5, linestyle=(0, (6, 3)), zorder=5)
            cross_cx_R = (bx + 2.0 * bxr + fn) / 4.0
            cross_cy_R = (fn - bx) / 4.0
            cross_fs_R = min(60.0, max(36.0 * n_scale,
                                       min(bxr - bx, fn - bxr) * 0.38))
            cross_fs_R = max(cross_fs_R, 30.0)
            ax.text(cross_cx_R, cross_cy_R, "SV × SNP", ha="center", va="center",
                    fontsize=cross_fs_R, color=cross_colour, style="italic", weight="bold",
                    zorder=7, path_effects=halo_effects)

            apex_w = fn - end_L
            if apex_w >= 14 and bx >= 14:
                apex_poly = [
                    (end_L / 2.0, end_L / 2.0),
                    (fn / 2.0, fn / 2.0),
                    ((bx + fn) / 2.0, (fn - bx) / 2.0),
                    ((bx + end_L) / 2.0, (end_L - bx) / 2.0),
                    (end_L / 2.0, end_L / 2.0),
                ]
                ax.plot([p[0] for p in apex_poly], [p[1] for p in apex_poly],
                        color="#0d0d0d", lw=3.5, linestyle=(0, (2, 2)), zorder=5)
                far_cx = sum(p[0] for p in apex_poly[:-1]) / 4.0
                far_cy = sum(p[1] for p in apex_poly[:-1]) / 4.0
                ax.text(far_cx, far_cy, "SNP$_L$ × SNP$_R$\ndistal cross-LD",
                        ha="center", va="center", fontsize=36, color="#0d0d0d",
                        style="italic", weight="bold",
                        zorder=7, path_effects=halo_effects)

        ax.axvline(float(boundary_idx), color=boundary_colour, lw=4.0, zorder=4)
        hor_y_lo = triangle_h + extra_top * 0.18
        hor_y_mid = triangle_h + extra_top * 0.50
        hor_y_hi = triangle_h + extra_top * 0.82
        hor_fs = max(28.0, 44.0 * n_scale)
        ax.annotate("Active HOR start",
                    xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, hor_y_lo),
                    ha="center", va="bottom", fontsize=hor_fs, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.5, color=boundary_colour),
                    path_effects=halo_effects)
        if has_right:
            ax.axvline(float(boundary_idx_right), color=boundary_colour, lw=4.0, zorder=4)
            ax.annotate("Active HOR end",
                        xy=(boundary_idx_right, triangle_h + extra_top * 0.02),
                        xytext=(boundary_idx_right, hor_y_hi),
                        ha="center", va="bottom", fontsize=hor_fs,
                        color=boundary_colour, weight="bold",
                        arrowprops=dict(arrowstyle="-", lw=2.5, color=boundary_colour),
                        path_effects=halo_effects)

        if has_right and payload.get("expanded_window_note_right"):
            bk_x = float(boundary_idx_right)
            break_fs = max(22.0, 36.0 * n_scale)
            ax.annotate("// AXIS BREAK //",
                        xy=(bk_x, triangle_h + extra_top * 0.05),
                        xytext=(bk_x, hor_y_mid),
                        ha="center", va="bottom", fontsize=break_fs, color="#7a1410",
                        weight="bold", style="italic", zorder=9,
                        arrowprops=dict(arrowstyle="-", lw=2.5, color="#b3211a",
                                        linestyle=(0, (4, 3))),
                        path_effects=halo_effects)

    upstream_col = None
    downstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if has_right and downstream_pos is not None:
        downstream_col = pos_to_column(variants, downstream_pos)

    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=4.0, zorder=4)
        ax.annotate(f"−{window_label} upstream",
                    xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.90),
                    ha="center", va="top", fontsize=42, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.5, color=upstream_colour),
                    path_effects=halo_effects)
    if downstream_col is not None:
        ax.axvline(downstream_col, color=downstream_colour, lw=4.0, zorder=4)
        ax.annotate(f"+{window_label_r} downstream",
                    xy=(downstream_col, -extra_bot * 0.05),
                    xytext=(downstream_col, -extra_bot * 0.90),
                    ha="center", va="top", fontsize=42, color=downstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.5, color=downstream_colour),
                    path_effects=halo_effects)

    if legend_handles:
        ncol = 1
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.04),
                        ncol=ncol, frameon=True, framealpha=0.96, edgecolor="0.35",
                        fontsize=32,
                        handlelength=2.6, handleheight=1.7,
                        borderpad=0.8, labelspacing=0.5)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    centromere_end = int(centromere_span["end"]) if centromere_span else None
    minimal_ticks(
        ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n,
        boundary_idx_right=boundary_idx_right,
        centromere_end=centromere_end,
        downstream_pos=downstream_pos,
        downstream_col=downstream_col,
    )
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    if legend_handles:
        title_pad = 180 if has_right else 140
    else:
        title_pad = 30
    ax.set_title(title, fontsize=50, pad=title_pad, weight="bold")
    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v05(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(22, min(42, n * 0.175))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(3.5, triangle_h * 0.30)
    extra_bot = max(4.0, triangle_h * 0.30)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    wedge_edge = "#404040"
    wedge_txt = "#2a2a2a"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.10)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
            Patch(facecolor="none", edgecolor=cross_colour, lw=3.0, linestyle="--"),
        ]
        legend_labels += [
            payload["snp_legend_label"].replace("  ·  ", "\n", 1),
            payload["sv_legend_label"].replace("  ·  ", "\n", 1),
            "SV × SNP\ncross-LD region",
        ]

        bx, fn = float(boundary_idx), float(n)
        for wedge in [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
            [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)],
        ]:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=3.0, linestyle=(0, (6, 4)), zorder=5)
        ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        ax.text((bx + fn) / 2.0, (fn - bx) * 0.22, "SV × SV", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        rhombus = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                   ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        ax.plot([p[0] for p in rhombus], [p[1] for p in rhombus],
                color=cross_colour, lw=3.0, linestyle=(0, (6, 3)), zorder=5)
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.55),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="upper center", bbox_to_anchor=(0.5, -0.16),
                        ncol=3, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=32, handlelength=2.2, handleheight=1.4,
                        borderpad=0.7, columnspacing=2.5)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    ax.set_title(title, fontsize=52, pad=30, weight="bold")

    cbar = fig.colorbar(coll, ax=ax, fraction=0.024, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v06(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(22, min(36, n * 0.110))
    fig_h = fig_w
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    R = np.triu(r) + np.triu(r, 1).T
    img = ax.imshow(R, cmap="coolwarm", vmin=-1, vmax=1, interpolation="nearest", origin="upper")

    boundary_colour = "#262626"
    cross_colour = "#6c3baa"
    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        ax.axhline(boundary_idx - 0.5, color=boundary_colour, lw=4.0)
        ax.axvline(boundary_idx - 0.5, color=boundary_colour, lw=4.0)
        ax.add_patch(Rectangle((-0.5, -0.5), boundary_idx, boundary_idx,
                               facecolor="none", edgecolor=snp_strip_edge, lw=4.0, zorder=4))
        ax.add_patch(Rectangle((boundary_idx - 0.5, boundary_idx - 0.5),
                               n - boundary_idx, n - boundary_idx,
                               facecolor="none", edgecolor=cen_strip_edge, lw=4.0, zorder=4))
        ax.add_patch(Rectangle((boundary_idx - 0.5, -0.5),
                               n - boundary_idx, boundary_idx,
                               facecolor="none", edgecolor=cross_colour, lw=4.0,
                               linestyle="--", zorder=4))
        ax.add_patch(Rectangle((-0.5, boundary_idx - 0.5),
                               boundary_idx, n - boundary_idx,
                               facecolor="none", edgecolor=cross_colour, lw=4.0,
                               linestyle="--", zorder=4))
        ax.text(boundary_idx / 2.0, boundary_idx / 2.0, "SNP × SNP",
                ha="center", va="center", fontsize=44, color="#222",
                style="italic", weight="bold")
        ax.text((boundary_idx + n) / 2.0, (boundary_idx + n) / 2.0, "SV × SV",
                ha="center", va="center", fontsize=44, color="#222",
                style="italic", weight="bold")
        ax.text((boundary_idx + n) / 2.0, boundary_idx / 2.0,
                "SV × SNP\ncross-LD", ha="center", va="center",
                fontsize=40, color=cross_colour, weight="bold")
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
            Patch(facecolor="none", edgecolor=cross_colour, lw=3.0, linestyle="--"),
        ]
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
            "SV × SNP cross-LD quadrants",
        ]

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.05),
                        ncol=3, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=30, handlelength=2.2, handleheight=1.4,
                        borderpad=0.7, columnspacing=2.5)
        leg.get_frame().set_linewidth(1.5)

    ax.set_xticks([0, boundary_idx, n - 1] if boundary_idx else [0, n - 1])
    centromere_start = int(centromere_span["start"]) if centromere_span else None
    xtick_labels = []
    if upstream_pos is not None:
        xtick_labels.append(f"{upstream_pos / 1e6:.3f}")
    elif n > 0:
        xtick_labels.append(f"{variants.iloc[0]['pos'] / 1e6:.3f}")
    if boundary_idx and centromere_start:
        xtick_labels.append(f"{centromere_start / 1e6:.3f}")
    xtick_labels.append(f"{variants.iloc[-1]['pos'] / 1e6:.3f}")
    ax.set_xticklabels(xtick_labels, fontsize=36)
    ax.set_yticks([])
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=44, weight="semibold")
    ax.set_ylabel("Genomic position (Mb, chm13v2.0)", fontsize=44, weight="semibold")
    title_pad = 80 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold")

    cbar = fig.colorbar(img, ax=ax, fraction=0.040, pad=0.020)
    cbar.set_label("LD (Pearson r)", fontsize=42, weight="semibold")
    cbar.ax.tick_params(labelsize=34)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v07(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(4.0, triangle_h * 0.30)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    wedge_edge = "#404040"
    wedge_txt = "#2a2a2a"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
            Patch(facecolor=cross_colour, alpha=0.25, edgecolor=cross_colour, lw=2.5),
        ]
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
            "SV × SNP cross-LD region",
        ]

        bx, fn = float(boundary_idx), float(n)
        for wedge in [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
            [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)],
        ]:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=3.0, linestyle=(0, (6, 4)), zorder=5)
        ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        ax.text((bx + fn) / 2.0, (fn - bx) * 0.22, "SV × SV", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        rhombus_pts = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                       ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0)]
        v07_alpha = 0.12 if n < 250 else (0.06 if n < 450 else 0.035)
        ax.add_patch(Polygon(rhombus_pts, closed=True, facecolor=cross_colour,
                             alpha=v07_alpha, edgecolor=cross_colour, lw=2.5, zorder=5))
        ax.text((bx + fn) / 2.0 - bx / 6.0, (fn - bx) / 2.0 + bx / 6.0, "SV × SNP cross-LD",
                ha="center", va="bottom", fontsize=68, color=cross_colour, weight="bold")
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=3, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=30, handlelength=2.2, handleheight=1.5,
                        borderpad=0.7, columnspacing=2.5)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold")
    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v08(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(5.5, triangle_h * 0.40)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
        ]
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
        ]

        bx, fn = float(boundary_idx), float(n)
        rhombus = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                   ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        ax.plot([p[0] for p in rhombus], [p[1] for p in rhombus],
                color=cross_colour, lw=3.0, linestyle=(0, (6, 3)), zorder=5)
        ax.text((bx + fn) / 2.0 - bx / 6.0, (fn - bx) / 2.0 + bx / 6.0, "SV × SNP cross-LD",
                ha="center", va="bottom", fontsize=68, color=cross_colour, weight="bold")
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))
        ax.annotate("", xy=(2, -extra_bot * 0.50), xytext=(boundary_idx - 2, -extra_bot * 0.50),
                    arrowprops=dict(arrowstyle="-", lw=3.5, color=snp_strip_edge))
        ax.text(boundary_idx / 2.0, -extra_bot * 0.72, "SNPs", ha="center", va="top",
                fontsize=46, color=snp_strip_edge, weight="bold")
        ax.annotate("", xy=(boundary_idx + 2, -extra_bot * 0.50), xytext=(n - 2, -extra_bot * 0.50),
                    arrowprops=dict(arrowstyle="-", lw=3.5, color=cen_strip_edge))
        ax.text((boundary_idx + n) / 2.0, -extra_bot * 0.72, "SVs", ha="center", va="top",
                fontsize=46, color=cen_strip_edge, weight="bold")

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.20),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=2, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=34, handlelength=2.5, handleheight=1.6,
                        borderpad=0.8, columnspacing=3.0)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=44, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold")
    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold")
    cbar.ax.tick_params(labelsize=36)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




def render_v09(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(20, min(38, n * 0.155))
    fig_h = max(12, min(22, n * 0.10))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="#1c1f26")
    ax.set_facecolor("#1c1f26")
    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(4.0, triangle_h * 0.30)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    fig.set_facecolor("#262a33")
    ax.set_facecolor("#262a33")
    text_colour = "#f5f5f5"
    snp_strip_fill = "#4a4a4a"
    snp_strip_edge = "#e0e0e0"
    cen_strip_fill = "#8a6630"
    cen_strip_edge = "#fadea0"
    boundary_colour = "#ffffff"
    upstream_colour = "#7ff0d8"
    cross_colour = "#cda3ff"
    wedge_edge = "#e0e0e0"
    wedge_txt = "#ffffff"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
        ]
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
        ]

        bx, fn = float(boundary_idx), float(n)
        for wedge in [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
            [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)],
        ]:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=3.0, linestyle=(0, (6, 4)), zorder=5)
        ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        ax.text((bx + fn) / 2.0, (fn - bx) * 0.22, "SV × SV", ha="center", va="center",
                fontsize=42, color=wedge_txt, style="italic", weight="bold", zorder=6)
        rhombus = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                   ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        ax.plot([p[0] for p in rhombus], [p[1] for p in rhombus],
                color=cross_colour, lw=3.5, linestyle=(0, (6, 3)), zorder=5)
        clx, cly = (bx + fn) / 2.0 - bx / 4.0, (fn - bx) / 2.0 + bx / 4.0
        ax.annotate("SV × SNP cross-LD", xy=(clx, cly),
                    xytext=(clx + max(6, n * 0.04), cly + max(4, triangle_h * 0.10)),
                    ha="left", va="bottom", fontsize=68, color=cross_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=cross_colour))
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=44, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=44, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=2, frameon=True, framealpha=0.92, edgecolor=text_colour,
                        facecolor="#1c1f26",
                        fontsize=34, handlelength=2.5, handleheight=1.6,
                        borderpad=0.8, columnspacing=3.0,
                        labelcolor=text_colour)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=46, weight="semibold", color=text_colour)
    ax.set_yticks([])
    ax.tick_params(axis="x", which="major", length=10, width=2.0, color=text_colour, labelcolor=text_colour, pad=10)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(text_colour)
    ax.spines["bottom"].set_linewidth(2.0)
    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=52, pad=title_pad, weight="bold", color=text_colour)
    cbar = fig.colorbar(coll, ax=ax, fraction=0.026, pad=0.012)
    cbar.set_label("LD (Pearson r)", fontsize=44, weight="semibold", color=text_colour)
    cbar.ax.tick_params(labelsize=36, colors=text_colour)
    cbar.outline.set_edgecolor(text_colour)

    fig.tight_layout()
    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35, facecolor=fig.get_facecolor())
    plt.close(fig)




def render_v10(args, payload) -> None:
    variants = payload["variants"]
    r = payload["r"]
    n = len(variants)
    boundary_idx = payload["boundary_idx"]
    centromere_span = payload["centromere_span"]
    upstream_pos = payload["upstream_pos"]
    title = payload["title"]
    out_jpg = payload["out_jpg"]

    fig_w = max(22, min(40, n * 0.165))
    fig_h = max(14, min(24, n * 0.115))
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.05, 0.08, 0.78, 0.82])

    polys, values = build_rotated_polygons(r)
    coll = PolyCollection(polys, array=values, cmap="coolwarm", edgecolors="none", clim=(-1, 1))
    ax.add_collection(coll)

    triangle_h = n / 2.0
    extra_top = max(5.0, triangle_h * 0.42)
    extra_bot = max(4.0, triangle_h * 0.30)
    ax.set_xlim(-max(4.0, 0.07 * n), n + max(1.0, 0.02 * n))
    ax.set_ylim(-extra_bot, triangle_h + extra_top)
    ax.set_aspect("equal")

    snp_strip_fill = "#dadada"
    snp_strip_edge = "#5a5a5a"
    cen_strip_fill = "#e0d2a8"
    cen_strip_edge = "#7a5a25"
    boundary_colour = "#262626"
    upstream_colour = "#2a9d8f"
    cross_colour = "#6c3baa"
    wedge_edge = "#404040"
    wedge_txt = "#2a2a2a"
    window_label = format_window_label(payload.get("actual_window_bp", args.upstream_extension_bp))
    legend_handles: list[Patch] = []
    legend_labels: list[str] = []

    cross_pts = None
    if boundary_idx is not None and 0 < boundary_idx < n:
        strip_y = triangle_h + extra_top * 0.06
        strip_h = max(1.5, extra_top * 0.07)
        ax.add_patch(Rectangle((0, strip_y), boundary_idx, strip_h,
                               facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0))
        ax.add_patch(Rectangle((boundary_idx, strip_y), n - boundary_idx, strip_h,
                               facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0))
        legend_handles += [
            Patch(facecolor=snp_strip_fill, edgecolor=snp_strip_edge, lw=2.0),
            Patch(facecolor=cen_strip_fill, edgecolor=cen_strip_edge, lw=2.0),
        ]
        cen_mb = (int(centromere_span["end"]) - int(centromere_span["start"])) / 1e6 if centromere_span else None
        legend_labels += [
            payload["snp_legend_label"],
            payload["sv_legend_label"],
        ]
        bx, fn = float(boundary_idx), float(n)
        for wedge in [
            [(0.0, 0.0), (bx / 2.0, bx / 2.0), (bx, 0.0), (0.0, 0.0)],
            [(bx, 0.0), ((bx + fn) / 2.0, (fn - bx) / 2.0), (fn, 0.0), (bx, 0.0)],
        ]:
            ax.plot([p[0] for p in wedge], [p[1] for p in wedge],
                    color=wedge_edge, lw=2.6, linestyle=(0, (6, 4)), zorder=5)
        ax.text(bx * 0.5, bx * 0.22, "SNP × SNP", ha="center", va="center",
                fontsize=38, color=wedge_txt, style="italic", weight="bold", zorder=6)
        ax.text((bx + fn) / 2.0, (fn - bx) * 0.22, "SV × SV", ha="center", va="center",
                fontsize=38, color=wedge_txt, style="italic", weight="bold", zorder=6)
        rhombus = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                   ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        ax.plot([p[0] for p in rhombus], [p[1] for p in rhombus],
                color=cross_colour, lw=3.0, linestyle=(0, (6, 3)), zorder=5)
        cross_pts = (bx, fn)
        ax.text((bx + fn) / 2.0 - bx / 4.0, (fn - bx) / 2.0 + bx / 4.0,
                "see inset →", ha="center", va="bottom",
                fontsize=32, color=cross_colour, weight="bold", style="italic")
        ax.axvline(float(boundary_idx), color=boundary_colour, lw=3.5, zorder=4)
        ax.annotate("Active HOR start", xy=(boundary_idx, triangle_h + extra_top * 0.02),
                    xytext=(boundary_idx, triangle_h + extra_top * 0.40),
                    ha="center", va="bottom", fontsize=40, color=boundary_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=boundary_colour))

    upstream_col = None
    if boundary_idx is not None and upstream_pos is not None:
        upstream_col = pos_to_column(variants, upstream_pos)
    if upstream_col is not None:
        ax.axvline(upstream_col, color=upstream_colour, lw=3.5, zorder=4)
        ax.annotate(f"−{window_label} upstream", xy=(upstream_col, -extra_bot * 0.05),
                    xytext=(upstream_col, -extra_bot * 0.85),
                    ha="center", va="top", fontsize=40, color=upstream_colour, weight="bold",
                    arrowprops=dict(arrowstyle="-", lw=2.2, color=upstream_colour))

    if legend_handles:
        leg = ax.legend(legend_handles, legend_labels,
                        loc="lower center", bbox_to_anchor=(0.5, 1.06),
                        ncol=2, frameon=True, framealpha=0.95, edgecolor="0.4",
                        fontsize=32, handlelength=2.5, handleheight=1.6,
                        borderpad=0.7, columnspacing=2.5)
        leg.get_frame().set_linewidth(1.5)

    centromere_start = int(centromere_span["start"]) if centromere_span else None
    minimal_ticks(ax, variants, boundary_idx, centromere_start, upstream_pos, upstream_col, n)
    ax.set_xlabel("Genomic position (Mb, chm13v2.0)", fontsize=42, weight="semibold")
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="x", length=10, width=2.0, pad=10)
    title_pad = 100 if legend_handles else 30
    ax.set_title(title, fontsize=48, pad=title_pad, weight="bold")

    cbar_ax = fig.add_axes([0.86, 0.30, 0.018, 0.45])
    cbar = fig.colorbar(coll, cax=cbar_ax)
    cbar.set_label("LD (Pearson r)", fontsize=40, weight="semibold")
    cbar.ax.tick_params(labelsize=32)

    if cross_pts is not None:
        bx, fn = cross_pts
        inset = fig.add_axes([0.08, 0.70, 0.13, 0.18])
        sub_polys = []
        sub_vals = []
        for i in range(int(bx)):
            for j in range(int(bx), int(fn)):
                corners = [
                    ((j + i) / 2.0, (j - i) / 2.0),
                    ((j + 1 + i) / 2.0, (j + 1 - i) / 2.0),
                    ((j + 1 + i + 1) / 2.0, (j + 1 - i - 1) / 2.0),
                    ((j + i + 1) / 2.0, (j - i - 1) / 2.0),
                ]
                sub_polys.append(corners)
                sub_vals.append(r[i, j])
        coll2 = PolyCollection(sub_polys, array=np.asarray(sub_vals), cmap="coolwarm",
                               edgecolors="none", clim=(-1, 1))
        inset.add_collection(coll2)
        inset.set_xlim(bx / 2.0 - 1, (bx + fn) / 2.0 + 1)
        inset.set_ylim(-1, fn / 2.0 + 1)
        inset.set_aspect("equal")
        rhombus_in = [(bx / 2.0, bx / 2.0), (fn / 2.0, fn / 2.0),
                      ((bx + fn) / 2.0, (fn - bx) / 2.0), (bx, 0.0), (bx / 2.0, bx / 2.0)]
        inset.plot([p[0] for p in rhombus_in], [p[1] for p in rhombus_in],
                   color=cross_colour, lw=2.5, linestyle=(0, (5, 3)))
        for s in ("top", "right", "left", "bottom"):
            inset.spines[s].set_color(cross_colour)
            inset.spines[s].set_linewidth(2.5)
        inset.set_xticks([])
        inset.set_yticks([])
        inset.set_title("SV × SNP cross-LD\n(zoom)", fontsize=30, color=cross_colour,
                        weight="bold", pad=6)

    fig.savefig(out_jpg, dpi=200, format="jpg", bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)




STYLE_RENDERERS = {
    "v01": render_v01,
    "v02": render_v02,
    "v03": render_v03,
    "v04": render_v04,
    "v05": render_v05,
    "v06": render_v06,
    "v07": render_v07,
    "v08": render_v08,
    "v09": render_v09,
    "v10": render_v10,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--style", choices=sorted(STYLE_RENDERERS), required=True)
    parser.add_argument("--chr-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sv-category", default="original_decomposed")
    parser.add_argument("--min-maf", type=float, default=0.0)
    parser.add_argument("--max-variants-for-plot", type=int, default=250)
    parser.add_argument("--snp-vcf", type=Path)
    parser.add_argument("--snp-region")
    parser.add_argument("--snp-min-maf", type=float, default=0.05)
    parser.add_argument("--cross-max-sv-plot", type=int, default=80)
    parser.add_argument("--cross-max-snp-plot", type=int, default=160)
    parser.add_argument(
        "--top-snp-fraction", type=float, default=0.20,
        help="After computing SNP-SV cross-r, keep the top X fraction of SNPs "
        "by mean |r| to SVs (drops the bottom 1-X). Default 0.20.",
    )
    parser.add_argument(
        "--min-snps-to-plot", type=int, default=20,
        help="Floor on number of SNPs kept after the fraction filter.",
    )
    parser.add_argument(
        "--max-snps-to-plot", type=int, default=500,
        help="Ceiling for plottability — at 1 Mb / 10 Mb windows the 20% "
        "subset can hit thousands of SNPs that won't render legibly. Default 500.",
    )
    parser.add_argument("--upstream-extension-bp", type=int, default=100_000)
    parser.add_argument("--active-bed", type=Path, required=True)
    parser.add_argument(
        "--snp-anchor-pos", type=int, default=None,
        help="GRCh38 centromere left-edge bp position for this chrom. When set, "
        "the -window upstream teal-line marker is drawn at "
        "(snp_anchor_pos - upstream_extension_bp) in GRCh38 coords, matching "
        "the 1KG SNP build. Falls back to CHM13 active-HOR-start - window if unset.",
    )
    parser.add_argument(
        "--snp-region-right", default=None,
        help="Optional second 1KG SNP region for the DOWNSTREAM (right) wing of "
        "the centromere. If provided, the figure shows three blocks along x: "
        "[SNP_left | SV | SNP_right]. Format chr:start-end (GRCh38).",
    )
    parser.add_argument(
        "--snp-anchor-pos-right", type=int, default=None,
        help="GRCh38 centromere right-edge bp position. When --snp-region-right "
        "is set, the +window downstream teal line is drawn at "
        "(snp_anchor_pos_right + actual_window_bp).",
    )
    parser.add_argument(
        "--active-bed-right-end", type=int, default=None,
        help="CHM13 active-HOR end (right edge) in bp. Used to label the right "
        "Active-HOR boundary on the figure. If unset, falls back to the active "
        "BED's chrom span end.",
    )
    parser.add_argument("--kinetochore-region")
    parser.add_argument("--allow-sv-only-fallback", action="store_true")
    parser.add_argument(
        "--sv-select", default="maf",
        choices=["maf", "maxr", "meanr", "stratified", "spaced"],
        help="How to choose the plotted SVs. maf = highest minor-allele frequency "
             "(shipped behaviour); maxr / meanr = strongest max / mean |r| to the "
             "candidate flanking SNPs; stratified = evenly spaced across the |r| "
             "ranking, so the full range is represented; spaced = evenly spaced by "
             "position.")
    parser.add_argument(
        "--balance-blocks", type=float, default=0.0,
        help="If > 0, cap each block's keep-count so no block exceeds this ratio "
             "times the middle block, reducing the width disparity between the SNP "
             "and SV blocks. 0 disables (shipped behaviour).")
    parser.add_argument(
        "--cache-payload", type=Path, default=None,
        help="Write the assembled payload (variants + r + geometry) to this .npz "
             "so figure styling can be iterated without re-reading the VCFs.")
    return parser.parse_args()


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

    upstream_pos = None
    snp_panel_loaded = False
    snp_genotypes = snp_variants = None
    overlap_samples: list[str] = []
    fallback_note = None

    expanded_window_note = None
    actual_window_bp = args.upstream_extension_bp
    if args.snp_vcf is not None and args.snp_region:
        chrom_q, coords_q = args.snp_region.split(":")
        start_q, end_q = (int(x) for x in coords_q.replace(",", "").split("-"))
        candidate_windows = [(start_q, end_q, args.upstream_extension_bp)]
        if args.upstream_extension_bp <= 100_000:
            for mult in (2, 5, 10):
                bigger = args.upstream_extension_bp * mult
                new_start = max(1, end_q - bigger + 1)
                candidate_windows.append((new_start, end_q, bigger))

        loaded = False
        for ws, we, wbp in candidate_windows:
            region_try = f"{chrom_q}:{ws}-{we}"
            try:
                snp_genotypes, snp_variants, overlap_samples = load_snp_panel_from_vcf(
                    vcf_path=args.snp_vcf, region=region_try,
                    target_samples=sv_samples, min_maf=args.snp_min_maf,
                )
            except ValueError:
                continue
            if len(snp_variants) >= args.min_snps_to_plot:
                snp_panel_loaded = True
                actual_window_bp = wbp
                if wbp != args.upstream_extension_bp:
                    expanded_window_note = (
                        f"requested {format_window_label(args.upstream_extension_bp)} window "
                        f"had a 1KG callability gap; auto-expanded upstream to "
                        f"{format_window_label(wbp)} to find SNPs"
                    )
                loaded = True
                break
        if not loaded:
            if not args.allow_sv_only_fallback:
                raise ValueError(f"No 1KG SNPs at any expansion of {args.snp_region}")
            fallback_note = (
                f"1KG SNP panel empty across all expansions of {args.snp_region}; "
                "SV-only triangle (acrocentric / GRCh38 callability gap)."
            )

    snp_genotypes_right = snp_variants_right = None
    overlap_samples_right: list[str] = []
    snp_right_loaded = False
    actual_window_bp_right = args.upstream_extension_bp
    expanded_window_note_right = None
    right_reanchored_pos = None
    if args.snp_region_right is not None:
        chrom_q_r, coords_q_r = args.snp_region_right.split(":")
        start_q_r, end_q_r = (int(x) for x in coords_q_r.replace(",", "").split("-"))

        region_try = f"{chrom_q_r}:{start_q_r}-{end_q_r}"
        try:
            snp_genotypes_right, snp_variants_right, overlap_samples_right = (
                load_snp_panel_from_vcf(
                    vcf_path=args.snp_vcf, region=region_try,
                    target_samples=sv_samples, min_maf=args.snp_min_maf,
                )
            )
            if len(snp_variants_right) >= args.min_snps_to_plot:
                snp_right_loaded = True
                actual_window_bp_right = args.upstream_extension_bp
        except ValueError:
            pass

        if not snp_right_loaded:
            broad_region = f"{chrom_q_r}:{start_q_r}-{start_q_r + 250_000_000}"
            first_pos = find_first_snp_pos(args.snp_vcf, broad_region)
            if first_pos is not None:
                new_start = first_pos
                new_end = first_pos + args.upstream_extension_bp - 1
                region_try = f"{chrom_q_r}:{new_start}-{new_end}"
                try:
                    snp_genotypes_right, snp_variants_right, overlap_samples_right = (
                        load_snp_panel_from_vcf(
                            vcf_path=args.snp_vcf, region=region_try,
                            target_samples=sv_samples, min_maf=args.snp_min_maf,
                        )
                    )
                    if len(snp_variants_right) >= args.min_snps_to_plot:
                        snp_right_loaded = True
                        actual_window_bp_right = args.upstream_extension_bp
                        right_reanchored_pos = first_pos
                        expanded_window_note_right = (
                            f"strict {format_window_label(args.upstream_extension_bp)} "
                            f"downstream of cent-right was empty in 1KG; "
                            f"re-anchored to first callable SNP at "
                            f"{chrom_q_r}:{first_pos:,} (gap "
                            f"{(first_pos - (start_q_r - 1)) / 1e6:.2f} Mb)"
                        )
                except ValueError:
                    pass

    n_snp_total = 0
    n_snp_kept = 0
    n_snp_total_right = 0
    n_snp_kept_right = 0
    n_sv_total = int(len(sv_variants))
    n_sv_kept = 0
    mean_cross_r_kept = None
    mean_cross_r_kept_right = None
    boundary_idx_right = None
    if snp_panel_loaded:
        sv_idx_overlap = np.array([sv_samples.index(s) for s in overlap_samples], dtype=int)
        sv_overlap_genotypes = sv_genotypes[:, sv_idx_overlap]
        n_snp_total = int(len(snp_variants))
        n_snp_total_right_pre = int(len(snp_variants_right)) if snp_right_loaded else 0

        def _frac_keep(total):
            if total <= 0:
                return 0
            k = max(args.min_snps_to_plot, int(round(total * args.top_snp_fraction)))
            return min(k, total, args.max_snps_to_plot)

        k_sv = min(args.cross_max_sv_plot, int(len(sv_variants)))
        k_l = _frac_keep(n_snp_total)
        k_r = _frac_keep(n_snp_total_right_pre)

        if args.balance_blocks and args.balance_blocks > 0:
            ratio = float(args.balance_blocks)
            desired = [k_l, k_sv, k_r] if snp_right_loaded else [k_l, k_sv]
            avail = ([n_snp_total, int(len(sv_variants)), n_snp_total_right_pre]
                     if snp_right_loaded else [n_snp_total, int(len(sv_variants))])
            mid = sorted(d for d in desired if d > 0)[len(
                [d for d in desired if d > 0]) // 2]
            lo, hi = mid / ratio, mid * ratio
            bal = [int(min(av, max(lo, min(hi, d)))) for d, av in zip(desired, avail)]
            bal = [max(1, b) for b in bal]
            if snp_right_loaded:
                k_l, k_sv, k_r = bal
            else:
                k_l, k_sv = bal
            print(f"[balance] ratio={ratio} desired={desired} -> keep={bal}")

        if args.sv_select == "maf":
            plot_sv, sv_plot_idx = select_plot_subset(sv_variants, k_sv)
        else:
            blocks = [np.abs(compute_cross_r(snp_genotypes, sv_overlap_genotypes))]
            if snp_right_loaded:
                blocks.append(np.abs(compute_cross_r(snp_genotypes_right,
                                                     sv_overlap_genotypes)))
            A = np.vstack(blocks)
            with np.errstate(invalid="ignore"):
                if args.sv_select == "meanr":
                    score = np.nanmean(A, axis=0)
                else:
                    score = np.nanmax(A, axis=0)
            score = np.nan_to_num(score, nan=-1.0)
            if args.sv_select == "spaced":
                pos = sv_variants["pos"].to_numpy()
                order = np.argsort(pos)
                take = np.linspace(0, len(order) - 1, min(k_sv, len(order))).round().astype(int)
                sel = order[np.unique(take)]
            elif args.sv_select == "stratified":
                order = np.argsort(-score)
                take = np.linspace(0, len(order) - 1, min(k_sv, len(order))).round().astype(int)
                sel = order[np.unique(take)]
            else:
                sel = np.argsort(-score)[:k_sv]
            sv_plot_idx = np.sort(np.asarray(sel, dtype=int))
            plot_sv = sv_variants.iloc[sv_plot_idx].reset_index(drop=True)
            print(f"[sv-select] {args.sv_select}: kept {len(sv_plot_idx)} of "
                  f"{len(sv_variants)}; score range "
                  f"{score[sv_plot_idx].min():.3f}-{score[sv_plot_idx].max():.3f}")
        sv_genos_for_cross = sv_overlap_genotypes[sv_plot_idx, :]
        n_sv_kept = int(len(plot_sv))
        cross_r_all = compute_cross_r(snp_genotypes, sv_genos_for_cross)
        snp_scores = np.mean(np.abs(cross_r_all), axis=1)
        n_keep = min(k_l, n_snp_total)
        score_order = np.argsort(-snp_scores)
        keep_idx = np.sort(score_order[:n_keep])
        plot_snp = snp_variants.iloc[keep_idx].reset_index(drop=True)
        plot_snp["cross_r_mean_abs"] = snp_scores[keep_idx]
        n_snp_kept = int(len(plot_snp))
        mean_cross_r_kept = float(np.mean(snp_scores[keep_idx]))

        plot_snp_right = None
        snp_genos_right_kept = None
        if snp_right_loaded:
            if overlap_samples_right != overlap_samples:
                raise RuntimeError(
                    "Left/right SNP panels have different sample overlaps with "
                    "the SV matrix; cannot build a combined LD matrix"
                )
            n_snp_total_right = int(len(snp_variants_right))
            cross_r_all_r = compute_cross_r(snp_genotypes_right, sv_genos_for_cross)
            snp_scores_r = np.mean(np.abs(cross_r_all_r), axis=1)
            n_keep_r = min(k_r, n_snp_total_right)
            order_r = np.argsort(-snp_scores_r)
            keep_idx_r = np.sort(order_r[:n_keep_r])
            plot_snp_right = snp_variants_right.iloc[keep_idx_r].reset_index(drop=True)
            plot_snp_right["cross_r_mean_abs"] = snp_scores_r[keep_idx_r]
            n_snp_kept_right = int(len(plot_snp_right))
            mean_cross_r_kept_right = float(np.mean(snp_scores_r[keep_idx_r]))
            snp_genos_right_kept = snp_genotypes_right[keep_idx_r, :]

        if plot_snp_right is not None:
            combined_variants = pd.concat(
                [plot_snp, plot_sv.reset_index(drop=True), plot_snp_right],
                ignore_index=True,
            )
            combined_genotypes = np.vstack([
                snp_genotypes[keep_idx, :],
                sv_genos_for_cross,
                snp_genos_right_kept,
            ])
            boundary_idx = len(plot_snp)
            boundary_idx_right = boundary_idx + len(plot_sv)
            title = f"{chrom}  peri-SNP (L + R) + centromeric SV LD"
        else:
            combined_variants = pd.concat(
                [plot_snp, plot_sv.reset_index(drop=True)], ignore_index=True,
            )
            combined_genotypes = np.vstack(
                [snp_genotypes[keep_idx, :], sv_genos_for_cross]
            )
            boundary_idx = len(plot_snp)
            title = f"{chrom}  peri-SNP + centromeric SV LD"
        combined_r = compute_r(combined_genotypes)
    else:
        plot_sv, sv_plot_idx = select_plot_subset(sv_variants, args.cross_max_sv_plot)
        combined_variants = plot_sv.reset_index(drop=True)
        combined_r = compute_r(sv_genotypes[sv_plot_idx, :])
        boundary_idx = None
        n_sv_kept = int(len(plot_sv))
        title = f"{chrom}  centromeric SV LD  (SV-only)"

    if args.snp_anchor_pos is not None:
        upstream_pos = int(args.snp_anchor_pos) - actual_window_bp
    elif active_span is not None:
        upstream_pos = int(active_span["start"]) - actual_window_bp

    downstream_pos = None
    if boundary_idx_right is not None:
        if right_reanchored_pos is not None:
            downstream_pos = right_reanchored_pos + actual_window_bp_right
        elif args.snp_anchor_pos_right is not None:
            downstream_pos = int(args.snp_anchor_pos_right) + actual_window_bp_right
        elif active_span is not None:
            downstream_pos = int(active_span["end"]) + actual_window_bp_right

    window_lbl = format_window_label(actual_window_bp)
    requested_lbl = format_window_label(args.upstream_extension_bp)
    window_phrase = (
        f"{window_lbl} upstream (auto-expanded from {requested_lbl})"
        if actual_window_bp != args.upstream_extension_bp else
        f"{window_lbl} upstream"
    )
    window_lbl_r = format_window_label(actual_window_bp_right)
    window_phrase_r = (
        f"{window_lbl_r} downstream (auto-expanded from {requested_lbl})"
        if actual_window_bp_right != args.upstream_extension_bp else
        f"{window_lbl_r} downstream"
    )
    cen_mb_val = ((int(active_span["end"]) - int(active_span["start"])) / 1e6
                  if active_span else None)
    if boundary_idx is not None:
        pct = int(round(args.top_snp_fraction * 100))
        snp_legend_label = (
            f"peri-centromeric SNPs (L)  ·  {window_phrase}  ·  "
            f"n = {n_snp_kept}/{n_snp_total} (top {pct}% by mean |r| to SVs, MAF ≥ {args.snp_min_maf:g})"
        )
        if cen_mb_val is not None:
            sv_legend_label = (
                f"Active HOR (centromere)  ·  {cen_mb_val:.2f} Mb  ·  "
                f"n = {n_sv_kept}/{n_sv_total} SVs"
            )
        else:
            sv_legend_label = f"Active HOR (centromeric SVs)  ·  n = {n_sv_kept}/{n_sv_total}"
    else:
        snp_legend_label = None
        if cen_mb_val is not None:
            sv_legend_label = (
                f"Active HOR (centromere)  ·  {cen_mb_val:.2f} Mb  ·  "
                f"n = {n_sv_kept} of {n_sv_total} SVs"
            )
        else:
            sv_legend_label = f"Active HOR (centromeric SVs)  ·  n = {n_sv_kept} of {n_sv_total}"

    snp_legend_label_right = None
    if boundary_idx_right is not None:
        pct = int(round(args.top_snp_fraction * 100))
        if right_reanchored_pos is not None:
            gap_mb = ((right_reanchored_pos - int(args.snp_anchor_pos_right or 0))
                      / 1e6 if args.snp_anchor_pos_right else 0)
            window_phrase_r = (
                f"{window_lbl_r} downstream (re-anchored, gap {gap_mb:.1f} Mb)"
            )
        snp_legend_label_right = (
            f"peri-centromeric SNPs (R)  ·  {window_phrase_r}  ·  "
            f"n = {n_snp_kept_right}/{n_snp_total_right} (top {pct}% by mean |r| to SVs, MAF ≥ {args.snp_min_maf:g})"
        )

    out_prefix = args.output_dir / f"{prefix}_{args.style}"
    payload = {
        "variants": combined_variants,
        "r": combined_r,
        "boundary_idx": boundary_idx,
        "boundary_idx_right": boundary_idx_right,
        "centromere_span": active_span,
        "upstream_pos": upstream_pos,
        "downstream_pos": downstream_pos,
        "actual_window_bp": actual_window_bp,
        "actual_window_bp_right": actual_window_bp_right,
        "title": title,
        "out_jpg": Path(f"{out_prefix}.jpg"),
        "fallback_note": fallback_note,
        "snp_legend_label": snp_legend_label,
        "snp_legend_label_right": snp_legend_label_right,
        "sv_legend_label": sv_legend_label,
        "expanded_window_note": expanded_window_note,
        "expanded_window_note_right": expanded_window_note_right,
    }
    if args.cache_payload is not None:
        args.cache_payload.parent.mkdir(parents=True, exist_ok=True)
        _scalars = {k: v for k, v in payload.items()
                    if k not in ("variants", "r", "out_jpg")}
        _scalars["out_jpg"] = str(payload["out_jpg"])
        np.savez_compressed(
            args.cache_payload,
            r=payload["r"],
            variant_pos=combined_variants["pos"].to_numpy(),
            variant_id=combined_variants.get(
                "id", pd.Series([""] * len(combined_variants))).astype(str).to_numpy(),
            variant_kind=combined_variants.get(
                "kind", pd.Series([""] * len(combined_variants))).astype(str).to_numpy(),
            scalars=json.dumps(_scalars),
        )
        print(f"[cache] wrote {args.cache_payload}")

    STYLE_RENDERERS[args.style](args, payload)

    summary = {
        "style": args.style,
        "prefix": prefix,
        "chrom": chrom,
        "snp_min_maf": args.snp_min_maf,
        "top_snp_fraction": args.top_snp_fraction,
        "sv_select": args.sv_select,
        "balance_blocks": args.balance_blocks,
        "n_snp_total": n_snp_total,
        "n_snp_kept": n_snp_kept,
        "n_snp_total_right": n_snp_total_right,
        "n_snp_kept_right": n_snp_kept_right,
        "n_sv_total": n_sv_total,
        "n_sv_kept": n_sv_kept,
        "mean_cross_r_of_kept_snps": mean_cross_r_kept,
        "mean_cross_r_of_kept_snps_right": mean_cross_r_kept_right,
        "upstream_extension_bp": args.upstream_extension_bp,
        "actual_window_bp": actual_window_bp,
        "actual_window_bp_right": actual_window_bp_right,
        "active_span": active_span,
        "fallback_note": fallback_note,
        "right_reanchored_pos": right_reanchored_pos,
        "expanded_window_note_right": expanded_window_note_right,
        "output_jpg": str(payload["out_jpg"]),
    }
    Path(f"{out_prefix}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
