#!/usr/bin/env python3
"""Render the scTHREAD graphical abstract in the Human-scATAC-Corpus layout language.

The figure is deliberately schematic: icon-like matrices, cell clusters and interface
cards explain database functions, while all quantitative labels are taken from the
2026-07-26 frozen catalog scope.
"""

from __future__ import annotations

import os

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import (
    Circle,
    Ellipse,
    FancyArrowPatch,
    FancyBboxPatch,
    PathPatch,
    Polygon,
    Rectangle,
)
from matplotlib.path import Path as MplPath


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STEM = ROOT / "figures" / "scTHREAD_graphical_abstract_scatac_style_v1"
ARIAL = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))

W_MM = 183.0
H_MM = 78.0

PURPLE = "#5B4B9A"
PURPLE_DARK = "#463879"
TEAL = "#0B7C86"
BLUE = "#2F6FB2"
CORAL = "#E34A33"
ORANGE = "#D88B1B"
GREEN = "#4C956C"
INK = "#202124"
MID = "#6D6E73"
LINE = "#C9CBD2"
WHITE = "#FFFFFF"
LAVENDER = "#F4F1FA"
PALE_BLUE = "#EFF6FB"
PALE_TEAL = "#EEF8F8"
PALE_PEACH = "#FFF3EA"
PALE_GOLD = "#FFF8E9"
PALE_GREEN = "#F0F7EF"


def setup_style() -> None:
    if ARIAL.exists():
        font_manager.fontManager.addfont(str(ARIAL))
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 450,
            "figure.dpi": 150,
            "axes.grid": False,
            "savefig.transparent": True,
        }
    )


def add_text(
    ax: plt.Axes,
    x: float,
    y: float,
    value: str,
    *,
    size: float = 6.0,
    color: str = INK,
    weight: str = "normal",
    ha: str = "center",
    va: str = "center",
    rotation: float = 0,
    zorder: int = 20,
) -> None:
    ax.text(
        x,
        y,
        value,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        rotation=rotation,
        family="Arial",
        clip_on=False,
        zorder=zorder,
    )


def rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face: str = WHITE,
    edge: str = PURPLE,
    lw: float = 0.9,
    radius: float = 2.2,
    zorder: int = 2,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.15,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = PURPLE,
    lw: float = 1.1,
    scale: float = 10.0,
    style: str = "-|>",
    zorder: int = 12,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=scale,
            color=color,
            linewidth=lw,
            shrinkA=0,
            shrinkB=0,
            zorder=zorder,
        )
    )


def line(
    ax: plt.Axes,
    xs: list[float],
    ys: list[float],
    *,
    color: str = INK,
    lw: float = 0.7,
    ls: str = "-",
    zorder: int = 10,
) -> None:
    ax.plot(xs, ys, color=color, linewidth=lw, linestyle=ls, zorder=zorder)


def draw_server(ax: plt.Axes, x: float, y: float, scale: float = 1.0) -> None:
    for i in range(3):
        yy = y + (2 - i) * 2.0 * scale
        ax.add_patch(
            FancyBboxPatch(
                (x, yy),
                8.0 * scale,
                1.45 * scale,
                boxstyle=f"round,pad=0.06,rounding_size={0.25 * scale}",
                facecolor=PURPLE,
                edgecolor=PURPLE,
                linewidth=0.5,
                zorder=8,
            )
        )
        ax.add_patch(
            Rectangle(
                (x + 0.75 * scale, yy + 0.47 * scale),
                0.52 * scale,
                0.52 * scale,
                facecolor=WHITE,
                edgecolor="none",
                zorder=9,
            )
        )
        ax.add_patch(
            Rectangle(
                (x + 1.65 * scale, yy + 0.47 * scale),
                0.52 * scale,
                0.52 * scale,
                facecolor=WHITE,
                edgecolor="none",
                zorder=9,
            )
        )
    ax.add_patch(
        Circle(
            (x + 7.8 * scale, y + 1.2 * scale),
            1.75 * scale,
            facecolor=WHITE,
            edgecolor=PURPLE,
            linewidth=1.1,
            zorder=10,
        )
    )
    line(
        ax,
        [x + 9.0 * scale, x + 10.8 * scale],
        [y + 0.0 * scale, y - 1.8 * scale],
        color=PURPLE,
        lw=1.3,
        zorder=10,
    )


def draw_human(ax: plt.Axes, x: float, y: float, scale: float = 1.0) -> None:
    ax.add_patch(
        Circle(
            (x, y + 3.4 * scale),
            1.25 * scale,
            facecolor=PURPLE_DARK,
            edgecolor="none",
            zorder=8,
        )
    )
    body = MplPath(
        [
            (x - 3.0 * scale, y - 2.6 * scale),
            (x - 2.6 * scale, y + 1.1 * scale),
            (x - 1.4 * scale, y + 2.2 * scale),
            (x, y + 2.4 * scale),
            (x + 1.4 * scale, y + 2.2 * scale),
            (x + 2.6 * scale, y + 1.1 * scale),
            (x + 3.0 * scale, y - 2.6 * scale),
            (x - 3.0 * scale, y - 2.6 * scale),
        ],
        [
            MplPath.MOVETO,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CURVE3,
            MplPath.CLOSEPOLY,
        ],
    )
    ax.add_patch(PathPatch(body, facecolor=PURPLE_DARK, edgecolor="none", zorder=8))


def draw_mouse(ax: plt.Axes, x: float, y: float, scale: float = 1.0) -> None:
    ax.add_patch(
        Ellipse(
            (x, y),
            6.2 * scale,
            3.3 * scale,
            facecolor=WHITE,
            edgecolor=CORAL,
            linewidth=0.9,
            zorder=8,
        )
    )
    ax.add_patch(
        Circle(
            (x + 3.2 * scale, y + 0.1 * scale),
            1.25 * scale,
            facecolor=WHITE,
            edgecolor=CORAL,
            linewidth=0.9,
            zorder=8,
        )
    )
    for dx, dy in [(2.8, 1.25), (3.55, 1.15)]:
        ax.add_patch(
            Circle(
                (x + dx * scale, y + dy * scale),
                0.55 * scale,
                facecolor=WHITE,
                edgecolor=CORAL,
                linewidth=0.8,
                zorder=8,
            )
        )
    ax.add_patch(
        Circle(
            (x + 4.05 * scale, y + 0.25 * scale),
            0.12 * scale,
            facecolor=INK,
            edgecolor="none",
            zorder=9,
        )
    )
    line(
        ax,
        [x - 3.0 * scale, x - 4.3 * scale, x - 5.3 * scale, x - 4.3 * scale],
        [y, y - 1.0 * scale, y - 0.7 * scale, y + 0.15 * scale],
        color=CORAL,
        lw=0.9,
        zorder=8,
    )


def draw_cell_cluster(
    ax: plt.Axes,
    x: float,
    y: float,
    *,
    color: str,
    scale: float = 1.0,
) -> None:
    offsets = [
        (-1.5, 0.0),
        (0.0, 0.0),
        (1.5, 0.0),
        (-0.75, 1.35),
        (0.75, 1.35),
        (0.0, 2.65),
    ]
    for dx, dy in offsets:
        ax.add_patch(
            Circle(
                (x + dx * scale, y + dy * scale),
                0.68 * scale,
                facecolor=WHITE,
                edgecolor=color,
                linewidth=0.8,
                zorder=9,
            )
        )
        ax.add_patch(
            Circle(
                (x + dx * scale, y + dy * scale),
                0.24 * scale,
                facecolor=color,
                edgecolor="none",
                alpha=0.75,
                zorder=10,
            )
        )


def draw_transcript(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    *,
    color: str = PURPLE,
    variant: int = 0,
    lw: float = 0.8,
) -> None:
    exon_w = width * 0.16
    gap = (width - 3 * exon_w) / 2
    xs = [x, x + exon_w + gap, x + 2 * (exon_w + gap)]
    if variant == 1:
        xs[1] += gap * 0.15
    for i, xx in enumerate(xs):
        height = 1.2 if not (variant == 1 and i == 1) else 0.75
        ax.add_patch(
            Rectangle(
                (xx, y - height / 2),
                exon_w,
                height,
                facecolor=WHITE if i == 1 else PALE_TEAL,
                edgecolor=color,
                linewidth=lw,
                zorder=9,
            )
        )
    line(
        ax,
        [xs[0] + exon_w, xs[1], xs[1] + exon_w, xs[2]],
        [y, y, y, y],
        color=color,
        lw=lw,
        zorder=8,
    )


def draw_matrix(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    rows: int = 4,
    cols: int = 8,
    palette: tuple[str, ...] = (PALE_BLUE, "#C9E0F0", "#A8CCE5"),
) -> None:
    cell_w = width / cols
    cell_h = height / rows
    for r in range(rows):
        for c in range(cols):
            color = palette[(r * 3 + c * 5 + (r + c) // 3) % len(palette)]
            ax.add_patch(
                Rectangle(
                    (x + c * cell_w, y + (rows - 1 - r) * cell_h),
                    cell_w,
                    cell_h,
                    facecolor=color,
                    edgecolor=INK,
                    linewidth=0.18,
                    zorder=8,
                )
            )


def draw_search_icon(ax: plt.Axes, x: float, y: float, color: str = PURPLE) -> None:
    ax.add_patch(
        Circle(
            (x, y),
            1.45,
            facecolor="none",
            edgecolor=color,
            linewidth=1.0,
            zorder=9,
        )
    )
    line(ax, [x + 1.0, x + 2.4], [y - 1.0, y - 2.4], color=color, lw=1.1, zorder=9)


def draw_folder(ax: plt.Axes, x: float, y: float, color: str) -> None:
    points = [
        (x, y),
        (x + 2.2, y),
        (x + 3.0, y + 1.0),
        (x + 7.0, y + 1.0),
        (x + 7.0, y - 3.4),
        (x, y - 3.4),
    ]
    ax.add_patch(
        Polygon(
            points,
            closed=True,
            facecolor=WHITE,
            edgecolor=color,
            linewidth=0.9,
            zorder=8,
        )
    )


def draw_download_icon(ax: plt.Axes, x: float, y: float, color: str = PURPLE) -> None:
    line(ax, [x, x], [y + 2.0, y - 1.0], color=color, lw=1.1, zorder=9)
    arrow(ax, (x, y + 0.2), (x, y - 2.0), color=color, lw=1.0, scale=8.0)
    line(ax, [x - 2.2, x - 2.2, x + 2.2, x + 2.2], [y - 1.6, y - 3.0, y - 3.0, y - 1.6], color=color, lw=0.9)


def draw_system_icon(ax: plt.Axes, x: float, y: float, kind: str, color: str) -> None:
    ax.add_patch(Circle((x, y), 2.1, facecolor="#F2F2F2", edgecolor="none", zorder=5))
    if kind == "blood":
        for dx, dy, angle in [(-0.7, 0.4, 20), (0.8, 0.3, -15), (0.0, -0.8, 8)]:
            ax.add_patch(
                Ellipse(
                    (x + dx, y + dy),
                    1.35,
                    0.68,
                    angle=angle,
                    facecolor="#F6A6A0",
                    edgecolor=color,
                    linewidth=0.7,
                    zorder=8,
                )
            )
    elif kind == "brain":
        for dx, dy in [(-0.65, 0.35), (0.15, 0.7), (0.75, 0.25), (-0.2, -0.4), (0.65, -0.55)]:
            ax.add_patch(
                Circle(
                    (x + dx, y + dy),
                    0.72,
                    facecolor="#F7B1C2",
                    edgecolor=color,
                    linewidth=0.55,
                    zorder=8,
                )
            )
    elif kind == "heart":
        vertices = [
            (x, y - 1.35),
            (x - 2.4, y + 0.1),
            (x - 1.7, y + 2.0),
            (x, y + 0.75),
            (x + 1.7, y + 2.0),
            (x + 2.4, y + 0.1),
            (x, y - 1.35),
            (x, y - 1.35),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor="#EE756A", edgecolor=color, linewidth=0.7, zorder=8))
    elif kind == "cancer":
        for dx, dy in [(-0.8, 0.6), (0.3, 0.8), (0.9, -0.1), (-0.2, -0.3), (-0.9, -0.6), (0.4, -0.9)]:
            ax.add_patch(
                Circle(
                    (x + dx, y + dy),
                    0.68,
                    facecolor="#D7B4D8",
                    edgecolor=color,
                    linewidth=0.6,
                    zorder=8,
                )
            )
    elif kind == "development":
        line(ax, [x - 1.3, x, x + 1.25], [y + 0.9, y, y + 1.0], color=color, lw=0.7)
        line(ax, [x, x + 1.25], [y, y - 1.0], color=color, lw=0.7)
        for dx, dy in [(-1.3, 0.9), (0, 0), (1.25, 1.0), (1.25, -1.0)]:
            ax.add_patch(Circle((x + dx, y + dy), 0.38, facecolor=WHITE, edgecolor=color, linewidth=0.7, zorder=8))
    else:
        ax.add_patch(Circle((x, y), 1.25, facecolor=WHITE, edgecolor=color, linewidth=0.8, zorder=8))
        line(ax, [x - 0.65, x - 0.1, x + 0.8], [y, y - 0.65, y + 0.65], color=color, lw=1.0, zorder=9)


def draw_coverage_panel(ax: plt.Axes) -> None:
    x, y, w, h = 3.0, 36.0, 45.0, 32.0
    rounded_box(ax, x, y, w, h, face=WHITE)
    add_text(ax, x + 2.2, y + h - 2.0, "Catalog coverage", size=8.2, weight="bold", ha="left", va="top")

    draw_server(ax, x + 3.5, y + 20.3, 0.62)
    add_text(ax, x + 7.1, y + 15.9, "469", size=10.0, color=PURPLE_DARK, weight="bold")
    add_text(ax, x + 7.1, y + 13.1, "samples", size=5.1, color=MID, weight="bold")
    add_text(ax, x + 16.8, y + 15.9, "31", size=10.0, color=TEAL, weight="bold")
    add_text(ax, x + 16.8, y + 13.1, "studies", size=5.1, color=MID, weight="bold")
    draw_human(ax, x + 6.6, y + 7.6, 0.55)
    draw_mouse(ax, x + 15.8, y + 7.8, 0.45)
    add_text(ax, x + 11.2, y + 3.7, "Human · mouse", size=5.1, color=PURPLE_DARK, weight="bold")
    add_text(ax, x + 11.2, y + 1.7, "ONT · PacBio", size=4.9, color=MID)

    icon_specs = [
        ("blood", "blood /\nmarrow", CORAL),
        ("brain", "brain", PURPLE),
        ("heart", "heart", CORAL),
        ("cancer", "cancer", PURPLE),
        ("development", "development", TEAL),
        ("methods", "methods", TEAL),
    ]
    icon_xs = [x + 25.7, x + 32.5, x + 39.3]
    icon_ys = [y + 22.2, y + 13.9]
    for idx, (kind, label, color) in enumerate(icon_specs):
        col, row = idx % 3, idx // 3
        ix, iy = icon_xs[col], icon_ys[row]
        draw_system_icon(ax, ix, iy, kind, color)
        add_text(ax, ix, iy - 3.25, label, size=3.55, color=INK)


def draw_database_panel(ax: plt.Axes) -> None:
    x, y, w, h = 50.0, 36.0, 36.0, 32.0
    rounded_box(ax, x, y, w, h, face=WHITE)

    cyl_x, cyl_y = x + 3.0, y + 3.0
    cyl_w, cyl_h = w - 6.0, h - 6.0
    ax.add_patch(
        Rectangle(
            (cyl_x, cyl_y + 2.6),
            cyl_w,
            cyl_h - 5.2,
            facecolor=PALE_PEACH,
            edgecolor=PURPLE_DARK,
            linewidth=0.75,
            zorder=5,
        )
    )
    ax.add_patch(
        Ellipse(
            (cyl_x + cyl_w / 2, cyl_y + cyl_h - 2.6),
            cyl_w,
            5.2,
            facecolor=PALE_PEACH,
            edgecolor=PURPLE_DARK,
            linewidth=0.75,
            zorder=7,
        )
    )
    ax.add_patch(
        Ellipse(
            (cyl_x + cyl_w / 2, cyl_y + 2.6),
            cyl_w,
            5.2,
            facecolor=PALE_PEACH,
            edgecolor=PURPLE_DARK,
            linewidth=0.75,
            zorder=7,
        )
    )
    # Hide the top half of the bottom ellipse so the cylinder reads as a solid object.
    ax.add_patch(
        Rectangle(
            (cyl_x + 0.2, cyl_y + 2.6),
            cyl_w - 0.4,
            2.8,
            facecolor=PALE_PEACH,
            edgecolor="none",
            zorder=8,
        )
    )
    add_text(ax, x + w / 2, y + h - 4.0, "Database", size=8.0, color=INK)
    add_text(ax, x + w / 2, y + 19.3, "scTHREAD", size=13.0, color=PURPLE, weight="bold")
    draw_transcript(ax, x + 10.6, y + 15.7, 14.8, color=TEAL, variant=1, lw=0.75)
    add_text(ax, x + w / 2, y + 11.0, "845,781", size=10.5, color=PURPLE_DARK, weight="bold")
    add_text(ax, x + w / 2, y + 8.4, "cells", size=6.4, color=PURPLE_DARK, weight="bold")
    add_text(ax, x + w / 2, y + 4.2, ">200k isoforms", size=6.2, color=TEAL, weight="bold")


def draw_evidence_icon(ax: plt.Axes, x: float, y: float, kind: str, color: str) -> None:
    if kind == "expression":
        heights = [1.2, 2.3, 3.2, 1.8]
        for i, hh in enumerate(heights):
            ax.add_patch(Rectangle((x + i * 1.25, y - 1.7), 0.78, hh, facecolor=color, edgecolor="none", zorder=9))
        line(ax, [x - 0.4, x + 5.0], [y - 1.7, y - 1.7], color=color, lw=0.7)
    elif kind == "isoform":
        draw_transcript(ax, x - 0.2, y + 0.7, 5.5, color=color, variant=0)
        draw_transcript(ax, x - 0.2, y - 1.1, 5.5, color=color, variant=1)
    elif kind == "apa":
        line(ax, [x, x + 4.3], [y, y], color=color, lw=0.8)
        ax.add_patch(Circle((x + 4.3, y), 0.42, facecolor=color, edgecolor="none", zorder=9))
        add_text(ax, x + 1.7, y + 1.0, "AAA", size=4.2, color=color, weight="bold")
    else:
        ax.add_patch(Polygon([(x, y - 1.4), (x + 1.1, y + 1.3), (x + 2.2, y - 1.4)], closed=True, facecolor="none", edgecolor=color, linewidth=0.8, zorder=9))
        ax.add_patch(Polygon([(x + 3.0, y - 1.4), (x + 4.1, y + 1.3), (x + 5.2, y - 1.4)], closed=True, facecolor="none", edgecolor=color, linewidth=0.8, zorder=9))
        line(ax, [x + 0.2, x + 5.0], [y - 1.55, y - 1.55], color=color, lw=0.6)


def draw_evidence_panel(ax: plt.Axes) -> None:
    x, y, w, h = 88.0, 36.0, 49.0, 32.0
    rounded_box(ax, x, y, w, h, face=WHITE)
    add_text(ax, x + w / 2, y + h - 2.0, "4 aligned RNA evidence layers", size=8.0, weight="bold", va="top")
    add_text(ax, x + w / 2, y + h - 5.4, "uniform IsoQuant reprocessing", size=5.0, color=MID)

    rows = [
        ("expression", "Gene expression", "cell × gene", BLUE, PALE_BLUE, (PALE_BLUE, "#C9E0F0", "#A8CCE5")),
        ("isoform", "Isoform usage (DIU)", "cell × isoform", PURPLE, LAVENDER, (LAVENDER, "#DCD5EE", "#BFB2DD")),
        ("apa", "poly(A) sites (APA)", "cell × PAS", ORANGE, PALE_GOLD, (PALE_GOLD, "#F4E1B9", "#E7C277")),
        ("ase", "Allelic balance + junctions", "cell × evidence", TEAL, PALE_TEAL, (PALE_TEAL, "#CDE8E5", "#9BD0CA")),
    ]
    row_h = 5.55
    row_y0 = y + h - 11.1
    for i, (kind, label, matrix_label, color, fill, palette) in enumerate(rows):
        yy = row_y0 - i * row_h
        rounded_box(ax, x + 2.0, yy - 2.35, w - 4.0, 4.8, face=fill, edge="#D5D6DC", lw=0.45, radius=0.9, zorder=4)
        draw_evidence_icon(ax, x + 4.0, yy, kind, color)
        add_text(ax, x + 11.1, yy + 0.45, label, size=5.45, weight="bold", ha="left")
        add_text(ax, x + 11.1, yy - 1.18, matrix_label, size=4.55, color=MID, ha="left")
        draw_matrix(
            ax,
            x + w - 11.0,
            yy - 1.65,
            8.0,
            3.3,
            rows=3,
            cols=6,
            palette=palette,
        )


def draw_analysis_panel(ax: plt.Axes) -> None:
    x, y, w, h = 140.0, 8.0, 40.0, 60.0
    rounded_box(ax, x, y, w, h, face=WHITE)
    add_text(ax, x + w / 2, y + h - 2.1, "Query & online analysis", size=8.5, weight="bold", va="top")

    # Search box
    rounded_box(ax, x + 4.0, y + h - 11.3, w - 8.0, 5.2, face=PALE_BLUE, edge=BLUE, lw=0.65, radius=1.1, zorder=4)
    draw_search_icon(ax, x + 7.3, y + h - 8.7, color=BLUE)
    add_text(ax, x + 11.0, y + h - 8.65, "PTPRC / CD45", size=5.6, color=INK, ha="left")
    arrow(ax, (x + w / 2, y + h - 12.0), (x + w / 2, y + h - 15.0), color=PURPLE, lw=0.8, scale=7.5)

    add_text(ax, x + w / 2, y + h - 17.0, "Multi-layer gene card", size=6.1, color=PURPLE_DARK, weight="bold")
    draw_matrix(
        ax,
        x + 6.0,
        y + h - 27.5,
        w - 12.0,
        7.2,
        rows=4,
        cols=10,
        palette=(PALE_BLUE, "#C6DFEE", PALE_GOLD, "#EED9A9", PALE_TEAL),
    )
    add_text(ax, x + 2.5, y + h - 23.8, "RNA\nlayers", size=4.6, color=MID, rotation=90)
    arrow(ax, (x + w / 2, y + h - 28.2), (x + w / 2, y + h - 31.0), color=PURPLE, lw=0.8, scale=7.5)

    add_text(ax, x + w / 2, y + h - 33.0, "Cell-type DIU · APA · ASE", size=5.7, color=INK, weight="bold")
    rng = np.random.default_rng(21)
    cluster_specs = [
        (x + 10.5, y + 19.2, BLUE),
        (x + 20.3, y + 17.0, CORAL),
        (x + 29.6, y + 20.0, ORANGE),
    ]
    for cx, cy, color in cluster_specs:
        pts = rng.normal(0.0, 1.0, size=(23, 2))
        pts[:, 0] *= 1.4
        pts[:, 1] *= 1.1
        ax.scatter(pts[:, 0] + cx, pts[:, 1] + cy, s=2.8, c=color, alpha=0.82, edgecolors="none", zorder=8)
    arrow(ax, (x + w / 2, y + 13.8), (x + w / 2, y + 11.0), color=PURPLE, lw=0.8, scale=7.5)

    rounded_box(ax, x + 4.0, y + 3.0, w - 8.0, 6.2, face=LAVENDER, edge=PURPLE, lw=0.6, radius=1.1, zorder=4)
    draw_download_icon(ax, x + 8.0, y + 7.1, color=PURPLE)
    add_text(ax, x + 12.0, y + 7.1, "Export", size=5.8, weight="bold", color=PURPLE_DARK, ha="left")
    add_text(ax, x + 12.0, y + 4.6, "CSV · JSON · API", size=5.3, color=INK, ha="left")


def draw_bottom_strip(ax: plt.Axes) -> None:
    x, y, w, h = 3.0, 8.0, 134.0, 26.0
    rounded_box(ax, x, y, w, h, face=WHITE)
    modules = [
        ("Browsing", "systems · species\nstudies · samples", TEAL),
        ("Searching", "gene · isoform\njunction · region", BLUE),
        ("Downloading", "catalog · sample\nlayered evidence", ORANGE),
        ("Gene detail", "expression · DIU\nAPA · ASE", PURPLE),
        ("Precomputed maps", "cell types · studies\ncross-layer views", CORAL),
    ]
    module_w = w / len(modules)
    for i, (title, body, color) in enumerate(modules):
        mx = x + i * module_w
        if i:
            line(ax, [mx, mx], [y + 2.2, y + h - 2.2], color=LINE, lw=0.45, zorder=5)
        add_text(ax, mx + module_w / 2, y + h - 2.2, title, size=7.0, weight="bold", va="top")

        icon_y = y + 14.0
        if i == 0:
            draw_human(ax, mx + 7.0, icon_y - 1.0, 0.52)
            draw_mouse(ax, mx + 15.0, icon_y - 0.8, 0.42)
            draw_cell_cluster(ax, mx + 21.8, icon_y - 1.7, color=TEAL, scale=0.48)
        elif i == 1:
            draw_search_icon(ax, mx + 7.0, icon_y, color=BLUE)
            rounded_box(ax, mx + 11.0, icon_y - 2.3, 12.0, 4.6, face=PALE_BLUE, edge=BLUE, lw=0.55, radius=0.8, zorder=4)
            line(ax, [mx + 12.5, mx + 21.5], [icon_y + 0.8, icon_y + 0.8], color=MID, lw=0.55)
            line(ax, [mx + 12.5, mx + 19.0], [icon_y - 0.7, icon_y - 0.7], color=MID, lw=0.55)
        elif i == 2:
            draw_folder(ax, mx + 4.0, icon_y + 1.2, PURPLE)
            draw_folder(ax, mx + 10.0, icon_y - 0.3, BLUE)
            draw_folder(ax, mx + 16.0, icon_y - 1.8, ORANGE)
        elif i == 3:
            draw_matrix(
                ax,
                mx + 4.0,
                icon_y - 3.0,
                12.0,
                6.0,
                rows=4,
                cols=7,
                palette=(PALE_BLUE, "#C7E0EF", PALE_GOLD, "#EACD91"),
            )
            draw_transcript(ax, mx + 18.0, icon_y + 1.2, 6.0, color=PURPLE, variant=0)
            draw_transcript(ax, mx + 18.0, icon_y - 1.2, 6.0, color=TEAL, variant=1)
        else:
            for j, (label, color) in enumerate([("DIU", BLUE), ("APA", ORANGE), ("ASE", TEAL)]):
                yy = icon_y + 2.3 - j * 2.4
                ax.add_patch(Circle((mx + 6.0, yy), 0.75, facecolor=color, edgecolor="none", zorder=8))
                add_text(ax, mx + 8.0, yy, label, size=4.7, color=color, weight="bold", ha="left")
            draw_cell_cluster(ax, mx + 20.5, icon_y - 1.2, color=CORAL, scale=0.52)
        add_text(ax, mx + module_w / 2, y + 4.2, body, size=4.65, color=MID, va="center")


def draw_header_and_footer(ax: plt.Axes) -> None:
    add_text(
        ax,
        W_MM / 2,
        74.1,
        "scTHREAD: a unified single-cell long-read isoform database",
        size=13.5,
        weight="bold",
    )
    add_text(ax, W_MM / 2, 70.5, "https://scthread.ai4sc.ac.cn", size=6.6, color=BLUE, weight="bold")

    rounded_box(ax, 3.0, 1.4, 177.0, 4.6, face=PALE_TEAL, edge=TEAL, lw=0.65, radius=1.2, zorder=3)
    add_text(
        ax,
        W_MM / 2,
        3.7,
        "Uniform IsoQuant reprocessing of public sc long-read RNA-seq · cell-resolved isoforms · registration-free access",
        size=5.7,
        color=TEAL,
        weight="bold",
    )


def render(stem: Path, dpi: int = 450) -> list[Path]:
    setup_style()
    stem.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4), facecolor="none")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")

    # A full-canvas white object keeps the transparent export readable on dark viewers.
    ax.add_patch(Rectangle((0, 0), W_MM, H_MM, facecolor=WHITE, edgecolor="none", zorder=0))

    draw_header_and_footer(ax)
    draw_coverage_panel(ax)
    draw_database_panel(ax)
    draw_evidence_panel(ax)
    draw_bottom_strip(ax)
    draw_analysis_panel(ax)

    # Reference-like directional grammar: corpus -> database -> evidence -> tools -> analysis.
    arrow(ax, (47.6, 58.0), (49.9, 58.0), color=PURPLE, lw=1.4, scale=13.0)
    arrow(ax, (85.6, 51.7), (87.9, 51.7), color=PURPLE, lw=1.1, scale=10.0)
    arrow(ax, (124.0, 36.0), (124.0, 34.2), color=PURPLE, lw=1.4, scale=12.0)
    arrow(ax, (137.0, 21.0), (139.8, 21.0), color=PURPLE, lw=1.4, scale=12.0)
    arrow(ax, (137.0, 52.0), (139.8, 52.0), color=PURPLE, lw=1.0, scale=9.0)

    outputs = [
        stem.with_suffix(".svg"),
        stem.with_suffix(".pdf"),
        stem.with_suffix(".png"),
    ]
    fig.savefig(outputs[0], format="svg", transparent=True)
    fig.savefig(outputs[1], format="pdf", transparent=True)
    fig.savefig(outputs[2], format="png", dpi=dpi, transparent=True)
    plt.close(fig)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=DEFAULT_STEM,
        help="Output stem; .svg, .pdf and .png are written.",
    )
    parser.add_argument("--dpi", type=int, default=450, help="PNG resolution.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = render(args.stem.resolve(), dpi=args.dpi)
    for output in outputs:
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
