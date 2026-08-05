#!/usr/bin/env python3
"""Render the website-centred scTHREAD NAR Figure 3 v4.

Figure contract
---------------
Core conclusion
    A PTPRC query connects cell context, significant cell-type-specific
    differential transcript usage, junction support and exportable evidence.
Evidence hierarchy
    Cell map with ENST00000367364 expression (visual hero) ->
    isoform usage (quantitative result) ->
    junction track (long-read support) -> export/API (reproducibility).
Guardrails
    The Cell map is exploratory; the displayed five isoforms retain the
    all-23-isoform denominator; 317 filtered junctions are distinct from
    the 8,994 database-wide junctions.

The renderer is offline and deterministic. It imports the already validated
v3 data loader, then combines three real portal crops with vector summaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

import nar_style as S
import render_nar_fig3_v3 as V3


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
WALKTHROUGH = FIGURES / "website_walkthrough" / "ptprc_views_v3"
ISOFORM_WALKTHROUGH = FIGURES / "website_walkthrough" / "ptprc_views_v4_isoform"

GENE = V3.GENE
GID = V3.GID
WIDTH_MM = V3.WIDTH_MM
HEIGHT_MM = V3.HEIGHT_MM

INK = V3.INK
SLATE = V3.SLATE
TEAL = V3.TEAL
BLUE = V3.BLUE
GOLD = V3.GOLD
GREY = V3.GREY
PALE = V3.PALE
PALE_BLUE = V3.PALE_BLUE
WHITE = V3.WHITE
GRID = V3.GRID

CELL_MAP_IMAGE = ISOFORM_WALKTHROUGH / "02_cell_map.png"
JUNCTION_IMAGE = WALKTHROUGH / "01_junctions.png"


def _panel_heading(
    fig: plt.Figure,
    letter: str,
    title: str,
    *,
    x: float,
    y: float,
    title_size: float = 7.0,
) -> None:
    """Draw an upright Nature-style panel label and compact claim heading."""
    fig.text(
        x - 0.030,
        y,
        letter,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        fontstyle="normal",
        fontfamily=S._FAM,
        color=INK,
    )
    fig.text(
        x,
        y,
        title,
        ha="left",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        fontstyle="normal",
        fontfamily=S._FAM,
        color=INK,
    )


def _rounded(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    *,
    facecolor: str,
    edgecolor: str = "none",
    linewidth: float = 0.55,
    radius: float = 0.025,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def _image_axes(
    fig: plt.Figure,
    rect: list[float],
    path: Path,
    *,
    y_slice: tuple[int, int],
    x_slice: tuple[int, int] | None = None,
) -> plt.Axes:
    if not path.is_file():
        raise FileNotFoundError(path)
    image = mpimg.imread(path)
    x0, x1 = x_slice if x_slice is not None else (0, image.shape[1])
    y0, y1 = y_slice
    if not (0 <= x0 < x1 <= image.shape[1] and 0 <= y0 < y1 <= image.shape[0]):
        raise ValueError(
            f"Invalid crop for {path.name}: x={x0}:{x1}, y={y0}:{y1}, "
            f"image={image.shape[1]}×{image.shape[0]}"
        )
    ax = fig.add_axes(rect)
    ax.imshow(
        image[y0:y1, x0:x1],
        interpolation="lanczos",
        aspect="auto",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#AEB8BF")
        spine.set_linewidth(0.55)
    return ax


def _status_chip(
    ax: plt.Axes,
    *,
    x: float,
    width: float,
    label: str,
    value: str,
    accent: str,
    pale: str,
) -> None:
    _rounded(
        ax,
        (x, 0.08),
        width,
        0.84,
        facecolor=pale,
        edgecolor=GRID,
        linewidth=0.5,
        radius=0.022,
    )
    ax.plot(
        [x + 0.013, x + 0.013],
        [0.22, 0.79],
        transform=ax.transAxes,
        color=accent,
        linewidth=2.2,
        solid_capstyle="round",
        clip_on=False,
    )
    ax.text(
        x + 0.028,
        0.69,
        label,
        transform=ax.transAxes,
        fontsize=4.3,
        fontweight="bold",
        color=accent,
        ha="left",
        va="center",
    )
    ax.text(
        x + 0.028,
        0.38,
        value,
        transform=ax.transAxes,
        fontsize=4.65,
        fontweight="bold",
        color=INK,
        ha="left",
        va="center",
    )


def _draw_query_panel(fig: plt.Figure, inputs: V3.Inputs) -> None:
    _panel_heading(
        fig,
        "a",
        "A live PTPRC query separates the current view from database-wide evidence",
        x=0.055,
        y=0.978,
    )
    fig.text(
        0.985,
        0.978,
        GID,
        ha="right",
        va="top",
        fontsize=4.65,
        color=SLATE,
        fontfamily=S._FAM,
    )
    # Search controls and genomic locus from the real live gene page. The
    # quantitative cards are rebuilt below as vector text for legibility.
    _image_axes(
        fig,
        [0.055, 0.824, 0.930, 0.119],
        inputs.portal_image,
        y_slice=(0, 280),
    )

    ax = fig.add_axes([0.055, 0.741, 0.930, 0.065])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _status_chip(
        ax,
        x=0.000,
        width=0.150,
        label="CURRENT VIEW",
        value="317 junctions · ≥10 molecules",
        accent=BLUE,
        pale="#EDF3F7",
    )
    _status_chip(
        ax,
        x=0.160,
        width=0.238,
        label="DATABASE-WIDE",
        value=(
            f"{int(inputs.api['coverage']['junctions']['features']):,} junctions · "
            f"{int(inputs.api['coverage']['junctions']['molecules']):,} molecules"
        ),
        accent=SLATE,
        pale="#F0F2F4",
    )
    _status_chip(
        ax,
        x=0.408,
        width=0.181,
        label="CELL-TYPE DTU · SIG.",
        value=f"q = {float(inputs.diu['qval']):.5f} · effect = {float(inputs.diu['effect']):.3f}",
        accent=TEAL,
        pale="#ECF4F2",
    )
    _status_chip(
        ax,
        x=0.599,
        width=0.181,
        label="APA · SIGNIFICANT",
        value=f"q = {float(inputs.apa['qval']):.5f} · effect = {float(inputs.apa['effect']):.3f}",
        accent=GOLD,
        pale="#F7F1E7",
    )
    _status_chip(
        ax,
        x=0.790,
        width=0.210,
        label="ASE · NOT SIGNIFICANT",
        value=f"q = {float(inputs.ase['qval']):.2f} · effect = {float(inputs.ase['effect']):.3f}",
        accent=GREY,
        pale="#F1F2F3",
    )


def _draw_cell_map_panel(fig: plt.Figure) -> None:
    _panel_heading(
        fig,
        "b",
        "Cell map places the leading PTPRC isoform in cell-type context",
        x=0.055,
        y=0.710,
        title_size=6.9,
    )
    _image_axes(
        fig,
        [0.055, 0.383, 0.555, 0.292],
        CELL_MAP_IMAGE,
        y_slice=(0, 1120),
    )
    fig.text(
        0.055,
        0.365,
        "71,913 sampled cells · ENST00000367364 expression (0–5 long-read molecules) "
        "· exploratory UMAP · no trajectory inference",
        fontsize=4.8,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_usage_panel(fig: plt.Figure, inputs: V3.Inputs) -> None:
    _panel_heading(
        fig,
        "c",
        "PTPRC exhibits cell-type-specific DTU",
        x=0.655,
        y=0.710,
        title_size=6.75,
    )
    matrix, totals, top_coverage = V3._usage_matrix(inputs.usage)
    fig.text(
        0.655,
        0.677,
        f"CELL-TYPE DTU · q = {float(inputs.diu['qval']):.5f} "
        f"· effect = {float(inputs.diu['effect']):.3f} · 23 isoforms "
        f"· Σ = {int(totals.sum()):,} molecules",
        fontsize=4.45,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax = fig.add_axes([0.710, 0.424, 0.175, 0.215])
    mesh = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1),
        np.arange(matrix.shape[0] + 1),
        matrix.to_numpy(),
        cmap="viridis_r",
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        linewidth=0.32,
        edgecolors=WHITE,
        rasterized=False,
    )
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    labels = ["…" + transcript[-6:] for transcript in matrix.columns]
    ax.set_xticks(np.arange(len(labels)) + 0.5)
    ax.set_xticklabels(labels, rotation=31, ha="left", rotation_mode="anchor")
    ax.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        length=0,
        pad=1.2,
        labelsize=4.65,
    )
    for index, tick in enumerate(ax.get_xticklabels()):
        if index == 0:
            tick.set_color(TEAL)
            tick.set_fontweight("bold")
        elif index == 1:
            tick.set_color(BLUE)
            tick.set_fontweight("bold")
    ax.set_yticks(np.arange(len(V3.CELL_TYPE_ORDER)) + 0.5)
    ax.set_yticklabels(V3.CELL_TYPE_ORDER, fontsize=5.0)
    ax.tick_params(axis="y", length=0, pad=1.7)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = float(matrix.iat[row, col])
            ax.text(
                col + 0.5,
                row + 0.5,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=4.25,
                fontweight="bold" if value >= 0.35 else "normal",
                color=WHITE if value >= 0.32 else INK,
            )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#BBC3C8")
        spine.set_linewidth(0.5)

    bar = fig.add_axes([0.908, 0.424, 0.052, 0.215])
    y = np.arange(len(V3.CELL_TYPE_ORDER))
    bar.barh(
        y,
        totals.to_numpy(),
        left=1.0,
        color=PALE_BLUE,
        edgecolor=BLUE,
        linewidth=0.42,
        height=0.58,
    )
    bar.set_xscale("log")
    bar.set_xlim(100, 100_000)
    bar.set_ylim(len(V3.CELL_TYPE_ORDER) - 0.5, -0.5)
    bar.set_xticks([])
    bar.set_yticks([])
    for spine in bar.spines.values():
        spine.set_visible(False)
    bar.text(
        0.0,
        1.13,
        "All 23",
        transform=bar.transAxes,
        fontsize=4.4,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
    )
    bar.text(
        0.0,
        1.055,
        "molecules",
        transform=bar.transAxes,
        fontsize=4.1,
        color=SLATE,
        ha="left",
        va="bottom",
    )
    for row, value in enumerate(totals.astype(int)):
        bar.text(
            1.04,
            row,
            f"{value:,}",
            transform=bar.get_yaxis_transform(),
            fontsize=4.15,
            color=INK,
            ha="left",
            va="center",
            clip_on=False,
        )

    cax = fig.add_axes([0.710, 0.395, 0.175, 0.0065])
    edges = np.linspace(0.0, 0.65, 66)
    values = (edges[:-1] + edges[1:]) / 2.0
    cax.pcolormesh(
        edges,
        [0.0, 1.0],
        values[np.newaxis, :],
        cmap=mesh.cmap,
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        rasterized=False,
    )
    cax.set_xlim(0.0, 0.65)
    cax.set_ylim(0.0, 1.0)
    cax.set_yticks([])
    cax.set_xticks([0.0, 0.3, 0.6])
    cax.tick_params(axis="x", labelsize=4.1, length=1.5, width=0.4, pad=0.6)
    for spine in cax.spines.values():
        spine.set_linewidth(0.4)
        spine.set_color(INK)
    fig.text(
        0.7975,
        0.376,
        "Usage fraction",
        fontsize=4.55,
        color=SLATE,
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.655,
        0.366,
        "DTU contrast: …367364 = 0.59 in Monocyte/Plasma; "
        "…697630 = 0.40 in Progenitor.",
        fontsize=4.45,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.655,
        0.348,
        f"Top five cover {top_coverage.min() * 100:.0f}–"
        f"{top_coverage.max() * 100:.0f}% per row; values not renormalized "
        "(all-23-isoform denominator).",
        fontsize=4.55,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_junction_access_panel(fig: plt.Figure, inputs: V3.Inputs) -> None:
    _panel_heading(
        fig,
        "d",
        "Junction evidence remains inspectable and exportable",
        x=0.055,
        y=0.328,
        title_size=6.9,
    )
    _image_axes(
        fig,
        [0.055, 0.105, 0.625, 0.190],
        JUNCTION_IMAGE,
        y_slice=(95, 780),
    )
    fig.text(
        0.055,
        0.081,
        "317 filtered junctions; the portal table renders its first 100 rows. "
        "Selecting an arc exposes molecules, reads, runs and studies.",
        fontsize=4.65,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    ax = fig.add_axes([0.715, 0.105, 0.270, 0.190])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _rounded(
        ax,
        (0.0, 0.0),
        1.0,
        1.0,
        facecolor=INK,
        edgecolor=INK,
        linewidth=0,
        radius=0.035,
    )
    _rounded(
        ax,
        (0.045, 0.785),
        0.280,
        0.125,
        facecolor="#2D5D58",
        edgecolor="#6EA49D",
        linewidth=0.5,
        radius=0.025,
    )
    ax.text(
        0.185,
        0.846,
        "EXPORT CSV",
        fontsize=4.75,
        fontweight="bold",
        color=WHITE,
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.365,
        0.872,
        "GET /api/gene/PTPRC/overview",
        fontsize=4.05,
        fontweight="bold",
        color="#B7D8D3",
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.365,
        0.818,
        "?species=human",
        fontsize=3.95,
        color="#B7D8D3",
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.plot([0.045, 0.955], [0.730, 0.730], color="#465063", linewidth=0.5)

    coverage = inputs.api["coverage"]
    rows = [
        (
            "ISOFORM EVIDENCE",
            f"{int(coverage['isoforms']['features']):,} isoforms",
            f"{int(coverage['isoforms']['molecules']):,} molecules",
        ),
        (
            "POLY(A) EVIDENCE",
            f"{int(coverage['pas']['features']):,} sites",
            f"{int(coverage['pas']['molecules']):,} molecules",
        ),
        (
            "JUNCTION EVIDENCE",
            f"{int(coverage['junctions']['features']):,} junctions",
            f"{int(coverage['junctions']['molecules']):,} molecules",
        ),
    ]
    for y, (label, features, molecules) in zip([0.635, 0.485, 0.335], rows):
        ax.text(
            0.055,
            y,
            label,
            fontsize=3.95,
            fontweight="bold",
            color="#AAB2C0",
            ha="left",
            va="center",
            fontfamily=S._FAM,
        )
        ax.text(
            0.055,
            y - 0.057,
            features,
            fontsize=4.45,
            color=WHITE,
            ha="left",
            va="center",
            fontfamily=S._FAM,
        )
        ax.text(
            0.945,
            y - 0.057,
            molecules,
            fontsize=4.45,
            color=WHITE,
            ha="right",
            va="center",
            fontfamily=S._FAM,
        )
    _rounded(
        ax,
        (0.045, 0.075),
        0.910,
        0.115,
        facecolor="#263B43",
        edgecolor="#55767B",
        linewidth=0.5,
        radius=0.020,
    )
    ax.text(
        0.500,
        0.132,
        "tables = snapshot = live API",
        fontsize=4.45,
        fontweight="bold",
        color="#8FD1C6",
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )


def render(inputs: V3.Inputs, stem: Path, dpi: int) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(
        {
            "font.family": S._FAM,
            "font.size": 6.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
            "axes.unicode_minus": False,
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
        }
    ):
        fig = plt.figure(figsize=(WIDTH_MM * S.MM, HEIGHT_MM * S.MM), facecolor="none")
        _draw_query_panel(fig, inputs)
        _draw_cell_map_panel(fig)
        _draw_usage_panel(fig, inputs)
        _draw_junction_access_panel(fig, inputs)

        outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
        for output in outputs:
            fig.savefig(
                output,
                dpi=dpi,
                transparent=True,
                bbox_inches=None,
                pad_inches=0,
            )
        plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v4",
        help="Output path without an extension.",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all frozen inputs and portal screenshots without rendering.",
    )
    args = parser.parse_args()

    for image in (CELL_MAP_IMAGE, JUNCTION_IMAGE):
        if not image.is_file():
            raise FileNotFoundError(image)
    inputs = V3.load_inputs()
    print(
        "INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "cell_types": inputs.usage["ct"].nunique(),
            "isoforms": inputs.usage["transcript_id"].nunique(),
            "matrix_molecules": int(round(inputs.usage["count"].sum())),
            "portal_images": 3,
        },
    )
    if args.validate_only:
        return
    outputs = render(inputs, args.stem.resolve(), args.dpi)
    for output in outputs:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
