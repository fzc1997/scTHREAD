#!/usr/bin/env python3
"""Composite the GPT Image 2 skeleton with exact text and real vector plots."""

from __future__ import annotations

import os

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
COMPONENTS = ROOT / "figures/ga_gpt2_components_v2"
DEFAULT_BACKGROUND = (
    COMPONENTS
    / "mcp_generation"
    / "gpt_image_20260726_163535_6a65b268_2.png"
)
DEFAULT_STEM = ROOT / "figures/scTHREAD_graphical_abstract_gpt2_hybrid_v2"
ARIAL = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))

W_MM = 183.0
H_MM = 78.0

PURPLE = "#5B4B9A"
PURPLE_DARK = "#443677"
TEAL = "#0B7C86"
BLUE = "#2F6FB2"
ORANGE = "#D88B1B"
CORAL = "#E34A33"
GREEN = "#4C956C"
INK = "#202124"
MID = "#6D6E73"
LIGHT = "#D7D9DF"
WHITE = "#FFFFFF"
PALE_TEAL = "#EEF8F8"
PALE_PEACH = "#FFF1E7"
PALE_LAV = "#F4F1FA"

SYSTEM_COLORS = {
    "Blood/immune": TEAL,
    "Neural/sensory": PURPLE,
    "Cancer": "#D55E00",
    "Endocrine": BLUE,
    "Heart/vascular": "#7A9E3A",
    "Development/embryo": ORANGE,
    "Reproductive": "#CC78BC",
    "Other tissues": "#8B7355",
}

SYSTEM_SHORT = {
    "Blood/immune": "Blood",
    "Neural/sensory": "Neural",
    "Cancer": "Cancer",
    "Endocrine": "Pancreas",
    "Heart/vascular": "Heart",
    "Development/embryo": "Embryo",
    "Reproductive": "Ovary",
    "Other tissues": "Other",
}

SYSTEM_LEGEND_ORDER = [
    "Blood/immune",
    "Neural/sensory",
    "Cancer",
    "Endocrine",
    "Heart/vascular",
    "Development/embryo",
    "Reproductive",
    "Other tissues",
]

UMAP_COLORS = {
    "Lymphoid": BLUE,
    "Myeloid": TEAL,
    "Erythroid": CORAL,
    "Progenitor": ORANGE,
    "Neural": PURPLE,
    "Cardiovascular": "#7A9E3A",
    "Stromal": "#B06A3C",
    "Other": "#A0A4AB",
}


def setup_style() -> None:
    if ARIAL.exists():
        font_manager.fontManager.addfont(str(ARIAL))
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 450,
            "figure.dpi": 150,
            "savefig.transparent": True,
            "axes.grid": False,
        }
    )


def text(
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
    zorder: int = 30,
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
        family="Arial",
        clip_on=False,
        zorder=zorder,
    )


def rounded(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face: str = WHITE,
    edge: str = LIGHT,
    lw: float = 0.55,
    radius: float = 0.8,
    alpha: float = 1.0,
    zorder: int = 15,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle=f"round,pad=0.10,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            alpha=alpha,
            zorder=zorder,
        )
    )


def add_axes_mm(
    fig: plt.Figure, x: float, y: float, width: float, height: float
) -> plt.Axes:
    return fig.add_axes([x / W_MM, y / H_MM, width / W_MM, height / H_MM])


def draw_transcript(ax: plt.Axes, x: float, y: float, width: float) -> None:
    exon_w = width * 0.15
    gap = (width - 3 * exon_w) / 2
    xs = [x, x + exon_w + gap, x + 2 * (exon_w + gap)]
    ax.plot(
        [xs[0] + exon_w, xs[1], xs[1] + exon_w, xs[2]],
        [y, y, y, y],
        color=TEAL,
        lw=1.0,
        zorder=24,
    )
    for idx, xx in enumerate(xs):
        ax.add_patch(
            Rectangle(
                (xx, y - 0.55),
                exon_w,
                1.1,
                facecolor=WHITE if idx == 1 else PALE_TEAL,
                edgecolor=TEAL,
                linewidth=0.8,
                zorder=25,
            )
        )


def overlay_catalog_strip(ax: plt.Axes) -> None:
    source = pd.read_csv(
        COMPONENTS / "catalog_system_composition_20260726.tsv", sep="\t"
    )
    source = source.sort_values("cells", ascending=False).reset_index(drop=True)
    rounded(ax, 1.7, 21.6, 43.6, 8.6, face=WHITE, edge="#AEB1B7", lw=0.55, radius=0.5)
    text(
        ax,
        23.5,
        29.1,
        "469 samples · 31 studies · Human + mouse · ONT + PacBio",
        size=3.6,
        weight="bold",
        color=PURPLE_DARK,
    )
    x0, y0, total_w, bar_h = 3.0, 26.1, 40.8, 1.45
    total = float(source["cells"].sum())
    cursor = x0
    centers: list[tuple[float, str, str]] = []
    for _, row in source.iterrows():
        width = total_w * float(row["cells"]) / total
        color = SYSTEM_COLORS.get(str(row["system"]), MID)
        ax.add_patch(
            Rectangle(
                (cursor, y0),
                width,
                bar_h,
                facecolor=color,
                edgecolor=WHITE,
                linewidth=0.25,
                zorder=22,
            )
        )
        centers.append(
            (
                cursor + width / 2,
                SYSTEM_SHORT.get(str(row["system"]), str(row["system"])),
                color,
            )
        )
        cursor += width
    # Two-row legend keeps the small but biologically important ovary and
    # embryo categories visible instead of reverting to method-based labels.
    available = {
        str(row["system"]): row for _, row in source.iterrows()
    }
    legend_items = [
        system for system in SYSTEM_LEGEND_ORDER if system in available
    ]
    for idx, system in enumerate(legend_items):
        row_idx, col_idx = divmod(idx, 4)
        xpos = [4.5, 14.4, 24.6, 34.8][col_idx]
        ypos = [24.15, 22.35][row_idx]
        label = SYSTEM_SHORT[system]
        color = SYSTEM_COLORS[system]
        ax.scatter([xpos - 0.75], [ypos], s=4.2, color=color, zorder=23)
        text(ax, xpos, ypos, label, size=2.55, color=INK, ha="left")


def overlay_featured_biology(ax: plt.Axes) -> None:
    rounded(
        ax,
        34.5,
        66.0,
        18.3,
        2.8,
        face="#FFF2F8",
        edge="#CC78BC",
        lw=0.45,
        radius=0.55,
        alpha=0.94,
        zorder=20,
    )
    text(
        ax,
        43.65,
        67.4,
        "Human ovary · scONT",
        size=3.8,
        color="#8F4B86",
        weight="bold",
    )
    rounded(
        ax,
        19.4,
        43.2,
        20.3,
        2.8,
        face="#FFF8E9",
        edge=ORANGE,
        lw=0.45,
        radius=0.55,
        alpha=0.94,
        zorder=20,
    )
    text(
        ax,
        29.55,
        44.6,
        "Mouse gastrula · scONT",
        size=3.65,
        color="#9A6511",
        weight="bold",
    )


def overlay_database(ax: plt.Axes) -> None:
    text(ax, 80.2, 68.9, "Database", size=6.9, weight="bold")
    text(ax, 80.2, 56.2, "scTHREAD", size=12.3, color=PURPLE, weight="bold")
    draw_transcript(ax, 73.7, 52.2, 13.0)
    text(ax, 80.2, 48.5, "845,781", size=9.6, color=PURPLE_DARK, weight="bold")
    text(ax, 80.2, 45.6, "cells", size=5.8, color=PURPLE_DARK, weight="bold")
    text(ax, 80.2, 41.5, ">200k isoforms", size=6.2, color=TEAL, weight="bold")
    rounded(
        ax,
        68.0,
        32.5,
        24.4,
        3.0,
        face=WHITE,
        edge="none",
        lw=0,
        radius=0.45,
        alpha=0.90,
        zorder=20,
    )
    text(
        ax,
        80.2,
        34.0,
        "uniform IsoQuant reprocessing",
        size=4.6,
        color=MID,
        weight="bold",
    )


def overlay_evidence(ax: plt.Axes) -> None:
    rounded(
        ax,
        101.0,
        68.4,
        36.2,
        3.5,
        face=WHITE,
        edge="none",
        lw=0,
        radius=0.4,
        alpha=0.94,
        zorder=18,
    )
    text(
        ax,
        119.1,
        70.1,
        "4 aligned RNA evidence layers",
        size=7.1,
        weight="bold",
    )

    rows = [
        (63.3, "Gene expression", PALE_TEAL),
        (51.8, "Isoform usage (DIU)", PALE_LAV),
        (40.1, "poly(A) sites (APA)", "#FFF8E9"),
        (28.2, "Allelic balance + junctions (ASE)", "#FFF2EC"),
    ]
    for y, label, fill in rows:
        rounded(
            ax,
            117.7,
            y - 1.8,
            18.9,
            3.6,
            face=fill,
            edge="none",
            lw=0,
            radius=0.45,
            alpha=0.94,
            zorder=19,
        )
        text(
            ax,
            127.15,
            y,
            label,
            size=4.25 if "Allelic" not in label else 3.55,
            weight="bold",
            color=INK,
        )
    rounded(
        ax,
        102.7,
        25.7,
        10.6,
        5.0,
        face="#FFF2EC",
        edge="none",
        lw=0,
        radius=0.4,
        alpha=0.98,
        zorder=20,
    )
    text(ax, 108.0, 29.2, "REF", size=3.6, color=PURPLE, weight="bold")
    text(ax, 108.0, 27.3, "ALT", size=3.6, color=GREEN, weight="bold")


def draw_ptprc_heatmap(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(COMPONENTS / "ptprc_two_isoform_switch.tsv", sep="\t")
    transcripts = ["ENST00000367364", "ENST00000697630"]
    cell_order = [
        "B cell",
        "Dendritic cell",
        "Erythroid",
        "Monocyte",
        "NK",
        "Plasma cell",
        "Progenitor",
        "T cell",
    ]
    matrix = (
        source.pivot(index="transcript_id", columns="cell_type", values="fraction")
        .reindex(index=transcripts, columns=cell_order)
        .fillna(0)
    )
    rounded(main_ax, 68.9, 4.8, 24.0, 11.2, face=WHITE, edge="#AEB1B7", lw=0.45, radius=0.35)
    ax = add_axes_mm(fig, 69.8, 6.0, 22.1, 8.4)
    cmap = LinearSegmentedColormap.from_list(
        "ptprc", ["#F7FAFC", "#BFDCEB", "#4E91BC", "#175B8D"]
    )
    ax.pcolormesh(
        np.arange(9),
        np.arange(3),
        matrix.to_numpy(),
        cmap=cmap,
        vmin=0,
        vmax=0.60,
        shading="flat",
        edgecolors=WHITE,
        linewidth=0.35,
    )
    ax.set_xlim(0, 8)
    ax.set_ylim(2, 0)
    ax.set_aspect("auto")
    ax.set_xticks(np.arange(8) + 0.5)
    ax.set_xticklabels(["B", "DC", "Ery", "Mono", "NK", "Pl", "Prog", "T"], fontsize=3.5)
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["…67364", "…97630"], fontsize=3.4)
    ax.tick_params(length=0, pad=0.8)
    for i in range(2):
        for j in range(8):
            value = matrix.iat[i, j]
            ax.text(
                j + 0.5,
                i + 0.5,
                f"{value * 100:.0f}",
                ha="center",
                va="center",
                fontsize=3.0,
                color=WHITE if value > 0.34 else INK,
                fontweight="bold" if value > 0.30 else "normal",
            )
    for spine in ax.spines.values():
        spine.set_color(WHITE)
        spine.set_linewidth(0.4)


def draw_inventory(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(COMPONENTS / "three_axis_inventory_source.tsv", sep="\t")
    order = ["DIU", "APA", "ASE"]
    source = source.set_index("axis").reindex(order).reset_index()
    rounded(main_ax, 99.2, 4.8, 38.2, 11.2, face=WHITE, edge="#AEB1B7", lw=0.45, radius=0.4)
    ax = add_axes_mm(fig, 100.8, 6.0, 34.8, 8.5)
    colors = {"DIU": BLUE, "APA": ORANGE, "ASE": TEAL}
    y = np.arange(3)
    fraction = source["sig_fraction"].to_numpy() * 100
    ax.barh(
        y,
        fraction,
        color=[colors[a] for a in source["axis"]],
        height=0.54,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=4.0, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 45)
    ax.set_xticks([0, 20, 40])
    ax.set_xticklabels(["0", "20", "40%"], fontsize=3.3)
    ax.tick_params(length=0, pad=0.8)
    for yi, row in source.iterrows():
        ax.text(
            max(float(row["sig_fraction"]) * 100 + 1.0, 1.0),
            yi,
            f"{int(row['n_sig_fdr05_effect_gate']):,}/{int(row['n_genes_tested']):,}",
            va="center",
            ha="left",
            fontsize=3.5,
            color=INK,
        )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(LIGHT)
    ax.spines["bottom"].set_linewidth(0.35)


def draw_umap(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(COMPONENTS / "atlas_umap_stratified_sample.tsv", sep="\t")
    rounded(main_ax, 144.0, 15.7, 35.5, 40.0, face=WHITE, edge="#AEB1B7", lw=0.55, radius=0.45)
    ax = add_axes_mm(fig, 145.2, 18.2, 32.8, 33.5)
    for label, color in UMAP_COLORS.items():
        sub = source[source["broad_class"] == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap1"],
            sub["umap2"],
            s=1.9,
            c=color,
            alpha=0.72,
            linewidths=0,
            rasterized=False,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    text(
        main_ax,
        161.5,
        53.8,
        "Harmonized cell-type atlas",
        size=5.5,
        weight="bold",
        color=INK,
    )
    text(
        main_ax,
        161.5,
        16.9,
        "74,906 portal cells · exploratory UMAP",
        size=3.8,
        color=MID,
    )


def overlay_access_labels(ax: plt.Axes) -> None:
    modules = [
        (13.0, "Browsing"),
        (34.3, "Searching"),
        (55.4, "Downloading"),
        (80.8, "Gene detail · PTPRC/CD45"),
        (118.3, "Precomputed maps"),
    ]
    widths = [18.0, 18.0, 18.0, 25.0, 24.0]
    for (x, heading), width in zip(modules, widths):
        rounded(
            ax,
            x - width / 2,
            17.2,
            width,
            3.0,
            face=WHITE,
            edge="none",
            lw=0,
            radius=0.35,
            alpha=0.92,
            zorder=20,
        )
        text(ax, x, 18.7, heading, size=5.0 if x < 70 else 4.4, weight="bold")


def overlay_query_labels(ax: plt.Axes) -> None:
    rounded(ax, 143.2, 67.4, 36.2, 3.2, face=WHITE, edge="none", lw=0, radius=0.4, alpha=0.94)
    text(ax, 161.3, 69.0, "Query & online analysis", size=6.8, weight="bold")
    rounded(ax, 147.0, 61.3, 28.0, 3.0, face="#F1F7FB", edge=BLUE, lw=0.6, radius=0.8)
    text(ax, 161.0, 62.8, "PTPRC / CD45", size=5.2, color=INK, weight="bold")
    text(ax, 161.3, 57.0, "Multi-layer gene card", size=4.5, color=PURPLE_DARK, weight="bold")
    rounded(ax, 148.8, 7.2, 25.0, 3.8, face=PALE_LAV, edge=PURPLE, lw=0.55, radius=0.75)
    text(ax, 161.3, 9.1, "Export · CSV · JSON · API", size=4.6, color=PURPLE_DARK, weight="bold")


def overlay_header_footer(ax: plt.Axes) -> None:
    ax.add_patch(Rectangle((0, 72.6), W_MM, 5.4, facecolor=WHITE, edgecolor="none", zorder=16))
    text(
        ax,
        W_MM / 2,
        76.2,
        "scTHREAD: a unified single-cell long-read isoform database",
        size=13.0,
        weight="bold",
    )
    text(
        ax,
        W_MM / 2,
        73.35,
        "https://scthread.ai4sc.ac.cn",
        size=6.3,
        color=BLUE,
        weight="bold",
    )
    rounded(ax, 3.0, 0.55, 135.0, 3.7, face=PALE_TEAL, edge=TEAL, lw=0.5, radius=0.7)
    text(
        ax,
        70.5,
        2.65,
        "Uniform reprocessing · cell-resolved isoforms · open access",
        size=4.5,
        color=TEAL,
        weight="bold",
    )


def validate_sources() -> None:
    required = [
        COMPONENTS / "catalog_system_composition_20260726.tsv",
        COMPONENTS / "catalog_study_biology_classification_20260726.tsv",
        COMPONENTS / "featured_project_cohorts_20260726.tsv",
        COMPONENTS / "ptprc_two_isoform_switch.tsv",
        COMPONENTS / "atlas_umap_stratified_sample.tsv",
        COMPONENTS / "three_axis_inventory_source.tsv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing real-data components: " + ", ".join(missing))
    systems = pd.read_csv(required[0], sep="\t")
    if int(systems["cells"].sum()) != 845_781:
        raise RuntimeError("Catalog composition does not sum to 845,781 cells.")
    expected_systems = set(SYSTEM_LEGEND_ORDER)
    observed_systems = set(systems["system"].astype(str))
    if observed_systems != expected_systems:
        raise RuntimeError(
            "Biological-system categories differ from the v2 contract: "
            f"{sorted(observed_systems)}"
        )
    forbidden = {"Benchmark", "Method", "Smart-seq2", "Differentiation"}
    if observed_systems & forbidden:
        raise RuntimeError("Technical categories leaked into the biological legend.")
    featured = pd.read_csv(required[2], sep="\t", dtype=str)
    if set(featured["display_label"]) != {"Human ovary", "Mouse gastrula"}:
        raise RuntimeError("Featured ovary/gastrula cohort contract is incomplete.")


def render(background: Path, stem: Path, dpi: int = 450) -> list[Path]:
    setup_style()
    validate_sources()
    if not background.exists():
        raise FileNotFoundError(background)
    stem.parent.mkdir(parents=True, exist_ok=True)

    bg = mpimg.imread(background)
    if bg.shape[1] / bg.shape[0] < 2.30 or bg.shape[1] / bg.shape[0] > 2.38:
        raise RuntimeError(f"Unexpected GPT skeleton aspect ratio: {bg.shape}")

    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4), facecolor="none")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.imshow(bg, extent=(0, W_MM, 0, H_MM), origin="upper", aspect="auto", zorder=0)

    overlay_header_footer(ax)
    overlay_catalog_strip(ax)
    overlay_featured_biology(ax)
    overlay_database(ax)
    overlay_evidence(ax)
    overlay_access_labels(ax)
    overlay_query_labels(ax)
    draw_ptprc_heatmap(fig, ax)
    draw_inventory(fig, ax)
    draw_umap(fig, ax)

    outputs = [stem.with_suffix(".svg"), stem.with_suffix(".pdf"), stem.with_suffix(".png")]
    fig.savefig(outputs[0], format="svg", transparent=True)
    fig.savefig(outputs[1], format="pdf", transparent=True)
    fig.savefig(outputs[2], format="png", dpi=dpi, transparent=True)
    plt.close(fig)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", type=Path, default=DEFAULT_BACKGROUND)
    parser.add_argument("--stem", type=Path, default=DEFAULT_STEM)
    parser.add_argument("--dpi", type=int, default=450)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = render(args.background.resolve(), args.stem.resolve(), dpi=args.dpi)
    for output in outputs:
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
