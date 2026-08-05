#!/usr/bin/env python3
"""Render the biologically focused, real-data overlay for graphical abstract v4."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

import render_ga_gpt2_hybrid as base


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STEM = (
    ROOT / "figures/ga_gpt2_vector_v4/scTHREAD_ga_gpt2_vector_v4_overlay"
)

SYSTEM_LABELS = {
    "Blood/immune": "Blood/immune",
    "Neural/sensory": "Neural",
    "Cancer": "Cancer",
    "Endocrine": "Endocrine",
    "Heart/vascular": "Heart",
    "Development/embryo": "Embryo",
    "Reproductive": "Reproductive",
    "Other tissues": "Other",
}

UMAP_LABELS = {
    "Lymphoid": "Lymphoid",
    "Myeloid": "Myeloid",
    "Erythroid": "Erythroid",
    "Progenitor": "Progenitor",
    "Neural": "Neural",
    "Cardiovascular": "Cardiovasc.",
    "Stromal": "Stromal",
    "Other": "Other",
}


def overlay_header_footer(ax: plt.Axes) -> None:
    base.text(
        ax,
        base.W_MM / 2,
        75.85,
        "scTHREAD: a single-cell long-read transcriptome atlas",
        size=9.6,
        weight="bold",
    )
    base.text(
        ax,
        base.W_MM / 2,
        73.55,
        "https://scthread.ai4sc.ac.cn",
        size=5.7,
        color=base.BLUE,
        weight="bold",
    )
    base.rounded(
        ax,
        3.0,
        0.55,
        135.0,
        3.65,
        face=base.PALE_TEAL,
        edge=base.TEAL,
        lw=0.5,
        radius=0.65,
    )
    base.text(
        ax,
        70.5,
        2.48,
        "Human and mouse tissues · cell types · isoforms · poly(A) · allelic regulation",
        size=4.45,
        color=base.TEAL,
        weight="bold",
    )


def overlay_catalog_strip(ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "catalog_system_composition_20260726.tsv", sep="\t"
    )
    source = source.sort_values("cells", ascending=False).reset_index(drop=True)
    base.rounded(
        ax,
        1.55,
        24.5,
        41.8,
        8.9,
        face=base.WHITE,
        edge="#AEB1B7",
        lw=0.55,
        radius=0.45,
    )
    base.text(
        ax,
        22.45,
        32.35,
        "Biological breadth · 469 samples · 31 studies",
        size=3.9,
        weight="bold",
        color=base.PURPLE_DARK,
    )
    base.text(
        ax,
        22.45,
        31.12,
        "Human + mouse · ONT + PacBio",
        size=3.15,
        color=base.MID,
        weight="bold",
    )
    x0, y0, total_w, bar_h = 3.0, 29.35, 39.0, 1.15
    total = float(source["cells"].sum())
    cursor = x0
    for _, row in source.iterrows():
        width = total_w * float(row["cells"]) / total
        system = str(row["system"])
        ax.add_patch(
            Rectangle(
                (cursor, y0),
                width,
                bar_h,
                facecolor=base.SYSTEM_COLORS[system],
                edgecolor=base.WHITE,
                linewidth=0.22,
                zorder=22,
            )
        )
        cursor += width

    available = set(source["system"].astype(str))
    legend_items = [
        system for system in base.SYSTEM_LEGEND_ORDER if system in available
    ]
    x_positions = [3.55, 13.8, 23.25, 32.55]
    for idx, system in enumerate(legend_items):
        row_idx, col_idx = divmod(idx, 4)
        x = x_positions[col_idx]
        y = [27.35, 25.62][row_idx]
        ax.scatter(
            [x - 0.58],
            [y],
            s=5.8,
            color=base.SYSTEM_COLORS[system],
            zorder=23,
        )
        base.text(
            ax,
            x,
            y,
            SYSTEM_LABELS[system],
            size=3.25,
            color=base.INK,
            ha="left",
        )


def overlay_featured_biology(ax: plt.Axes) -> None:
    base.rounded(
        ax,
        34.0,
        66.45,
        14.5,
        2.55,
        face="#FFF2F8",
        edge="#CC78BC",
        lw=0.45,
        radius=0.5,
        alpha=0.96,
        zorder=20,
    )
    base.text(
        ax,
        41.25,
        67.73,
        "Human ovary (scONT)",
        size=3.6,
        color="#8F4B86",
        weight="bold",
    )
    base.rounded(
        ax,
        14.9,
        44.1,
        15.8,
        2.1,
        face="#FFF8E9",
        edge=base.ORANGE,
        lw=0.45,
        radius=0.5,
        alpha=0.96,
        zorder=20,
    )
    base.text(
        ax,
        22.8,
        45.15,
        "Mouse gastrula (scONT)",
        size=3.25,
        color="#9A6511",
        weight="bold",
    )


def overlay_database(ax: plt.Axes) -> None:
    base.text(ax, 78.2, 68.9, "Unified database", size=6.25, weight="bold")
    base.text(
        ax,
        78.2,
        56.9,
        "scTHREAD",
        size=10.8,
        color=base.PURPLE,
        weight="bold",
    )
    base.draw_transcript(ax, 71.7, 52.5, 13.0)
    base.text(
        ax,
        78.2,
        48.2,
        "845,781",
        size=8.8,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        78.2,
        45.35,
        "cells",
        size=5.3,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        78.2,
        41.45,
        ">200k isoforms",
        size=5.7,
        color=base.TEAL,
        weight="bold",
    )
    base.rounded(
        ax,
        66.1,
        33.0,
        24.2,
        2.9,
        face=base.WHITE,
        edge="none",
        lw=0,
        radius=0.4,
        alpha=0.94,
        zorder=20,
    )
    base.text(
        ax,
        78.2,
        34.45,
        "cell-resolved transcriptomes",
        size=4.35,
        color=base.MID,
        weight="bold",
    )


def overlay_evidence(ax: plt.Axes) -> None:
    base.rounded(
        ax,
        98.35,
        69.2,
        35.6,
        3.3,
        face=base.WHITE,
        edge="none",
        lw=0,
        radius=0.35,
        alpha=0.96,
        zorder=18,
    )
    base.text(
        ax,
        116.15,
        70.85,
        "Four linked RNA evidence layers",
        size=5.8,
        weight="bold",
    )
    rows = [
        (66.31, "Gene expression", base.PALE_TEAL),
        (55.63, "Isoform usage (DIU)", base.PALE_LAV),
        (44.23, "poly(A) site usage (APA)", "#FFF8E9"),
        (32.68, "Allelic expression (ASE)", "#FFF2EC"),
    ]
    for y, label, fill in rows:
        base.rounded(
            ax,
            117.96,
            y - 1.65,
            15.09,
            3.3,
            face=fill,
            edge="none",
            lw=0,
            radius=0.35,
            alpha=0.97,
            zorder=19,
        )
        if "Allelic" in label:
            base.text(
                ax,
                125.505,
                y + 0.65,
                "Allelic expression (ASE)",
                size=2.85,
                weight="bold",
                color=base.INK,
            )
            base.text(
                ax,
                125.505,
                y - 0.65,
                "+ splice junctions",
                size=2.95,
                weight="bold",
                color=base.INK,
            )
        else:
            sizes = {
                "Gene expression": 3.8,
                "Isoform usage (DIU)": 3.35,
                "poly(A) site usage (APA)": 2.95,
            }
            base.text(
                ax,
                125.505,
                y,
                label,
                size=sizes[label],
                weight="bold",
                color=base.INK,
            )
    base.rounded(
        ax,
        103.0,
        27.7,
        8.9,
        4.0,
        face="#FFF2EC",
        edge="none",
        lw=0,
        radius=0.35,
        alpha=0.98,
        zorder=20,
    )
    base.text(ax, 107.45, 30.55, "REF", size=3.6, color=base.PURPLE, weight="bold")
    base.text(ax, 107.45, 28.75, "ALT", size=3.6, color=base.GREEN, weight="bold")


def draw_ptprc_heatmap(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(base.COMPONENTS / "ptprc_two_isoform_switch.tsv", sep="\t")
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
    main_ax.add_patch(
        Rectangle(
            (66.6, 6.0),
            26.0,
            10.9,
            facecolor=base.WHITE,
            edgecolor="none",
            zorder=20,
        )
    )
    ax = base.add_axes_mm(fig, 70.3, 6.8, 21.0, 8.9)
    cmap = LinearSegmentedColormap.from_list(
        "ptprc-v4", ["#F7FAFC", "#BFDCEB", "#4E91BC", "#175B8D"]
    )
    ax.pcolormesh(
        np.arange(9),
        np.arange(3),
        matrix.to_numpy(),
        cmap=cmap,
        vmin=0,
        vmax=0.60,
        shading="flat",
        edgecolors=base.WHITE,
        linewidth=0.35,
    )
    ax.set_xlim(0, 8)
    ax.set_ylim(2, 0)
    ax.set_xticks(np.arange(8) + 0.5)
    ax.set_xticklabels(
        ["B", "DC", "Ery", "Mo", "NK", "Pl", "Pro", "T"], fontsize=4.0
    )
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["…67364", "…97630"], fontsize=3.6)
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
                fontsize=3.25,
                color=base.WHITE if value > 0.34 else base.INK,
                fontweight="bold" if value > 0.30 else "normal",
            )
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_inventory(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "three_axis_inventory_source.tsv", sep="\t"
    )
    order = ["DIU", "APA", "ASE"]
    source = source.set_index("axis").reindex(order).reset_index()
    main_ax.add_patch(
        Rectangle(
            (99.45, 6.0),
            34.95,
            10.9,
            facecolor=base.WHITE,
            edgecolor="none",
            zorder=20,
        )
    )
    ax = base.add_axes_mm(fig, 101.0, 7.0, 31.8, 8.6)
    colors = {"DIU": base.PURPLE, "APA": base.ORANGE, "ASE": base.TEAL}
    y = np.arange(3)
    fraction = source["sig_fraction"].to_numpy() * 100
    ax.barh(
        y,
        fraction,
        color=[colors[axis] for axis in source["axis"]],
        height=0.50,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=4.15, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 42)
    ax.set_xticks([0, 20, 40])
    ax.set_xticklabels(["0", "20", "40%"], fontsize=3.8)
    ax.tick_params(length=0, pad=0.8)
    for yi, row in source.iterrows():
        ax.text(
            max(float(row["sig_fraction"]) * 100 + 0.9, 0.9),
            yi,
            f"{int(row['n_sig_fdr05_effect_gate']):,}/{int(row['n_genes_tested']):,}",
            va="center",
            ha="left",
            fontsize=3.8,
            color=base.INK,
        )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(base.LIGHT)
    ax.spines["bottom"].set_linewidth(0.35)


def draw_umap(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "atlas_umap_stratified_sample.tsv", sep="\t"
    )
    ax = base.add_axes_mm(fig, 143.2, 21.1, 31.5, 29.8)
    for label, color in base.UMAP_COLORS.items():
        sub = source[source["broad_class"] == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap1"],
            sub["umap2"],
            s=1.75,
            c=color,
            alpha=0.72,
            linewidths=0,
            rasterized=False,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    base.text(
        main_ax,
        158.98,
        56.1,
        "Cell-type landscape",
        size=5.0,
        weight="bold",
    )
    legend_x = [145.45, 153.25, 161.15, 168.65]
    legend_y = [54.25, 52.45]
    labels = list(base.UMAP_COLORS)
    for idx, label in enumerate(labels):
        row_idx, col_idx = divmod(idx, 4)
        x = legend_x[col_idx]
        y = legend_y[row_idx]
        main_ax.scatter(
            [x - 0.52],
            [y],
            s=4.8,
            color=base.UMAP_COLORS[label],
            zorder=30,
        )
        base.text(
            main_ax,
            x,
            y,
            UMAP_LABELS[label],
            size=2.9,
            color=base.INK,
            ha="left",
        )
    base.text(
        main_ax,
        158.98,
        19.35,
        "74,906 portal cells · exploratory view",
        size=3.65,
        color=base.MID,
    )


def overlay_access_labels(ax: plt.Axes) -> None:
    modules = [
        (12.32, 18.0, "Browse"),
        (33.19, 18.0, "Search"),
        (53.01, 18.0, "Download"),
        (79.56, 26.0, "PTPRC isoform usage"),
        (116.95, 25.5, "Genome-wide events"),
    ]
    for x, width, heading in modules:
        base.rounded(
            ax,
            x - width / 2,
            17.35,
            width,
            2.9,
            face=base.WHITE,
            edge="none",
            lw=0,
            radius=0.3,
            alpha=0.94,
            zorder=20,
        )
        base.text(
            ax,
            x,
            18.8,
            heading,
            size=4.75 if x < 70 else 4.25,
            weight="bold",
        )


def overlay_query_labels(ax: plt.Axes) -> None:
    base.rounded(
        ax,
        141.58,
        69.25,
        34.8,
        3.1,
        face=base.WHITE,
        edge="none",
        lw=0,
        radius=0.35,
        alpha=0.96,
    )
    base.text(
        ax,
        158.98,
        70.8,
        "Explore genes & cell types",
        size=5.8,
        weight="bold",
    )
    base.text(
        ax,
        158.98,
        65.8,
        "PTPRC / CD45",
        size=4.8,
        color=base.INK,
        weight="bold",
    )
    base.text(
        ax,
        158.98,
        59.35,
        "Linked gene evidence",
        size=4.0,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.rounded(
        ax,
        147.03,
        9.9,
        23.9,
        3.8,
        face=base.PALE_LAV,
        edge=base.PURPLE,
        lw=0.55,
        radius=0.7,
    )
    base.text(
        ax,
        158.98,
        11.8,
        "Export · CSV · JSON · API",
        size=4.1,
        color=base.PURPLE_DARK,
        weight="bold",
    )


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v4"
    base.validate_sources()
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(
        figsize=(base.W_MM / 25.4, base.H_MM / 25.4),
        facecolor="none",
    )
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, base.W_MM)
    ax.set_ylim(0, base.H_MM)
    ax.set_aspect("equal")
    ax.axis("off")

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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, default=DEFAULT_STEM)
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()
    for output in render(args.stem.resolve(), dpi=args.dpi):
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
