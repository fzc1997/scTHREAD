#!/usr/bin/env python3
"""Render the evidence-first, reduced-icon graphical abstract v5 overlay."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle

import render_ga_gpt2_vector_v4 as v4


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STEM = (
    ROOT / "figures/ga_gpt2_vector_v5/scTHREAD_ga_gpt2_vector_v5_overlay"
)
base = v4.base


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


def overlay_catalog(ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "catalog_system_composition_20260726.tsv", sep="\t"
    )
    source = source.sort_values("cells", ascending=False).reset_index(drop=True)

    base.text(ax, 23.0, 70.75, "Biological coverage", size=5.8, weight="bold")
    base.rounded(
        ax,
        1.55,
        28.3,
        42.2,
        14.2,
        face=base.WHITE,
        edge="#AEB1B7",
        lw=0.5,
        radius=0.45,
    )
    base.text(
        ax,
        22.65,
        40.85,
        "469 samples · 31 studies · 845,781 cells",
        size=4.05,
        weight="bold",
        color=base.PURPLE_DARK,
    )
    base.text(
        ax,
        22.65,
        39.22,
        "Human + mouse · ONT + PacBio",
        size=3.15,
        color=base.MID,
        weight="bold",
    )

    x0, y0, total_w, bar_h = 3.05, 36.95, 39.15, 1.15
    total = float(source["cells"].sum())
    cursor = x0
    for _, row in source.iterrows():
        system = str(row["system"])
        width = total_w * float(row["cells"]) / total
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
        y = [34.42, 31.78][row_idx]
        ax.scatter(
            [x - 0.58],
            [y],
            s=5.4,
            color=base.SYSTEM_COLORS[system],
            zorder=23,
        )
        base.text(
            ax,
            x,
            y,
            v4.SYSTEM_LABELS[system],
            size=3.15,
            color=base.INK,
            ha="left",
        )


def overlay_featured_biology(ax: plt.Axes) -> None:
    base.rounded(
        ax,
        34.1,
        56.05,
        14.5,
        2.45,
        face="#FFF2F8",
        edge="#CC78BC",
        lw=0.42,
        radius=0.48,
        alpha=0.97,
        zorder=20,
    )
    base.text(
        ax,
        41.35,
        57.28,
        "Human ovary (scONT)",
        size=3.35,
        color="#8F4B86",
        weight="bold",
    )
    base.rounded(
        ax,
        20.7,
        47.15,
        16.4,
        2.25,
        face="#FFF8E9",
        edge=base.ORANGE,
        lw=0.42,
        radius=0.48,
        alpha=0.97,
        zorder=20,
    )
    base.text(
        ax,
        28.9,
        48.28,
        "Mouse gastrula (scONT)",
        size=3.2,
        color="#9A6511",
        weight="bold",
    )


def overlay_query(ax: plt.Axes) -> None:
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
        59.25,
        "PTPRC-linked transcriptome evidence",
        size=3.55,
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


def overlay_bottom_cards(ax: plt.Axes) -> None:
    cards = [
        (2.6, 6.1, 25.0, 11.5, base.PALE_LAV, "#B39DDA"),
        (30.0, 6.1, 58.5, 11.5, base.PALE_LAV, "#C0AEDF"),
        (91.0, 6.1, 46.5, 11.5, "#FFF8ED", "#EFB868"),
    ]
    for x, y, width, height, face, edge in cards:
        base.rounded(
            ax,
            x,
            y,
            width,
            height,
            face=face,
            edge=edge,
            lw=0.5,
            radius=0.7,
            alpha=0.98,
            zorder=18,
        )

    base.text(ax, 15.1, 14.55, "Open data access", size=4.35, weight="bold")
    base.text(
        ax,
        15.1,
        11.45,
        "Browse · Search · Download",
        size=3.35,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        15.1,
        8.75,
        "CSV · JSON · API",
        size=3.15,
        color=base.MID,
        weight="bold",
    )
    base.text(
        ax,
        59.25,
        16.25,
        "PTPRC / CD45 isoform usage",
        size=4.25,
        weight="bold",
    )
    base.text(
        ax,
        114.25,
        16.25,
        "Genome-wide RNA events",
        size=4.25,
        weight="bold",
    )

    ax.plot(
        [78.2, 78.2],
        [31.1, 20.25],
        color=base.PURPLE_DARK,
        lw=0.75,
        zorder=17,
    )
    ax.plot(
        [15.1, 114.25],
        [20.25, 20.25],
        color=base.PURPLE_DARK,
        lw=0.75,
        zorder=17,
    )
    for x in (15.1, 59.25, 114.25):
        ax.add_patch(
            FancyArrowPatch(
                (x, 20.25),
                (x, 17.82),
                arrowstyle="-|>",
                mutation_scale=6.5,
                linewidth=0.75,
                color=base.PURPLE_DARK,
                zorder=17,
            )
        )


def draw_ptprc_heatmap(fig: plt.Figure) -> None:
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
    ax = base.add_axes_mm(fig, 36.0, 7.0, 50.2, 7.75)
    cmap = LinearSegmentedColormap.from_list(
        "ptprc-v5", ["#F7FAFC", "#BFDCEB", "#4E91BC", "#175B8D"]
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
        ["B", "DC", "Ery", "Mo", "NK", "Pl", "Pro", "T"], fontsize=3.9
    )
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["ENST…67364", "ENST…97630"], fontsize=3.25)
    ax.tick_params(length=0, pad=0.65)
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


def draw_inventory(fig: plt.Figure) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "three_axis_inventory_source.tsv", sep="\t"
    )
    source = source.set_index("axis").reindex(["DIU", "APA", "ASE"]).reset_index()
    ax = base.add_axes_mm(fig, 94.4, 7.0, 40.6, 7.8)
    colors = {"DIU": base.PURPLE, "APA": base.ORANGE, "ASE": base.TEAL}
    y = np.arange(3)
    fraction = source["sig_fraction"].to_numpy() * 100
    ax.barh(
        y,
        fraction,
        color=[colors[axis] for axis in source["axis"]],
        height=0.48,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(["DIU", "APA", "ASE"], fontsize=3.9, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 45)
    ax.set_xticks([0, 20, 40])
    ax.set_xticklabels(["0", "20", "40%"], fontsize=3.65)
    ax.tick_params(length=0, pad=0.7)
    for yi, row in source.iterrows():
        ax.text(
            max(float(row["sig_fraction"]) * 100 + 0.9, 0.9),
            yi,
            f"{int(row['n_sig_fdr05_effect_gate']):,}/{int(row['n_genes_tested']):,}",
            va="center",
            ha="left",
            fontsize=3.45,
            color=base.INK,
        )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(base.LIGHT)
    ax.spines["bottom"].set_linewidth(0.35)


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v5"
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
    overlay_catalog(ax)
    overlay_featured_biology(ax)
    v4.overlay_database(ax)
    v4.overlay_evidence(ax)
    overlay_query(ax)
    overlay_bottom_cards(ax)
    draw_ptprc_heatmap(fig)
    draw_inventory(fig)
    v4.draw_umap(fig, ax)

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
