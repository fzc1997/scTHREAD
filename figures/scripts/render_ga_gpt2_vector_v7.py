#!/usr/bin/env python3
"""Render the denser, overlap-free graphical abstract v7 overlay."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

import render_ga_gpt2_vector_v4 as v4
import render_ga_gpt2_vector_v5 as v5
import render_ga_gpt2_vector_v6 as v6


ROOT = Path(__file__).resolve().parents[2]
V6_DATA = ROOT / "figures/ga_gpt2_vector_v6"
V7_DATA = ROOT / "figures/ga_gpt2_vector_v7"
DEFAULT_STEM = V7_DATA / "scTHREAD_ga_gpt2_vector_v7_overlay"
base = v6.base


def format_cells(value: int) -> str:
    return f"{value / 1000:.1f}k"


def overlay_catalog(ax: plt.Axes) -> None:
    composition = pd.read_csv(
        base.COMPONENTS / "catalog_system_composition_20260726.tsv", sep="\t"
    )
    composition = composition.sort_values(
        "cells", ascending=False
    ).reset_index(drop=True)
    tissues = pd.read_csv(V6_DATA / "registry_tissue_examples.tsv", sep="\t")
    if len(tissues) != 10:
        raise AssertionError(f"Expected 10 registry examples, observed {len(tissues)}")
    tissues = tissues.sort_values(
        ["cells", "display"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)

    base.text(ax, 23.0, 70.75, "Biological coverage", size=5.8, weight="bold")
    base.text(
        ax,
        23.0,
        68.65,
        "Tissues · disease · development · two species",
        size=3.25,
        color=base.MID,
        weight="bold",
    )
    base.rounded(
        ax,
        1.55,
        24.55,
        42.2,
        28.0,
        face="#FBFBFD",
        edge="#AEB1B7",
        lw=0.52,
        radius=0.55,
    )
    base.text(
        ax,
        22.65,
        50.95,
        "469 samples · 31 studies · 845,781 cells",
        size=4.05,
        weight="bold",
        color=base.PURPLE_DARK,
    )
    base.text(
        ax,
        22.65,
        49.25,
        "Human + mouse · ONT + PacBio",
        size=3.15,
        color=base.MID,
        weight="bold",
    )

    x0, y0, total_w, bar_h = 3.05, 47.15, 39.15, 1.18
    total = float(composition["cells"].sum())
    cursor = x0
    for _, row in composition.iterrows():
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

    available = set(composition["system"].astype(str))
    legend_items = [
        system for system in base.SYSTEM_LEGEND_ORDER if system in available
    ]
    x_positions = [3.55, 13.8, 23.25, 32.55]
    for idx, system in enumerate(legend_items):
        row_idx, col_idx = divmod(idx, 4)
        x = x_positions[col_idx]
        y = [45.05, 43.1][row_idx]
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

    ax.plot(
        [3.0, 42.25],
        [41.45, 41.45],
        color="#D8DAE0",
        lw=0.35,
        zorder=22,
    )
    base.text(
        ax,
        22.65,
        40.15,
        "Representative registry datasets · cells",
        size=3.4,
        color=base.PURPLE_DARK,
        weight="bold",
    )

    y_centres = [38.15, 35.65, 33.15, 30.65, 28.15]
    columns = [(2.75, 19.15), (22.55, 19.45)]
    for idx, row in tissues.iterrows():
        col_idx, row_idx = divmod(idx, 5)
        x, width = columns[col_idx]
        y = y_centres[row_idx]
        system = str(row["system"])
        color = base.SYSTEM_COLORS[system]
        base.rounded(
            ax,
            x,
            y - 0.88,
            width,
            1.76,
            face=base.WHITE,
            edge="#D9DBE1",
            lw=0.28,
            radius=0.32,
            zorder=21,
        )
        ax.add_patch(
            Rectangle(
                (x + 0.28, y - 0.56),
                0.62,
                1.12,
                facecolor=color,
                edgecolor="none",
                zorder=23,
            )
        )
        label = str(row["display"])
        label_size = 2.85 if len(label) <= 15 else 2.6
        base.text(
            ax,
            x + 1.25,
            y,
            label,
            size=label_size,
            color=base.INK,
            weight="bold",
            ha="left",
        )
        base.text(
            ax,
            x + width - 0.55,
            y,
            format_cells(int(row["cells"])),
            size=2.7,
            color=base.MID,
            weight="bold",
            ha="right",
        )


def overlay_database(ax: plt.Axes) -> None:
    base.text(ax, 78.2, 68.9, "Unified database", size=6.25, weight="bold")

    # The incoming ribbons terminate at x=68.60 mm.  Starting the opaque
    # brand plate at x=69.0 mm keeps every ribbon outside the wordmark.
    base.rounded(
        ax,
        69.0,
        56.25,
        18.4,
        5.1,
        face=base.PURPLE_DARK,
        edge=base.PURPLE_DARK,
        lw=0,
        radius=1.35,
        zorder=27,
    )
    base.text(
        ax,
        78.2,
        58.8,
        "scTHREAD",
        size=7.8,
        color=base.WHITE,
        weight="bold",
        zorder=31,
    )
    base.draw_transcript(ax, 72.0, 53.25, 12.4)
    base.text(
        ax,
        78.2,
        48.35,
        "845,781",
        size=8.8,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        78.2,
        45.45,
        "cells",
        size=5.3,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.rounded(
        ax,
        70.0,
        39.8,
        16.4,
        4.25,
        face=base.TEAL,
        edge=base.TEAL,
        lw=0,
        radius=0.8,
        zorder=27,
    )
    base.text(
        ax,
        78.2,
        42.85,
        ">200k isoforms",
        size=3.85,
        color=base.WHITE,
        weight="bold",
        zorder=31,
    )
    base.text(
        ax,
        78.2,
        41.2,
        "cell resolved",
        size=2.8,
        color=base.WHITE,
        weight="bold",
        zorder=31,
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
        (66.31, "Gene expression", "across cell types", 3.55, 2.75),
        (
            55.63,
            "Isoform usage (DIU)",
            "cell-type transcript shifts",
            3.05,
            2.45,
        ),
        (
            44.23,
            "poly(A) site usage",
            "alternative 3′ ends",
            2.9,
            2.55,
        ),
        (
            32.68,
            "Allelic expression",
            "ASE + splice junctions",
            2.85,
            2.45,
        ),
    ]
    for y, title, subtitle, title_size, subtitle_size in rows:
        base.text(
            ax,
            125.505,
            y + 0.82,
            title,
            size=title_size,
            color=base.WHITE,
            weight="bold",
            zorder=31,
        )
        base.text(
            ax,
            125.505,
            y - 1.02,
            subtitle,
            size=subtitle_size,
            color=base.WHITE,
            weight="bold",
            zorder=31,
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


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v7"
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

    v5.overlay_header_footer(ax)
    overlay_catalog(ax)
    overlay_database(ax)
    overlay_evidence(ax)
    v6.overlay_query(ax)
    v6.overlay_bottom_cards(ax)
    v5.draw_ptprc_heatmap(fig)
    v6.draw_inventory(fig)
    v6.draw_two_isoform_umaps(fig, ax)

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
