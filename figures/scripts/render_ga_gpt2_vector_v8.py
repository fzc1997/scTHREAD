#!/usr/bin/env python3
"""Render the reading-order and annotation-polished graphical abstract v8."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

import render_ga_gpt2_vector_v4 as v4
import render_ga_gpt2_vector_v5 as v5
import render_ga_gpt2_vector_v6 as v6
import render_ga_gpt2_vector_v7 as v7


ROOT = Path(__file__).resolve().parents[2]
V6_DATA = ROOT / "figures/ga_gpt2_vector_v6"
V8_DATA = ROOT / "figures/ga_gpt2_vector_v8"
DEFAULT_STEM = V8_DATA / "scTHREAD_ga_gpt2_vector_v8_overlay"
base = v7.base


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

    # Lightweight rules organize the pictograms by species without adding
    # cohort-specific callouts or competing with the biological icons.
    base.text(
        ax,
        2.7,
        67.0,
        "HUMAN TISSUES",
        size=2.55,
        color=base.PURPLE_DARK,
        weight="bold",
        ha="left",
    )
    ax.plot(
        [10.25, 42.2],
        [67.0, 67.0],
        color="#D8CEF0",
        lw=0.35,
        zorder=25,
    )
    base.text(
        ax,
        2.7,
        57.65,
        "MOUSE / DEVELOPMENT",
        size=2.55,
        color="#A66A12",
        weight="bold",
        ha="left",
    )
    ax.plot(
        [14.15, 42.2],
        [57.65, 57.65],
        color="#F1D6A9",
        lw=0.35,
        zorder=25,
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
        row_idx, col_idx = divmod(idx, 2)
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
            v7.format_cells(int(row["cells"])),
            size=2.7,
            color=base.MID,
            weight="bold",
            ha="right",
        )


def draw_inventory(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "three_axis_inventory_source.tsv", sep="\t"
    )
    source = source.set_index("axis").reindex(["DIU", "APA", "ASE"]).reset_index()
    ax = base.add_axes_mm(fig, 94.25, 8.15, 40.75, 5.35)
    colors = {"DIU": base.PURPLE, "APA": base.ORANGE, "ASE": base.TEAL}
    y = np.arange(3)
    fraction = source["sig_fraction"].to_numpy() * 100
    ax.barh(
        y,
        fraction,
        color=[colors[axis] for axis in source["axis"]],
        height=0.42,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(["DIU", "APA", "ASE"], fontsize=3.7, fontweight="bold")
    ax.set_xlim(0, 45)
    ax.set_ylim(2.45, -0.50)
    ax.set_xticks([])
    ax.tick_params(length=0, pad=0.6)
    for yi, row in source.iterrows():
        ax.text(
            max(float(row["sig_fraction"]) * 100 + 0.75, 0.75),
            yi,
            (
                f"{int(row['n_sig_fdr05_effect_gate']):,}/"
                f"{int(row['n_genes_tested']):,} passed"
            ),
            va="center",
            ha="left",
            fontsize=3.3,
            color=base.INK,
        )
    for spine in ax.spines.values():
        spine.set_visible(False)
    base.text(
        main_ax,
        94.65,
        7.05,
        "ASE: 673 raw P<0.05 · min q=0.639",
        size=2.65,
        color=base.MID,
        ha="left",
    )


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v8"
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
    v7.overlay_database(ax)
    v7.overlay_evidence(ax)
    v6.overlay_query(ax)
    v6.overlay_bottom_cards(ax)
    v5.draw_ptprc_heatmap(fig)
    draw_inventory(fig, ax)
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
