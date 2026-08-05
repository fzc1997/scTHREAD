#!/usr/bin/env python3
"""Render the evidence-led graphical abstract v6 overlay from frozen data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle

import render_ga_gpt2_vector_v4 as v4
import render_ga_gpt2_vector_v5 as v5


ROOT = Path(__file__).resolve().parents[2]
V6_DATA = ROOT / "figures/ga_gpt2_vector_v6"
DEFAULT_STEM = V6_DATA / "scTHREAD_ga_gpt2_vector_v6_overlay"
base = v5.base


def overlay_catalog(ax: plt.Axes) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "catalog_system_composition_20260726.tsv", sep="\t"
    )
    source = source.sort_values("cells", ascending=False).reset_index(drop=True)
    tissues = pd.read_csv(V6_DATA / "registry_tissue_examples.tsv", sep="\t")
    if len(tissues) != 10:
        raise AssertionError(f"Expected 10 registry examples, observed {len(tissues)}")

    base.text(ax, 23.0, 70.75, "Biological coverage", size=5.8, weight="bold")
    base.rounded(
        ax,
        1.55,
        24.8,
        42.2,
        18.2,
        face=base.WHITE,
        edge="#AEB1B7",
        lw=0.5,
        radius=0.45,
    )
    base.text(
        ax,
        22.65,
        41.45,
        "469 samples · 31 studies · 845,781 cells",
        size=4.05,
        weight="bold",
        color=base.PURPLE_DARK,
    )
    base.text(
        ax,
        22.65,
        39.88,
        "Human + mouse · ONT + PacBio",
        size=3.15,
        color=base.MID,
        weight="bold",
    )

    x0, y0, total_w, bar_h = 3.05, 37.8, 39.15, 1.12
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
        y = [35.55, 33.55][row_idx]
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
            size=3.2,
            color=base.INK,
            ha="left",
        )

    labels = tissues["display"].astype(str).tolist()
    base.text(
        ax,
        22.65,
        30.95,
        "Representative tissues & contexts",
        size=3.55,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        22.65,
        28.95,
        " · ".join(labels[:5]),
        size=3.15,
        color=base.INK,
    )
    base.text(
        ax,
        22.65,
        27.15,
        " · ".join(labels[5:]),
        size=3.15,
        color=base.INK,
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
        59.2,
        "PTPRC isoform expression",
        size=4.7,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        158.98,
        57.2,
        "same 71,913-cell embedding",
        size=3.25,
        color=base.MID,
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
        "Genes with cell-type effects",
        size=4.15,
        weight="bold",
    )
    base.text(
        ax,
        114.25,
        14.45,
        "FDR < 0.05 + prespecified effect gate",
        size=2.85,
        color=base.MID,
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


def expression_cmap(name: str, colors: list[str]) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        name,
        [
            (0.00, "#D9DDE3"),
            (0.18, colors[0]),
            (0.52, colors[1]),
            (1.00, colors[2]),
        ],
    )


def draw_color_scale(
    fig: plt.Figure,
    main_ax: plt.Axes,
    x: float,
    width: float,
    cap: int,
    cmap: LinearSegmentedColormap,
) -> None:
    scale_ax = base.add_axes_mm(fig, x, 27.9, width, 0.72)
    edges = np.linspace(0, cap, 33)
    values = ((edges[:-1] + edges[1:]) / 2)[None, :]
    scale_ax.pcolormesh(
        edges,
        np.array([0.0, 1.0]),
        values,
        cmap=cmap,
        vmin=0,
        vmax=cap,
        shading="flat",
        linewidth=0,
        rasterized=False,
    )
    scale_ax.set_xlim(0, cap)
    scale_ax.set_ylim(0, 1)
    scale_ax.axis("off")
    base.text(
        main_ax,
        x + width / 2,
        26.7,
        f"0–{cap} molecules/cell",
        size=2.85,
        color=base.MID,
        weight="bold",
    )


def draw_two_isoform_umaps(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(
        V6_DATA / "ptprc_two_isoform_umap_points.tsv", sep="\t"
    )
    stats = json.loads(
        (V6_DATA / "ptprc_two_isoform_umap_stats.json").read_text()
    )
    if len(source) != stats["retained_cells"] or len(source) != 11_472:
        raise AssertionError("Frozen PTPRC UMAP sample count changed")

    xmin, xmax = source["umap1"].min(), source["umap1"].max()
    ymin, ymax = source["umap2"].min(), source["umap2"].max()
    xpad = (xmax - xmin) * 0.025
    ypad = (ymax - ymin) * 0.025
    panels = [
        (
            142.65,
            "ENST…67364",
            "expr_67364",
            int(stats["ENST00000367364"]["display_cap"]),
            expression_cmap(
                "ptprc-67364-v6",
                ["#C9BCE6", "#7F66BA", "#3E2B75"],
            ),
        ),
        (
            160.0,
            "ENST…97630",
            "expr_97630",
            int(stats["ENST00000697630"]["display_cap"]),
            expression_cmap(
                "ptprc-97630-v6",
                ["#A7DAD6", "#39A39F", "#07565E"],
            ),
        ),
    ]
    panel_width = 15.25
    for x, title, column, cap, cmap in panels:
        base.text(
            main_ax,
            x + panel_width / 2,
            55.1,
            title,
            size=3.45,
            weight="bold",
            color=base.INK,
        )
        map_ax = base.add_axes_mm(fig, x, 30.0, panel_width, 23.8)
        zero = source[column].eq(0)
        map_ax.scatter(
            source.loc[zero, "umap1"],
            source.loc[zero, "umap2"],
            s=0.52,
            c="#D9DDE3",
            alpha=0.42,
            linewidths=0,
            rasterized=False,
        )
        positive = source.loc[~zero].sort_values(column, kind="mergesort")
        map_ax.scatter(
            positive["umap1"],
            positive["umap2"],
            s=0.82,
            c=np.minimum(positive[column].to_numpy(), cap),
            cmap=cmap,
            vmin=0,
            vmax=cap,
            alpha=0.84,
            linewidths=0,
            rasterized=False,
        )
        map_ax.set_xlim(xmin - xpad, xmax + xpad)
        map_ax.set_ylim(ymin - ypad, ymax + ypad)
        map_ax.set_aspect("equal", adjustable="box")
        map_ax.set_xticks([])
        map_ax.set_yticks([])
        for spine in map_ax.spines.values():
            spine.set_visible(False)
        draw_color_scale(fig, main_ax, x, panel_width, cap, cmap)

    base.text(
        main_ax,
        158.98,
        24.55,
        "all positive cells shown",
        size=2.95,
        color=base.INK,
        weight="bold",
    )
    base.text(
        main_ax,
        158.98,
        22.85,
        "+ stratified double-zero background",
        size=2.8,
        color=base.MID,
    )


def draw_inventory(fig: plt.Figure) -> None:
    source = pd.read_csv(
        base.COMPONENTS / "three_axis_inventory_source.tsv", sep="\t"
    )
    source = source.set_index("axis").reindex(["DIU", "APA", "ASE"]).reset_index()
    ax = base.add_axes_mm(fig, 94.25, 6.75, 40.75, 6.65)
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
    ax.set_ylim(2.65, -0.50)
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
    ax.text(
        0.75,
        2.38,
        "673 raw P<0.05 · min q=0.639",
        va="center",
        ha="left",
        fontsize=2.75,
        color=base.MID,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v6"
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
    v4.overlay_database(ax)
    v4.overlay_evidence(ax)
    overlay_query(ax)
    overlay_bottom_cards(ax)
    v5.draw_ptprc_heatmap(fig)
    draw_inventory(fig)
    draw_two_isoform_umaps(fig, ax)

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
