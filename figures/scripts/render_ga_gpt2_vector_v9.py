#!/usr/bin/env python3
"""Render the scientifically corrected and geometry-polished GA v9 overlay."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle

import render_ga_gpt2_vector_v6 as v6


ROOT = Path(__file__).resolve().parents[2]
V9_DATA = ROOT / "figures/ga_gpt2_vector_v9"
DEFAULT_STEM = V9_DATA / "scTHREAD_ga_gpt2_vector_v9_overlay"
base = v6.base

SYSTEM_LABELS = {
    "Blood/immune": "Blood/immune",
    "Neural/sensory": "Neural",
    "Cancer": "Cancer/models",
    "Endocrine": "Endocrine",
    "Heart/vascular": "Heart",
    "Development/embryo": "Development",
}
SYSTEM_ORDER = [
    "Blood/immune",
    "Neural/sensory",
    "Cancer",
    "Endocrine",
    "Heart/vascular",
    "Development/embryo",
]


def validate_sources() -> None:
    required = [
        V9_DATA / "catalog_system_composition_v9.tsv",
        V9_DATA / "registry_context_examples_v9.tsv",
        V9_DATA / "three_axis_inventory_9999_v9.tsv",
        base.COMPONENTS / "ptprc_two_isoform_switch.tsv",
        v6.V6_DATA / "ptprc_two_isoform_umap_points.tsv",
        v6.V6_DATA / "ptprc_two_isoform_umap_stats.json",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing GA v9 source: " + ", ".join(missing))
    systems = pd.read_csv(required[0], sep="\t")
    if set(systems["system"]) != set(SYSTEM_ORDER):
        raise AssertionError("GA v9 biological-system classes changed")
    if int(systems["cells"].sum()) != 923_389:
        raise AssertionError("GA v9 system composition does not close")
    inventory = pd.read_csv(required[2], sep="\t").set_index("axis")
    expected = {"DIU": (2_008, 8_092), "APA": (2_558, 10_531), "ASE": (0, 6_930)}
    for axis, (passed, tested) in expected.items():
        if (
            int(inventory.loc[axis, "n_sig_fdr05_effect_gate"]) != passed
            or int(inventory.loc[axis, "n_genes_tested"]) != tested
        ):
            raise AssertionError(f"GA v9 {axis} inventory changed")


def overlay_header_footer(ax: plt.Axes) -> None:
    base.text(
        ax,
        base.W_MM / 2,
        75.20,
        "scTHREAD: a single-cell long-read transcriptome atlas",
        size=9.0,
        weight="bold",
    )
    base.text(
        ax,
        base.W_MM / 2,
        72.95,
        "https://scthread.ai4sc.ac.cn",
        size=5.25,
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
        (
            "Human & mouse tissues · cell types · isoforms · poly(A) sites · "
            "allele-aware expression"
        ),
        size=4.15,
        color=base.TEAL,
        weight="bold",
    )


def overlay_catalog(ax: plt.Axes) -> None:
    composition = pd.read_csv(
        V9_DATA / "catalog_system_composition_v9.tsv", sep="\t"
    )
    composition = composition.set_index("system").reindex(SYSTEM_ORDER).reset_index()
    contexts = pd.read_csv(
        V9_DATA / "registry_context_examples_v9.tsv", sep="\t"
    ).sort_values(
        ["cells", "display"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    if len(contexts) != 10:
        raise AssertionError(f"Expected 10 registry contexts, observed {len(contexts)}")

    base.text(ax, 23.0, 70.75, "Biological coverage", size=5.8, weight="bold")
    base.text(
        ax,
        23.0,
        68.62,
        "Tissues · cell lines · disease · development · two species",
        size=2.95,
        color=base.MID,
        weight="bold",
    )

    base.text(
        ax,
        2.7,
        67.0,
        "HUMAN TISSUES / CELLS",
        size=2.45,
        color=base.PURPLE_DARK,
        weight="bold",
        ha="left",
    )
    ax.plot([13.4, 42.2], [67.0, 67.0], color="#D8CEF0", lw=0.35, zorder=25)
    base.text(
        ax,
        2.7,
        57.65,
        "MOUSE",
        size=2.45,
        color="#A66A12",
        weight="bold",
        ha="left",
    )
    ax.plot([7.2, 13.2], [57.65, 57.65], color="#F1D6A9", lw=0.35, zorder=25)
    base.text(
        ax,
        16.0,
        57.65,
        "DEVELOPMENT",
        size=2.45,
        color="#A66A12",
        weight="bold",
        ha="left",
    )
    ax.plot([24.2, 42.2], [57.65, 57.65], color="#F1D6A9", lw=0.35, zorder=25)

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
        "453 runs · 34 datasets · 923,389 cells",
        size=3.85,
        weight="bold",
        color=base.PURPLE_DARK,
    )
    base.text(
        ax,
        22.65,
        49.25,
        "Human + mouse · ONT + PacBio",
        size=3.1,
        color=base.MID,
        weight="bold",
    )

    x0, y0, total_w, bar_h = 3.05, 47.15, 39.15, 1.18
    total = float(composition["cells"].sum())
    cursor = x0
    for row in composition.itertuples(index=False):
        width = total_w * float(row.cells) / total
        ax.add_patch(
            Rectangle(
                (cursor, y0),
                width,
                bar_h,
                facecolor=base.SYSTEM_COLORS[row.system],
                edgecolor=base.WHITE,
                linewidth=0.24,
                zorder=22,
            )
        )
        cursor += width

    x_positions = [3.55, 16.0, 29.0]
    for idx, system in enumerate(SYSTEM_ORDER):
        row_idx, col_idx = divmod(idx, 3)
        x = x_positions[col_idx]
        y = [45.05, 43.1][row_idx]
        ax.scatter(
            [x - 0.58],
            [y],
            s=5.5,
            color=base.SYSTEM_COLORS[system],
            edgecolors=base.INK,
            linewidths=0.16,
            zorder=23,
        )
        base.text(
            ax,
            x,
            y,
            SYSTEM_LABELS[system],
            size=2.95,
            color=base.INK,
            ha="left",
        )

    ax.plot([3.0, 42.25], [41.45, 41.45], color="#D8DAE0", lw=0.35, zorder=22)
    base.text(
        ax,
        22.65,
        40.15,
        "Registry-backed tissues & contexts · cells",
        size=3.3,
        color=base.PURPLE_DARK,
        weight="bold",
    )

    y_centres = [38.15, 35.65, 33.15, 30.65, 28.15]
    columns = [(2.75, 19.15), (22.55, 19.45)]
    for idx, row in contexts.iterrows():
        row_idx, col_idx = divmod(idx, 2)
        x, width = columns[col_idx]
        y = y_centres[row_idx]
        color = base.SYSTEM_COLORS[str(row["system"])]
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
        label_size = 3.0 if len(label) <= 16 else 2.7
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
            f"{int(row['cells']) / 1000:.1f}k",
            size=2.75,
            color=base.MID,
            weight="bold",
            ha="right",
        )


def overlay_database(ax: plt.Axes) -> None:
    base.text(ax, 78.2, 68.9, "Unified database", size=6.25, weight="bold")

    base.rounded(
        ax,
        70.1,
        56.45,
        16.2,
        4.8,
        face=base.PURPLE_DARK,
        edge=base.PURPLE_DARK,
        lw=0,
        radius=1.25,
        zorder=27,
    )
    base.text(
        ax,
        78.2,
        58.85,
        "scTHREAD",
        size=7.15,
        color=base.WHITE,
        weight="bold",
        zorder=31,
    )
    base.draw_transcript(ax, 72.0, 53.25, 12.4)
    base.text(
        ax,
        78.2,
        48.35,
        "923,389",
        size=8.8,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.text(
        ax,
        78.2,
        45.45,
        "called cells",
        size=4.65,
        color=base.PURPLE_DARK,
        weight="bold",
    )
    base.rounded(
        ax,
        70.1,
        39.8,
        16.2,
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
        size=3.75,
        color=base.WHITE,
        weight="bold",
        zorder=31,
    )
    base.text(
        ax,
        78.2,
        41.18,
        "cell-resolved",
        size=2.95,
        color=base.WHITE,
        weight="bold",
        zorder=31,
    )

    # Explicit provenance replaces the old unlabeled hierarchy. Draw the
    # connector first so the filled pill masks the line behind its text.
    ax.plot(
        [78.2, 78.2],
        [39.72, 20.25],
        color=base.PURPLE_DARK,
        lw=0.72,
        zorder=17,
    )
    base.rounded(
        ax,
        67.8,
        29.85,
        20.8,
        3.75,
        face=base.PALE_LAV,
        edge=base.PURPLE,
        lw=0.5,
        radius=0.75,
        zorder=25,
    )
    base.text(
        ax,
        78.2,
        31.72,
        "study → run → cell → transcript",
        size=3.2,
        color=base.PURPLE_DARK,
        weight="bold",
        zorder=31,
    )


def overlay_evidence(ax: plt.Axes) -> None:
    base.text(
        ax,
        116.15,
        71.85,
        "Transcriptome evidence",
        size=5.05,
        weight="bold",
    )
    rows = [
        (66.31, "Gene expression", "across cell types", 3.8, 3.0),
        (55.63, "Isoform usage", "cell-type shifts · DIU", 3.7, 2.85),
        (44.23, "poly(A) sites", "alternative 3′ ends", 3.7, 2.9),
        (32.68, "Allele-aware evidence", "ASE counts · junctions", 3.2, 2.8),
    ]
    for y, title, subtitle, title_size, subtitle_size in rows:
        base.text(
            ax,
            125.05,
            y + 0.82,
            title,
            size=title_size,
            color=base.WHITE,
            weight="bold",
            zorder=31,
        )
        base.text(
            ax,
            125.05,
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
        28.15,
        8.9,
        3.35,
        face="#FFF2EC",
        edge="none",
        lw=0,
        radius=0.35,
        alpha=0.98,
        zorder=20,
    )
    base.text(ax, 107.45, 30.25, "REF", size=3.45, color=base.PURPLE, weight="bold")
    base.text(ax, 107.45, 28.85, "ALT", size=3.45, color=base.GREEN, weight="bold")


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
        60.05,
        "PTPRC isoform expression",
        size=4.35,
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
        10.25,
        23.9,
        3.35,
        face=base.PALE_LAV,
        edge=base.PURPLE,
        lw=0.55,
        radius=0.7,
    )
    base.text(
        ax,
        158.98,
        11.93,
        "Export · CSV · JSON · API",
        size=4.0,
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

    base.text(ax, 15.1, 14.65, "Open data access", size=4.3, weight="bold")
    access = [
        (3.7, 7.0, "Browse", base.PURPLE),
        (11.25, 7.0, "Search", "#7456AE"),
        (18.8, 7.0, "Download", base.TEAL),
    ]
    for x, y, label, color in access:
        base.rounded(
            ax,
            x,
            y + 2.1,
            6.65,
            2.45,
            face=color,
            edge=color,
            lw=0,
            radius=0.55,
            zorder=23,
        )
        base.text(
            ax,
            x + 3.325,
            y + 3.325,
            label,
            size=2.65 if label == "Download" else 2.85,
            color=base.WHITE,
            weight="bold",
            zorder=31,
        )
    base.text(
        ax,
        15.1,
        7.8,
        "CSV · JSON · API",
        size=2.9,
        color=base.MID,
        weight="bold",
    )
    base.text(
        ax,
        59.25,
        16.2,
        "PTPRC isoform usage (% of all 23)",
        size=4.0,
        weight="bold",
    )
    base.text(
        ax,
        114.25,
        16.2,
        "Cell-type-associated genes",
        size=4.1,
        weight="bold",
    )
    base.text(
        ax,
        114.25,
        14.45,
        "FDR < 0.05 + effect-size criterion",
        size=2.85,
        color=base.MID,
        weight="bold",
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
    source = pd.read_csv(
        base.COMPONENTS / "ptprc_two_isoform_switch.tsv", sep="\t"
    )
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
    ax = base.add_axes_mm(fig, 36.0, 7.9, 50.2, 6.55)
    cmap = LinearSegmentedColormap.from_list(
        "ptprc-v9", ["#F7FAFC", "#BFDCEB", "#4E91BC", "#175B8D"]
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
        ["B", "DC", "Ery", "Mo", "NK", "Pl", "Pro", "T"], fontsize=3.8
    )
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["ENST…67364", "ENST…97630"], fontsize=3.2)
    ax.tick_params(length=0, pad=0.55)
    for i in range(2):
        for j in range(8):
            value = matrix.iat[i, j]
            ax.text(
                j + 0.5,
                i + 0.5,
                f"{value * 100:.0f}",
                ha="center",
                va="center",
                fontsize=3.2,
                color=base.WHITE if value > 0.34 else base.INK,
                fontweight="bold" if value > 0.30 else "normal",
            )
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_inventory(fig: plt.Figure, main_ax: plt.Axes) -> None:
    source = pd.read_csv(
        V9_DATA / "three_axis_inventory_9999_v9.tsv", sep="\t"
    )
    source = source.set_index("axis").reindex(["DIU", "APA", "ASE"]).reset_index()
    ax = base.add_axes_mm(fig, 94.25, 8.25, 40.75, 5.0)
    colors = {"DIU": base.PURPLE, "APA": base.ORANGE, "ASE": base.CORAL}
    y = np.arange(3)
    fraction = source["sig_fraction"].to_numpy() * 100
    ax.barh(
        y,
        fraction,
        color=[colors[axis] for axis in source["axis"]],
        height=0.42,
    )
    # A hollow marker makes the zero ASE row visually explicit without
    # suggesting that the database contains no ASE evidence.
    ax.scatter(
        [0.35],
        [2],
        s=7.5,
        facecolors=base.WHITE,
        edgecolors=base.CORAL,
        linewidths=0.8,
        zorder=4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(["DIU", "APA", "ASE"], fontsize=3.7, fontweight="bold")
    ax.set_xlim(0, 45)
    ax.set_ylim(2.45, -0.50)
    ax.set_xticks([])
    ax.tick_params(length=0, pad=0.6)
    for yi, row in source.iterrows():
        ax.text(
            max(float(row["sig_fraction"]) * 100 + 0.75, 1.25),
            yi,
            (
                f"{int(row['n_sig_fdr05_effect_gate']):,} / "
                f"{int(row['n_genes_tested']):,} genes"
            ),
            va="center",
            ha="left",
            fontsize=3.3,
            color=base.INK,
        )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ase = source.set_index("axis").loc["ASE"]
    base.text(
        main_ax,
        94.65,
        7.05,
        (
            f"ASE: {int(ase['n_raw_p_lt_0.05']):,} nominal P < 0.05; "
            f"0 after FDR + effect filter (min q = {float(ase['min_q']):.3f})"
        ),
        size=2.35,
        color=base.MID,
        ha="left",
    )


def render(stem: Path, dpi: int = 450) -> list[Path]:
    base.setup_style()
    base.mpl.rcParams["svg.hashsalt"] = "scTHREAD-ga-v9"
    validate_sources()
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
    overlay_database(ax)
    overlay_evidence(ax)
    overlay_query(ax)
    overlay_bottom_cards(ax)
    draw_ptprc_heatmap(fig)
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
