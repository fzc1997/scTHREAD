#!/usr/bin/env python3
"""Build a dense nine-panel scTHREAD biological-application main figure."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))

from make_biology_main_figure import (
    CELL_LABELS,
    CELL_ORDER,
    DATASETS,
    GENES,
    draw_junction_icon,
    load_cell_profiles,
    load_effects,
    panel_label,
    setup_style,
)


PANEL_IDS = tuple("abcdefghi")
FIGURE_SIZE_MM = (183, 185)
ROBUSTNESS_COLUMNS = (
    "all",
    "low_cycle",
    "low_common_activation",
    "low_both",
    "classical_like",
    "nonclassical_like",
)
GENE_COLORS = {"CCND3": "#0173B2", "CUTA": "#029E73", "PSMB8": "#9C4E97"}


def load_robustness(path: Path) -> pd.DataFrame:
    raw = json.loads(path.read_text(encoding="utf-8"))["genes"]
    rows = []
    for gene in GENES:
        baseline = raw[gene]["all"]["median_delta"]
        values = {
            "all": 1.0,
            "low_cycle": raw[gene]["low_cycle"]["effect_retention_vs_all"],
            "low_common_activation": raw[gene]["low_common_activation"][
                "effect_retention_vs_all"
            ],
            "low_both": raw[gene]["low_both"]["effect_retention_vs_all"],
            "classical_like": (
                raw[gene]["monocyte_subtype_sensitivity"]["classical_like"][
                    "median_delta"
                ]
                / baseline
            ),
            "nonclassical_like": (
                raw[gene]["monocyte_subtype_sensitivity"]["nonclassical_like"][
                    "median_delta"
                ]
                / baseline
            ),
        }
        for condition, relative_effect in values.items():
            rows.append(
                {
                    "gene": gene,
                    "condition": condition,
                    "relative_effect": relative_effect,
                }
            )
    frame = pd.DataFrame(rows)
    if not frame["relative_effect"].between(0.5, 1.6).all():
        raise ValueError("Unexpected robustness effect outside plotting contract")
    return frame


def load_breadth(path: Path, effects: pd.DataFrame | None = None) -> pd.DataFrame:
    """Panel i: the monocyte-minus-T contrast against the pre-declared extension.

    The extension summary has no monocyte-versus-T key; its `myeloid_vs_lymphoid`
    entry pools neutrophils into the myeloid side and B, NK, CD4 and CD8 into the
    lymphoid side. Labelling that "Monocyte - T" made panel i disagree with
    panels e and f by 0.28 at CCND3 in the same cohort, so the monocyte-minus-T
    series is now taken from the short-read per-donor effects that panel e plots.
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    if effects is not None:
        sr = effects[effects["dataset"] == "Independent SR"]
        if sr.empty:
            raise ValueError("Short-read effects missing; cannot build panel i")
        for gene in GENES:
            g = sr[sr["gene"] == gene]["effect"]
            rows.append({
                "comparison": "Monocyte − T",
                "gene": gene,
                "median_effect": float(g.median()),
                "estimable_donors": int(g.notna().sum()),
                "positive_donors": int((g > 0).sum()),
            })
    for label, key in (
        (None, None),
        ("Neutrophil − B/NK", "primary_neutrophil_vs_bnk"),
    ):
        if label is None:
            continue
        for gene in GENES:
            value = raw[key][gene]
            rows.append(
                {
                    "comparison": label,
                    "gene": gene,
                    "median_effect": value["median_effect"],
                    "estimable_donors": value["estimable_donors"],
                    "positive_donors": value["positive_donors"],
                }
            )
    frame = pd.DataFrame(rows)
    cuta = frame.loc[
        (frame["comparison"] == "Neutrophil − B/NK") & (frame["gene"] == "CUTA")
    ].iloc[0]
    if int(cuta["estimable_donors"]) != 1 or float(cuta["median_effect"]) >= 0:
        raise ValueError("Frozen CUTA scope boundary changed")
    return frame


def draw_panel_a(ax) -> None:
    panel_label(ax, "a", x=-0.02, y=1.03)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    draw_junction_icon(
        ax,
        0.01,
        "CCND3",
        "ENST00000372991",
        "ENST00000372988 / ENST00000415497",
        "two T-associated junctions",
    )
    draw_junction_icon(ax, 0.34, "CUTA", "ENST00000482684", "ENST00000488034")
    draw_junction_icon(ax, 0.67, "PSMB8", "ENST00000374882", "ENST00000374881")
    ax.text(
        0.01,
        0.00,
        "Exact-junction carrier classes; symbolic exon diagrams, not to scale",
        fontsize=4.8,
        color="#666666",
        va="bottom",
    )


def draw_panel_b(ax) -> None:
    panel_label(ax, "b", x=-0.08, y=1.03)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    stages = (
        (0.16, "Discovery", "long read", "13 donors", "#E7F2F8"),
        (0.50, "Held back", "long read", "4 donors", "#E7F2F8"),
        (0.84, "Independent", "short read", "20 donors", "#FFF0D9"),
    )
    for index, (x, title, platform, ntext, color) in enumerate(stages):
        ax.add_patch(
            Rectangle(
                (x - 0.13, 0.38),
                0.26,
                0.42,
                facecolor=color,
                edgecolor="#555555",
                lw=0.55,
            )
        )
        ax.text(x, 0.68, title, fontsize=6.2, fontweight="bold", ha="center")
        ax.text(x, 0.54, platform, fontsize=5.5, ha="center")
        ax.text(x, 0.43, ntext, fontsize=5.1, color="#555555", ha="center")
        if index < 2:
            ax.annotate(
                "",
                xy=(stages[index + 1][0] - 0.145, 0.59),
                xytext=(x + 0.145, 0.59),
                arrowprops={"arrowstyle": "->", "lw": 0.6, "color": "#555555"},
            )
    ax.text(
        0.5,
        0.17,
        r"endpoint: $\Delta$ within-gene exact-junction usage",
        fontsize=5.6,
        ha="center",
    )
    ax.text(
        0.5,
        0.04,
        "short reads do not establish full-length linkage",
        fontsize=4.8,
        color="#666666",
        ha="center",
    )


def draw_donor_panel(
    ax: mpl.axes.Axes,
    effects: pd.DataFrame,
    dataset: str,
    color: str,
    label: str,
    title: str,
    statistic: str,
    show_ylabel: bool,
) -> None:
    panel_label(ax, label, x=-0.15, y=1.04)
    ax.axhline(0, color="#777777", lw=0.55, ls=(0, (2, 2)), zorder=0)
    # Use rank-ordered deterministic offsets rather than random jitter. This
    # keeps the donor distribution legible in the narrow panel and makes
    # repeated renders pixel-stable.
    donor_count = None
    for index, gene in enumerate(GENES):
        values = effects.loc[
            (effects["dataset"] == dataset) & (effects["gene"] == gene), "effect"
        ].to_numpy()
        donor_count = len(values)
        jitter = np.linspace(-0.085, 0.085, len(values)) if len(values) > 1 else np.array([0.0])
        order = np.argsort(values)
        x_positions = np.empty(len(values))
        x_positions[order] = index + jitter
        marker = "s" if dataset == "Independent SR" else "o"
        ax.scatter(
            x_positions,
            values,
            s=11,
            marker=marker,
            facecolor=color,
            edgecolor="white",
            linewidth=0.35,
            alpha=0.92,
            zorder=3,
        )
        median = float(np.median(values))
        ax.plot([index - 0.14, index + 0.14], [median, median], color="black", lw=1.15, solid_capstyle="round", zorder=4)
    ax.set_title(title, loc="left", fontsize=6.25, fontweight="bold", pad=7)
    ax.text(
        0.98,
        1.075,
        statistic,
        transform=ax.transAxes,
        fontsize=4.7,
        color="#555555",
        ha="right",
        va="bottom",
    )
    ax.set_xticks(range(3), GENES, fontweight="bold")
    ax.set_xlim(-0.45, 2.45)
    ax.set_ylim(-0.08, 1.08)
    ax.set_yticks([0, 0.5, 1.0])
    if show_ylabel:
        ax.set_ylabel("Monocyte − T\njunction usage")
    else:
        ax.set_ylabel("")


def draw_panel_f(ax, effects: pd.DataFrame) -> pd.DataFrame:
    panel_label(ax, "f", x=-0.14, y=1.04)
    medians = (
        effects.groupby(["dataset", "gene"], sort=False)["effect"]
        .median()
        .rename("median_effect")
        .reset_index()
    )
    xlabels = [item[0] for item in DATASETS]
    for gene in GENES:
        values = [
            float(
                medians.loc[
                    (medians["dataset"] == dataset) & (medians["gene"] == gene),
                    "median_effect",
                ].iloc[0]
            )
            for dataset in xlabels
        ]
        ax.plot(
            range(3),
            values,
            marker="o",
            ms=3.3,
            lw=1.1,
            color=GENE_COLORS[gene],
            label=gene,
        )
    ax.set_xticks(range(3), ["Discovery\nLR", "Held-back\nLR", "Independent\nSR"])
    ax.set_ylabel("Median effect")
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_title("Median effect by cohort", loc="left", fontsize=6.4, fontweight="bold")
    ax.legend(frameon=False, fontsize=4.8, loc="upper left", ncol=1)
    return medians


def draw_panel_g(ax, robustness: pd.DataFrame) -> None:
    panel_label(ax, "g", x=-0.09, y=1.04)
    matrix = np.empty((3, len(ROBUSTNESS_COLUMNS)))
    for i, gene in enumerate(GENES):
        for j, condition in enumerate(ROBUSTNESS_COLUMNS):
            matrix[i, j] = robustness.loc[
                (robustness["gene"] == gene)
                & (robustness["condition"] == condition),
                "relative_effect",
            ].iloc[0]
    cmap = mpl.colormaps["viridis_r"].copy()
    mesh = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1) - 0.5,
        np.arange(matrix.shape[0] + 1) - 0.5,
        matrix,
        cmap=cmap,
        norm=Normalize(0.6, 1.5),
        shading="flat",
        rasterized=False,
    )
    for i in range(3):
        for j in range(len(ROBUSTNESS_COLUMNS)):
            color = "white" if matrix[i, j] > 1.2 else "black"
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=4.8, color=color)
    ax.set_xticks(
        range(len(ROBUSTNESS_COLUMNS)),
        ["All", "Low\ncycle", "Low\nactivation", "Low\nboth", "Classical\nmono.", "Nonclassical\nmono."],
    )
    ax.set_yticks(range(3), GENES, fontweight="bold")
    ax.set_xlim(-0.5, len(ROBUSTNESS_COLUMNS) - 0.5)
    ax.set_ylim(2.5, -0.5)
    ax.tick_params(length=0)
    ax.set_title("Pre-declared state and subtype sensitivities", loc="left", fontsize=6.4, fontweight="bold")
    colorbar = ax.figure.colorbar(
        mesh,
        ax=ax,
        orientation="horizontal",
        fraction=0.08,
        pad=0.23,
        boundaries=np.linspace(0.6, 1.5, 37),
    )
    colorbar.set_ticks([0.6, 1.0, 1.5])
    colorbar.set_label("Effect relative to all cells", fontsize=5.2)
    colorbar.ax.tick_params(labelsize=4.6, width=0.4, length=2)


def draw_panel_h(ax, profiles: pd.DataFrame) -> None:
    panel_label(ax, "h", x=-0.045, y=1.04)
    matrix = np.empty((3, len(CELL_ORDER)))
    counts = np.empty_like(matrix)
    for i, gene in enumerate(GENES):
        for j, cell_type in enumerate(CELL_ORDER):
            row = profiles.loc[
                (profiles["gene"] == gene) & (profiles["cell_type"] == cell_type)
            ].iloc[0]
            matrix[i, j] = row["median_usage"]
            counts[i, j] = row["estimable_donors"]
    mesh = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1) - 0.5,
        np.arange(matrix.shape[0] + 1) - 0.5,
        matrix,
        cmap=mpl.colormaps["viridis_r"],
        norm=Normalize(0, 1),
        shading="flat",
        rasterized=False,
    )
    for i in range(3):
        for j in range(len(CELL_ORDER)):
            value = matrix[i, j]
            n_value = int(counts[i, j])
            color = "white" if value >= 0.55 else "black"
            ax.text(j, i - 0.10, f"{value:.2f}", ha="center", fontsize=4.7, color=color)
            ax.text(j, i + 0.23, f"n={n_value}", ha="center", fontsize=4.0, color=color)
            if n_value < 5:
                ax.add_patch(
                    Rectangle(
                        (j - 0.49, i - 0.49),
                        0.98,
                        0.98,
                        fill=False,
                        hatch="////",
                        edgecolor="#CC3311",
                        lw=0.65,
                    )
                )
    ax.set_xticks(range(len(CELL_ORDER)), CELL_LABELS, rotation=42, ha="right")
    ax.set_yticks(range(3), GENES, fontweight="bold")
    ax.set_xlim(-0.5, len(CELL_ORDER) - 0.5)
    ax.set_ylim(2.5, -0.5)
    ax.tick_params(length=0)
    ax.set_title("Cell-class topology", loc="left", fontsize=6.4, fontweight="bold")
    ax.text(
        1,
        1.035,
        "hatched: n<5",
        transform=ax.transAxes,
        fontsize=4.6,
        color="#666666",
        ha="right",
    )
    colorbar = ax.figure.colorbar(
        mesh,
        ax=ax,
        orientation="horizontal",
        fraction=0.07,
        pad=0.28,
        boundaries=np.linspace(0, 1, 41),
    )
    colorbar.set_ticks([0, 0.5, 1])
    colorbar.set_label("Median monocyte-associated carrier usage", fontsize=5.2)
    colorbar.ax.tick_params(labelsize=4.6, width=0.4, length=2)


def draw_panel_i(ax, breadth: pd.DataFrame) -> None:
    panel_label(ax, "i", x=-0.08, y=1.04)
    ybase = np.arange(3)
    styles = (
        ("Monocyte − T", 0.12, "#0173B2", "o"),
        ("Neutrophil − B/NK", -0.12, "#666666", "s"),
    )
    ax.axvline(0, color="#777777", lw=0.55, ls=(0, (2, 2)))
    for comparison, offset, color, marker in styles:
        subset = breadth.loc[breadth["comparison"] == comparison].set_index("gene")
        for i, gene in enumerate(GENES):
            row = subset.loc[gene]
            is_boundary = comparison.startswith("Neutrophil") and gene == "CUTA"
            ax.scatter(
                row["median_effect"],
                i + offset,
                s=18,
                marker=marker,
                facecolor="none" if is_boundary else color,
                edgecolor="#CC3311" if is_boundary else color,
                linewidth=0.9 if is_boundary else 0.5,
                label=comparison if i == 0 else None,
                zorder=3,
            )
            ax.text(
                row["median_effect"] + 0.035,
                i + offset,
                f'{int(row["positive_donors"])}/{int(row["estimable_donors"])}',
                fontsize=4.3,
                va="center",
                color="#CC3311" if is_boundary else color,
            )
    ax.set_yticks(ybase, GENES, fontweight="bold")
    ax.set_ylim(2.55, -0.55)
    ax.set_xlim(-0.36, 1.02)
    ax.set_xticks([-0.25, 0, 0.5, 1.0])
    ax.set_xlabel("Median exact-junction effect")
    ax.set_title("Locus-specific breadth", loc="left", fontsize=6.4, fontweight="bold")
    ax.legend(
        frameon=False,
        fontsize=4.5,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        ncol=1,
        handletextpad=0.3,
        borderaxespad=0,
    )
    ax.text(
        0.02,
        0.02,
        "labels: positive / estimable donors",
        transform=ax.transAxes,
        fontsize=4.4,
        color="#666666",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery", required=True, type=Path)
    parser.add_argument("--sealed", required=True, type=Path)
    parser.add_argument("--short-read", required=True, type=Path)
    parser.add_argument("--extension-summary", required=True, type=Path)
    parser.add_argument("--state-sensitivity", required=True, type=Path)
    parser.add_argument("--font", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", default="v0.2")
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)

    setup_style(args.font)
    effects = load_effects(args.discovery, args.sealed, args.short_read)
    profiles = load_cell_profiles(args.extension_summary)
    robustness = load_robustness(args.state_sensitivity)
    breadth = load_breadth(args.extension_summary, effects)

    fig = plt.figure(figsize=(FIGURE_SIZE_MM[0] / 25.4, FIGURE_SIZE_MM[1] / 25.4))
    grid = fig.add_gridspec(
        4,
        12,
        height_ratios=[0.72, 1.0, 0.83, 1.15],
        left=0.065,
        right=0.985,
        top=0.985,
        bottom=0.105,
        hspace=0.50,
        wspace=0.95,
    )
    ax_a = fig.add_subplot(grid[0, :8])
    ax_b = fig.add_subplot(grid[0, 8:])
    ax_c = fig.add_subplot(grid[1, :4])
    ax_d = fig.add_subplot(grid[1, 4:8])
    ax_e = fig.add_subplot(grid[1, 8:])
    ax_f = fig.add_subplot(grid[2, :4])
    ax_g = fig.add_subplot(grid[2, 4:])
    ax_h = fig.add_subplot(grid[3, :8])
    ax_i = fig.add_subplot(grid[3, 8:])

    draw_panel_a(ax_a)
    draw_panel_b(ax_b)
    draw_donor_panel(
        ax_c,
        effects,
        "Discovery LR",
        "#0173B2",
        "c",
        "Discovery long reads",
        r"exact q=5.49×$10^{-4}$",
        True,
    )
    draw_donor_panel(
        ax_d,
        effects,
        "Sealed LR",
        "#56B4E9",
        "d",
        "Held-back long reads",
        "gene P=0.125",
        False,
    )
    draw_donor_panel(
        ax_e,
        effects,
        "Independent SR",
        "#DE8F05",
        "e",
        "Independent short reads",
        r"sign P=1.91×$10^{-6}$",
        False,
    )
    medians = draw_panel_f(ax_f, effects)
    draw_panel_g(ax_g, robustness)
    draw_panel_h(ax_h, profiles)
    draw_panel_i(ax_i, breadth)

    stem = args.output_dir / f"Figure2_biological_application_rich_{args.version}"
    fig.savefig(f"{stem}.pdf", transparent=True)
    fig.savefig(f"{stem}.svg", transparent=True)
    fig.savefig(f"{stem}.png", dpi=450, transparent=True)
    fig.savefig(f"{stem}.tiff", dpi=600, transparent=True)
    plt.close(fig)

    effects.to_csv(args.output_dir / "Figure2_source_donor_effects.tsv", sep="\t", index=False)
    medians.to_csv(args.output_dir / "Figure2_source_cohort_medians.tsv", sep="\t", index=False)
    profiles.to_csv(args.output_dir / "Figure2_source_celltype_profiles.tsv", sep="\t", index=False)
    robustness.to_csv(args.output_dir / "Figure2_source_state_sensitivity.tsv", sep="\t", index=False)
    breadth.to_csv(args.output_dir / "Figure2_source_breadth_boundary.tsv", sep="\t", index=False)
    audit = {
        "schema_version": f"scTHREAD-rich-biology-figure-{args.version}",
        "panel_ids": list(PANEL_IDS),
        "figure_size_mm": list(FIGURE_SIZE_MM),
        "core_genes": list(GENES),
        "claim_boundary": (
            "Exact-junction usage, descriptive state sensitivity and locus-specific "
            "breadth are shown; no full-length short-read linkage, shared regulation "
            "or mechanism is claimed."
        ),
    }
    (args.output_dir / f"Figure2_audit_{args.version}.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
