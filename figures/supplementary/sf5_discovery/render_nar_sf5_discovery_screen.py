#!/usr/bin/env python3
"""Render the complete nine-locus discovery screen for Supplementary Figure S5."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "style"))
import nar_style as S


GENES = ("CCND3", "CUTA", "PSMB8", "RPL41", "SMAP2", "EMP3", "EVI2B", "IRF1-AS1", "S100A6")
SELECTED = {"CCND3": S.BLUE, "CUTA": S.TEAL, "PSMB8": "#815A92"}
DISCOVERY_ONLY = {"RPL41", "SMAP2"}
OTHER = S.GREY
FIGURE_SIZE_MM = (183, 118)


def setup_style(font_path: Path) -> FontProperties:
    fontManager.addfont(str(font_path))
    font = FontProperties(fname=str(font_path))
    family = font.get_name()
    mpl.rcParams.update(
        {
            "font.family": family,
            "font.sans-serif": [family],
            "font.size": 7.0,
            "axes.labelsize": 7.0,
            "axes.titlesize": 7.5,
            "xtick.labelsize": 6.0,
            "ytick.labelsize": 6.0,
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.transparent": True,
            "figure.facecolor": "none",
            "axes.facecolor": "none",
            "axes.edgecolor": S.INK,
            "text.color": S.INK,
            "axes.labelcolor": S.INK,
            "xtick.color": S.INK,
            "ytick.color": S.INK,
        }
    )
    return font


def panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.08, label.lower(), transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="right",
            color=S.INK, clip_on=False)


def load_inputs(effects_path: Path, test_path: Path) -> tuple[pd.DataFrame, dict]:
    effects = pd.read_csv(effects_path, sep="\t")
    if {"cohort", "donor", "effect"}.issubset(effects.columns):
        effects = effects.loc[effects["cohort"].eq("Discovery LR")].copy()
    elif {"biological_unit", "monocyte_minus_t_cell_usage"}.issubset(effects.columns):
        effects = effects.rename(
            columns={"biological_unit": "donor", "monocyte_minus_t_cell_usage": "effect"}
        ).copy()
    else:
        raise ValueError("Discovery source lacks the required donor/effect columns")
    effects = effects.loc[effects["gene"].isin(GENES)].copy()
    effects["cohort"] = "Discovery LR"
    if set(effects["gene"]) != set(GENES):
        raise ValueError("Discovery source does not contain the frozen nine-locus panel")
    counts = effects.groupby("gene").size()
    if not counts.eq(13).all():
        raise ValueError(f"Expected 13 discovery donors per locus; observed {counts.to_dict()}")
    audit = json.loads(test_path.read_text(encoding="utf-8"))
    if audit["biological_units"] != 13 or audit["estimable_genes"] != 9:
        raise ValueError("Discovery exact-test contract changed")
    if audit["observed_global_genes"] != 6:
        raise ValueError("Unexpected global discovery summary")
    return effects, audit


def gene_color(gene: str) -> str:
    if gene in SELECTED:
        return SELECTED[gene]
    if gene in DISCOVERY_ONLY:
        return "#5B5B5B"
    return OTHER


def draw_panel_a(ax: mpl.axes.Axes, effects: pd.DataFrame, audit: dict) -> None:
    panel_label(ax, "a")
    rng = np.random.default_rng(20260828)
    y_positions = np.arange(len(GENES))[::-1]
    ax.axvline(0, color="#777777", lw=0.55, ls=(0, (2, 2)), zorder=0)
    for y, gene in zip(y_positions, GENES):
        values = effects.loc[effects["gene"].eq(gene), "effect"].to_numpy()
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        color = gene_color(gene)
        ax.scatter(values, np.full(len(values), y) + jitter, s=11, color=color,
                   edgecolor="white", linewidth=0.3, alpha=0.88, zorder=2)
        median = float(np.median(values))
        ax.plot([median - 0.025, median + 0.025], [y, y], color="#111111", lw=1.15, zorder=3)
        entry = audit["genes"][gene]
        q = entry["benjamini_hochberg_q"]
        q_text = f"q={q:.3g}" if q >= 0.001 else f"q={q:.2g}"
        ax.text(1.08, y, f"{int(entry['positive_units'])}/13  {q_text}",
                va="center", ha="right", fontsize=5.0, color=color)
    ax.set_yticks(y_positions, GENES, fontweight="bold")
    for tick in ax.get_yticklabels():
        tick.set_color(gene_color(tick.get_text()))
    ax.set_xlim(-0.12, 1.18)
    ax.set_ylim(-0.7, len(GENES) - 0.25)
    ax.set_xlabel("Discovery donor effect\n(monocyte − T exact-junction usage)")
    ax.set_title("Donor-level effects", loc="left", fontweight="bold", pad=4)
    ax.text(1.0, 1.01, "positive / 13 donors   BH q", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.0, color="#555555")


def draw_panel_b(ax: mpl.axes.Axes, effects: pd.DataFrame) -> None:
    panel_label(ax, "b")
    matrix = effects.pivot(index="gene", columns="donor", values="effect").loc[list(GENES)]
    # Keep donors in source order, which is the order used by the frozen table.
    matrix = matrix.loc[:, effects["donor"].drop_duplicates().tolist()]
    image = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1) - 0.5,
        np.arange(matrix.shape[0] + 1) - 0.5,
        matrix.to_numpy(),
        cmap=mpl.colormaps["RdBu_r"],
        norm=TwoSlopeNorm(vmin=-0.15, vcenter=0.0, vmax=0.9),
        shading="flat",
        rasterized=False,
    )
    ax.set_yticks(range(len(GENES)), GENES, fontweight="bold")
    for tick in ax.get_yticklabels():
        tick.set_color(gene_color(tick.get_text()))
    donor_labels = [d.replace("mm", "m").replace("omb", "o").replace("orb", "r").replace("prom", "p")
                    for d in matrix.columns]
    ax.set_xticks(range(len(donor_labels)), donor_labels, rotation=62, ha="left")
    ax.set_xlim(-0.5, matrix.shape[1] - 0.5)
    ax.set_ylim(matrix.shape[0] - 0.5, -0.5)
    ax.tick_params(length=0)
    ax.set_title("Effect matrix across donors", loc="left", fontweight="bold", pad=4)
    ax.set_xlabel("Biological unit (13 donors)")
    colorbar = ax.figure.colorbar(image, ax=ax, orientation="horizontal", fraction=0.08, pad=0.27)
    colorbar.solids.set_rasterized(False)
    colorbar.outline.set_rasterized(False)
    colorbar.set_label("Monocyte − T effect", fontsize=5.4)
    colorbar.set_ticks([-0.1, 0, 0.4, 0.8])
    colorbar.ax.tick_params(labelsize=5.0, width=0.4, length=2)


def box(ax, xy, width, height, title, subtitle, facecolor, edgecolor="#555555") -> None:
    x, y = xy
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.015",
        facecolor=facecolor, edgecolor=edgecolor, lw=0.6,
    ))
    ax.text(x + width / 2, y + height * 0.62, title, ha="center", va="center",
            fontsize=7.0, fontweight="bold")
    ax.text(x + width / 2, y + height * 0.30, subtitle, ha="center", va="center",
            fontsize=5.2, color="#555555")


def draw_panel_c(ax: mpl.axes.Axes, audit: dict) -> None:
    panel_label(ax, "c")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.98, "Selection path fixed before cross-platform follow-up", ha="left", va="top",
            fontsize=7.5, fontweight="bold", color=S.INK)
    steps = (
        (0.03, "11", "nominated\ncandidates", "#E8F1F6"),
        (0.27, "9", "estimable\nloci", "#E8F1F6"),
        (0.51, "5", "BH q < 0.05", "#E5F3EC"),
        (0.75, "3", "independent\nshort-read loci", "#FFF0D9"),
    )
    for idx, (x, title, subtitle, color) in enumerate(steps):
        box(ax, (x, 0.49), 0.17, 0.27, title, subtitle, color)
        if idx < len(steps) - 1:
            ax.add_patch(FancyArrowPatch(
                (x + 0.18, 0.625), (steps[idx + 1][0] - 0.015, 0.625),
                arrowstyle="->", mutation_scale=8, lw=0.7, color="#555555",
            ))
    ax.text(0.02, 0.30, "Excluded before testing: SSU72 and TRAPPC1, both lacking direct observation of both carrier directions.",
            ha="left", va="center", fontsize=5.1, color="#555555")
    ax.text(0.02, 0.16, "Five corrected discovery loci: CCND3, CUTA, PSMB8, RPL41 and SMAP2.",
            ha="left", va="center", fontsize=5.3)
    ax.text(0.02, 0.07, "RPL41 and SMAP2 remain long-read observations; CCND3, CUTA and PSMB8 form the cross-platform panel.",
            ha="left", va="center", fontsize=5.1, color="#555555")
    ax.text(0.98, 0.01, "Two-sided exact sign-flip test; 8,192 assignments; BH correction over nine loci",
            ha="right", va="bottom", fontsize=5.0, color="#666666")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effects", required=True, type=Path)
    parser.add_argument("--test-json", required=True, type=Path)
    parser.add_argument("--font", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True)
    setup_style(args.font)
    effects, audit = load_inputs(args.effects, args.test_json)

    fig = plt.figure(figsize=(FIGURE_SIZE_MM[0] * S.MM, FIGURE_SIZE_MM[1] * S.MM))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 0.56],
        width_ratios=[1.12, 1.0],
        left=0.085,
        right=0.985,
        top=0.895,
        bottom=0.115,
        hspace=0.62,
        wspace=0.42,
    )
    draw_panel_a(fig.add_subplot(grid[0, 0]), effects, audit)
    draw_panel_b(fig.add_subplot(grid[0, 1]), effects)
    draw_panel_c(fig.add_subplot(grid[1, :]), audit)
    fig.text(0.085, 0.025,
             "All points are donor-level effects from the 13-donor discovery long-read cohort; values are not pooled across donors.",
             ha="left", va="bottom", fontsize=4.9, color="#666666")

    fig.text(0.085, 0.965,
             "Supplementary Figure 5 · Complete nine-locus discovery screen",
             ha="left", va="top", fontsize=9.5, fontweight="bold", color=S.INK)

    stem = args.output_dir / "NAR_SF5_discovery_screen_20260828"
    with mpl.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0}):
        fig.savefig(f"{stem}.pdf", transparent=True)
        fig.savefig(f"{stem}.svg", transparent=True)
        fig.savefig(f"{stem}.png", dpi=450, transparent=True)
        fig.savefig(f"{stem}.tiff", dpi=600, transparent=True)
    plt.close(fig)

    effects.to_csv(args.output_dir / "NAR_SF5_source_discovery_effects.tsv", sep="\t", index=False)
    (args.output_dir / "NAR_SF5_source_exact_test.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    contract = {
        "figure": "Supplementary Figure S5",
        "figure_size_mm": list(FIGURE_SIZE_MM),
        "panels": {"a": "nine-locus donor effects and BH q-values", "b": "13-donor effect matrix", "c": "fixed selection path"},
        "source_effect_rows": int(len(effects)),
        "source_genes": list(GENES),
        "biological_units": 13,
        "claim_boundary": "Complete discovery-screen audit; the figure does not establish mechanism or causal regulation.",
    }
    (args.output_dir / "NAR_SF5_audit_20260828.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
