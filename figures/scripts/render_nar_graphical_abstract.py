#!/usr/bin/env python3
"""Render a simplified NAR-compliant scTHREAD graphical abstract."""

from __future__ import annotations

import os

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STEM = ROOT / "figures/scTHREAD_graphical_abstract_NAR"
WIDTH_MM = 150.0
HEIGHT_MM = 60.0
MM = 1 / 25.4
FONT_PATH = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))

NAVY = "#17324D"
BLUE = "#0173B2"
ORANGE = "#DE8F05"
GREEN = "#029E73"
PALE_BLUE = "#E7F2F8"
PALE_ORANGE = "#FBEED9"
PALE_GREEN = "#E3F3EC"
GREY = "#5C6873"
WHITE = "#FFFFFF"


def configure() -> None:
    fontManager.addfont(str(FONT_PATH))
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 600,
            "savefig.transparent": True,
            "figure.facecolor": "none",
            "axes.facecolor": "none",
        }
    )


def card(ax, x: float, color: str, fill: str, title: str) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, 7),
            40,
            34,
            boxstyle="round,pad=0.8,rounding_size=2.4",
            facecolor=fill,
            edgecolor=color,
            linewidth=1.6,
        )
    )
    ax.text(
        x + 20,
        35.5,
        title,
        ha="center",
        va="center",
        color=NAVY,
        fontsize=12,
        fontweight="bold",
    )


def source_icon(ax, x: float) -> None:
    for cx, cy, radius, color in (
        (x + 12, 22, 4.8, BLUE),
        (x + 20, 17, 4.8, ORANGE),
        (x + 28, 22, 4.8, GREEN),
    ):
        ax.add_patch(Circle((cx, cy), radius, facecolor=color, edgecolor=WHITE, linewidth=1.2))
    ax.text(x + 20, 12.0, "Human + mouse\nONT + PacBio", ha="center", va="center", fontsize=12, color=NAVY)


def evidence_icon(ax, x: float) -> None:
    labels = (("Isoforms", BLUE), ("Junctions", ORANGE), ("Poly(A)", GREEN))
    for index, (label, color) in enumerate(labels):
        y = 26 - index * 7
        ax.plot([x + 8, x + 14, x + 20], [y, y + 2.2, y], color=color, linewidth=2.4)
        ax.text(x + 24, y + 0.5, label, ha="left", va="center", fontsize=12, color=NAVY)


def access_icon(ax, x: float) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x + 9, 16),
            22,
            14,
            boxstyle="round,pad=0.4,rounding_size=1.2",
            facecolor=WHITE,
            edgecolor=GREEN,
            linewidth=1.5,
        )
    )
    ax.plot([x + 12, x + 28], [26, 26], color=GREEN, linewidth=1.5)
    ax.text(x + 20, 21.5, "PTPRC", ha="center", va="center", fontsize=12, color=NAVY, fontweight="bold")
    ax.text(x + 20, 12.0, "Search · inspect\nAPI · download", ha="center", va="center", fontsize=12, color=NAVY)


def build(stem: Path) -> None:
    configure()
    fig = plt.figure(figsize=(WIDTH_MM * MM, HEIGHT_MM * MM), facecolor="none")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH_MM)
    ax.set_ylim(0, HEIGHT_MM)
    ax.axis("off")

    ax.text(
        75,
        54,
        "scTHREAD",
        ha="center",
        va="center",
        fontsize=16,
        color=NAVY,
        fontweight="bold",
    )
    ax.text(
        75,
        47.5,
        "Provenance-linked single-cell long-read transcript evidence",
        ha="center",
        va="center",
        fontsize=12,
        color=GREY,
    )

    card(ax, 4, BLUE, PALE_BLUE, "Multi-study data")
    card(ax, 55, ORANGE, PALE_ORANGE, "Evidence layers")
    card(ax, 106, GREEN, PALE_GREEN, "Open retrieval")
    source_icon(ax, 4)
    evidence_icon(ax, 55)
    access_icon(ax, 106)

    for start, end in ((44.5, 54.5), (95.5, 105.5)):
        ax.add_patch(
            FancyArrowPatch(
                (start, 24),
                (end, 24),
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.8,
                color=NAVY,
            )
        )

    ax.text(
        75,
        2.5,
        "https://scthread.ai4sc.ac.cn",
        ha="center",
        va="bottom",
        fontsize=12,
        color=BLUE,
        fontweight="bold",
    )

    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), transparent=True)
    fig.savefig(stem.with_suffix(".svg"), transparent=True)
    fig.savefig(stem.with_suffix(".png"), dpi=600, transparent=True)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, transparent=True)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, default=DEFAULT_STEM)
    args = parser.parse_args()
    build(args.stem.resolve())


if __name__ == "__main__":
    main()
