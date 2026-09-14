#!/usr/bin/env python3
"""Render a publication-grade CD45 read-versus-UMI sensitivity figure.

Figure contract
---------------
Claim: the CD45 exon-3 splice-choice pattern is robust to UMI deduplication
within the two tested studies, while the primary published measurement remains
explicitly read-level.
Evidence: eight study-by-cell-type strata, 31 successfully annotated runs,
312,919 informative reads and 114,610 UMI class counts from the frozen TSV.
Archetype: quantitative grid with paired validation and effect-size boundary.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.font_manager import fontManager

import argparse
import os
ARIAL = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))
import sys
REPO = Path(__file__).resolve().parents[3]
ROOT = REPO
INPUT = REPO / "evidence_layers/cd45_sensitivity/cd45_umi_sensitivity_20260826.tsv"
OUT = S.OUTDIR / "NAR_SF7"
sys.path.insert(0, str(REPO / "figures/style"))
import nar_style as S


def panel_label(ax, label: str) -> None:
    ax.text(-0.08, 1.08, label.lower(), transform=ax.transAxes,
            fontsize=8, fontweight="bold", va="top", ha="right",
            color=S.INK, clip_on=False)

if ARIAL.exists():
    fontManager.addfont(str(ARIAL))
    FONT = "Arial"
else:
    FONT = "DejaVu Sans"

mpl.rcParams.update(
    {
        "font.family": FONT,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "axes.facecolor": "none",
        "axes.edgecolor": S.INK,
        "text.color": S.INK,
        "axes.labelcolor": S.INK,
        "xtick.color": S.INK,
        "ytick.color": S.INK,
    }
)


STUDY_COLORS = {"GSE276974": S.BLUE, "GSE307660": S.CORAL}
CELL_COLORS = {
    "B cell": S.BLUE,
    "NK": "#815A92",
    "Monocyte": S.TEAL,
    "T cell": S.CORAL,
}


def main() -> None:
    df = pd.read_csv(INPUT, sep="\t")
    required = {
        "study", "cell_type", "eligible_runs", "informative_reads",
        "informative_umi_molecules", "ro_read_fraction", "ro_umi_fraction",
        "absolute_delta_percentage_points",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if len(df) != 8:
        raise ValueError(f"Expected 8 study-by-cell-type strata, found {len(df)}")

    df = df.copy()
    df["read_pct"] = df["ro_read_fraction"] * 100
    df["umi_pct"] = df["ro_umi_fraction"] * 100
    df["label"] = df["study"] + " · " + df["cell_type"]
    max_delta = float(df["absolute_delta_percentage_points"].max())
    # These totals are the independently audited run-level summary reported in
    # cd45_umi_sensitivity_20260826.json and S15.  The eight displayed strata
    # overlap runs across cell types, so summing their row counts would double
    # count reads and UMI classes.
    total_reads = 312_919
    total_umi = 114_610
    total_runs = 31

    fig = plt.figure(figsize=(183 * S.MM, 100 * S.MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.08, 1.0, 1.15], wspace=0.42)
    ax_a, ax_b, ax_c = [fig.add_subplot(gs[0, i]) for i in range(3)]

    fig.text(
        0.085, 0.965,
        "Supplementary Figure 7 · CD45 read-versus-UMI sensitivity",
        ha="left", va="top", fontsize=9.5, fontweight="bold", color=S.INK,
    )
    panel_label(ax_a, "a")
    panel_label(ax_b, "b")
    panel_label(ax_c, "c")
    fig.text(
        0.085, 0.915,
        "Orthogonal read-level and UMI-class measurements preserve the within-study usage pattern",
        ha="left", va="top", fontsize=5.7, color="#526577",
    )

    # Panel A: paired read/UMI slopes. Offset the two studies to keep every
    # observed stratum visible while retaining a common percentage axis.
    study_offsets = {"GSE276974": -0.07, "GSE307660": 0.07}
    x_base = {"read": 0.0, "umi": 1.0}
    for _, row in df.iterrows():
        off = study_offsets[row["study"]]
        x = [x_base["read"] + off, x_base["umi"] + off]
        y = [row["read_pct"], row["umi_pct"]]
        ax_a.plot(x, y, color=CELL_COLORS[row["cell_type"]], lw=1.1, alpha=0.82, zorder=2)
        ax_a.scatter(x, y, s=22, color=CELL_COLORS[row["cell_type"]], edgecolor="white", lw=0.45, zorder=3)
    ax_a.set_xticks([0, 1], ["Read\nlevel", "UMI\nclass"])
    ax_a.set_ylabel("CD45RO fraction (%)", fontsize=6.7)
    ax_a.set_ylim(-1, 44)
    ax_a.set_title("Paired estimates", loc="left", fontsize=7.5, fontweight="bold", pad=4)
    ax_a.grid(axis="y", color="#E6EBF0", lw=0.45)
    ax_a.tick_params(labelsize=6.1, length=2.5)
    ax_a.text(0.02, 0.98, "lines: cell type", transform=ax_a.transAxes, fontsize=5.3, color="#526577", va="top")

    # Panel B: cross-measurement agreement and study-specific rank preservation.
    for study, g in df.groupby("study", sort=True):
        ax_b.scatter(
            g["read_pct"], g["umi_pct"], s=34, color=STUDY_COLORS[study],
            edgecolor="white", lw=0.55, label=study, zorder=3,
        )
        for _, row in g.iterrows():
            ax_b.text(row["read_pct"] + 0.35, row["umi_pct"] + 0.15, row["cell_type"].replace(" cell", ""),
                      fontsize=5.0, color=STUDY_COLORS[study], va="center")
    lim = 44
    ax_b.plot([0, lim], [0, lim], ls="--", lw=0.7, color="#8A99A8", zorder=1)
    ax_b.set_xlim(-1, lim)
    ax_b.set_ylim(-1, lim)
    ax_b.set_xlabel("Read-level fraction (%)", fontsize=6.7)
    ax_b.set_ylabel("UMI-class fraction (%)", fontsize=6.7)
    ax_b.set_title("Agreement", loc="left", fontsize=7.5, fontweight="bold", pad=4)
    ax_b.legend(frameon=False, fontsize=5.2, loc="upper left", handletextpad=0.3, borderpad=0.2)
    ax_b.text(0.04, 0.08, "within-study Spearman ρ = 1.00\nfor both cohorts", transform=ax_b.transAxes,
              fontsize=5.4, color="#344454", va="bottom")
    ax_b.tick_params(labelsize=6.1, length=2.5)

    # Panel C: absolute deltas with the real sample-size context.
    plot_df = df.sort_values(["study", "absolute_delta_percentage_points"]).reset_index(drop=True)
    ys = list(range(len(plot_df)))
    for y, (_, row) in zip(ys, plot_df.iterrows()):
        c = STUDY_COLORS[row["study"]]
        ax_c.hlines(y, 0, row["absolute_delta_percentage_points"], color=c, lw=1.3, alpha=0.72)
        ax_c.scatter(row["absolute_delta_percentage_points"], y, s=28, color=c, edgecolor="white", lw=0.5, zorder=3)
    ax_c.set_yticks(ys, [f"{r.study[-6:]} · {r.cell_type}" for r in plot_df.itertuples()], fontsize=5.0)
    ax_c.set_xlabel("Absolute change (percentage points)", fontsize=6.7)
    ax_c.set_xlim(0, 2.55)
    ax_c.set_title("Absolute sensitivity", loc="left", fontsize=7.5, fontweight="bold", pad=4)
    ax_c.grid(axis="x", color="#E6EBF0", lw=0.45)
    ax_c.tick_params(axis="x", labelsize=6.1, length=2.5)
    ax_c.text(0.98, 0.97, f"max = {max_delta:.2f} pp", transform=ax_c.transAxes,
              ha="right", va="top", fontsize=5.6, color="#344454", fontweight="bold")

    fig.text(
        0.02, 0.015,
        f"31 annotated runs · {total_reads:,} informative reads · {total_umi:,} informative UMI classes · conflicting CB/UB keys excluded",
        ha="left", va="bottom", fontsize=5.8, color="#526577",
    )
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.19, top=0.84)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0}):
        fig.savefig(f"{OUT}.svg", bbox_inches=None, pad_inches=0, transparent=True)
        fig.savefig(f"{OUT}.pdf", bbox_inches=None, pad_inches=0, transparent=True)
        fig.savefig(f"{OUT}.png", dpi=450, bbox_inches=None, pad_inches=0, transparent=True)
        fig.savefig(f"{OUT}.tiff", dpi=600, bbox_inches=None, pad_inches=0, transparent=True)
    plt.close(fig)
    print(OUT.with_suffix(".svg"))
    print(OUT.with_suffix(".pdf"))
    print(OUT.with_suffix(".tiff"))


if __name__ == "__main__":
    main()
