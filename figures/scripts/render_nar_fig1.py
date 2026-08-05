#!/usr/bin/env python3
"""NAR Fig. 1 — Database content, scope and architecture (Database Issue, not GB)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

OUT = "NAR_Fig1"


def _round_box(ax, x, y, w, h, text, fc, ec, fs=6.2, fw="normal", tc=S.INK):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.012,rounding_size=0.04",
            facecolor=fc, edgecolor=ec, linewidth=0.8, mutation_aspect=0.6,
            transform=ax.transAxes, clip_on=False,
        )
    )
    ax.text(
        x + w / 2, y + h / 2, text,
        transform=ax.transAxes, ha="center", va="center",
        fontsize=fs, fontweight=fw, color=tc, linespacing=1.25,
    )


def panel_a_workflow(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    S.panel_label(ax, "a", x=-0.02, y=1.05)
    ax.set_title("Architecture: public long reads → multi-layer database", pad=6)

    boxes = [
        (0.02, 0.28, 0.18, 0.48, "Public sc\nlong-read\ndatasets", S.CREAM, S.GREY),
        (0.27, 0.28, 0.20, 0.48, "Uniform\nreprocessing\n(IsoQuant)", S.SOFT, S.TEAL),
        (0.54, 0.18, 0.20, 0.68,
         "Evidence layers\n• isoform\n• junction\n• poly(A)\n• allele", "#EAF3F0", S.TEAL),
        (0.80, 0.28, 0.18, 0.48, "Web portal\n+ REST API\nscthread.ai4sc.ac.cn",
         "#F8EFE8", S.CORAL),
    ]
    for x, y, w, h, t, fc, ec in boxes:
        _round_box(ax, x, y, w, h, t, fc, ec, fs=6.0, fw="bold" if "Evidence" in t or "Web" in t else "normal")

    for x0, x1 in ((0.20, 0.27), (0.47, 0.54), (0.74, 0.80)):
        ax.annotate(
            "", xy=(x1, 0.52), xytext=(x0, 0.52),
            xycoords=ax.transAxes, textcoords=ax.transAxes,
            arrowprops=dict(arrowstyle="-|>", color=S.SLATE, lw=0.9, mutation_scale=8),
        )
    ax.text(
        0.5, 0.06,
        "Registration-free · gene-centric multi-evidence browser · versioned tables",
        transform=ax.transAxes, ha="center", va="center", fontsize=5.8, color=S.SLATE,
    )


def panel_b_dual_scope(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    S.panel_label(ax, "b", x=-0.02, y=1.05)
    ax.set_title("Two scopes (do not mix denominators)", pad=6)

    # Catalog card
    _round_box(
        ax, 0.04, 0.12, 0.44, 0.76,
        "",
        "#F7F3EC", S.GOLD, fs=1,
    )
    ax.text(0.26, 0.78, "Public catalog", transform=ax.transAxes,
            ha="center", fontsize=7, fontweight="bold", color=S.GOLD)
    ax.text(0.26, 0.58, "≈450 runs\n≈3.0M cells\n>200k isoforms",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=7.5, fontweight="bold", color=S.INK, linespacing=1.35)
    ax.text(0.26, 0.22, "Portal About / API\nresource scale",
            transform=ax.transAxes, ha="center", fontsize=5.6, color=S.SLATE)

    # Cohort card
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())
    _round_box(
        ax, 0.52, 0.12, 0.44, 0.76,
        "",
        S.SOFT, S.TEAL, fs=1,
    )
    ax.text(0.74, 0.78, "Uniform cohort", transform=ax.transAxes,
            ha="center", fontsize=7, fontweight="bold", color=S.TEAL)
    ax.text(
        0.74, 0.58,
        f"{ch['n_studies']} studies\n{ch['n_runs']} runs\n{ch['n_annotated_cells']:,} cells",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=7.5, fontweight="bold", color=S.INK, linespacing=1.35,
    )
    ax.text(0.74, 0.22, "All reprocessed\nfeatured tables",
            transform=ax.transAxes, ha="center", fontsize=5.6, color=S.SLATE)


def panel_c_cells(ax):
    cc = pd.read_csv(S.FIGDATA / "cohort_cellcount.tsv", sep="\t").sort_values("n_cells")
    S.panel_label(ax, "c")
    colors = [S.TEAL if n >= 20_000 else S.BLUE if n >= 5_000 else S.GREY for n in cc.n_cells]
    ax.barh(cc.gse, cc.n_cells / 1000, color=colors, edgecolor="white", linewidth=0.4, height=0.72)
    ax.set_xlabel("Annotated cells (×10³)")
    ax.set_title("Uniform cohort: annotated cells per study", pad=4)
    ax.set_xlim(0, cc.n_cells.max() / 1000 * 1.12)
    total = cc.n_cells.sum()
    ax.text(
        0.98, 0.02, f"sum = {total:,}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=S.SLATE,
    )
    S.style_ax(ax)


def panel_d_classes(ax):
    cls = pd.read_csv(S.FIGDATA / "fig1b_classification_all.tsv", sep="\t")
    # aggregate by system
    order = ["FSM", "ISM", "NIC", "NNC"]
    systems = list(cls.groupby("system")["novel_frac"].mean().sort_values(ascending=False).index)
    # rebuild stacked fractions per system
    rows = []
    for sys in systems:
        sub = cls[cls.system == sys]
        # average class fractions across studies in system
        means = {c: sub[c].mean() for c in order}
        s = sum(means.values())
        means = {c: means[c] / s for c in order}
        rows.append((sys, means))

    y = np.arange(len(rows))
    left = np.zeros(len(rows))
    for c in order:
        vals = np.array([r[1][c] for r in rows])
        ax.barh(y, vals, left=left, color=S.CLASS[c], edgecolor="white",
                linewidth=0.3, height=0.7, label=c)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of spliced molecules")
    ax.set_title("Content composition by biological system", pad=4)
    ax.legend(loc="lower right", ncol=4, fontsize=5.5, handlelength=0.9,
              columnspacing=0.6, borderaxespad=0.2)
    ax.axvline(0.5, color=S.GREY, lw=0.4, ls=":", zorder=0)
    S.panel_label(ax, "d")
    S.style_ax(ax)


def main():
    fig = plt.figure(figsize=(183 * S.MM, 145 * S.MM))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 1.15],
        hspace=0.42, wspace=0.32,
        left=0.08, right=0.98, top=0.91, bottom=0.07,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    panel_a_workflow(ax_a)
    panel_b_dual_scope(ax_b)
    panel_c_cells(ax_c)
    panel_d_classes(ax_d)

    fig.suptitle(
        "scTHREAD database content and architecture",
        x=0.07, y=0.99, ha="left", fontsize=9, fontweight="bold", color=S.INK,
    )
    path = S.save(fig, OUT)
    print("wrote", path)


if __name__ == "__main__":
    main()
