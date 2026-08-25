#!/usr/bin/env python3
"""NAR Database Issue–style Fig.1 (content + construction).

Does NOT overwrite NAR_Fig1.* — writes NAR_Fig1_narstyle.{pdf,png}.

Design (aligned with CIRCpedia / PolyA_DB / EVmiRNA Fig.1 habits):
  a  Construction schematic + live catalog KPIs
  b  Species × platform composition
  c  Biological systems by isoquant_cells
  d  Dual scope (catalog vs labeled analysis subset)

All catalog numbers from sample_registry isoquant_status=done.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S
from render_nar_bio import catalog_headline, load_catalog_bio

INK, SLATE, TEAL, CORAL, GOLD, BLUE, GREY, CREAM = (
    S.INK, S.SLATE, S.TEAL, S.CORAL, S.GOLD, S.BLUE, S.GREY, S.CREAM
)
GREYL = "#D9DCE0"
PURPLE = "#6B5B95"


def tag(ax, letter, x=-0.02, y=1.06):
    S.panel_label(ax, letter, x=x, y=y)


def title(ax, s, pad=3.0):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=7.4, color=INK, pad=pad)


def card(ax, x, y, w, h, fc=CREAM, ec=GREYL, lw=0.9, rs=0.035):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.008,rounding_size={rs}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=ax.transAxes, clip_on=False,
        zorder=1,
    ))


def render():
    cat = catalog_headline()
    studies = load_catalog_bio()
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())

    # species for cells: fix ? as unknown
    sp = studies.copy()
    sp["sp_n"] = sp.species.fillna("?").str.lower().map(
        lambda s: "human" if "human" in s else ("mouse" if "mouse" in s or s == "mm" else "other")
    )
    # re-attribute cells for ? species by platform-heavy ONT nonstandard as other
    cells_by_sp = sp.groupby("sp_n")["n_cells"].sum()
    runs_by_sp = sp.groupby("sp_n")["n_runs"].sum()
    cells_by_plat = sp.groupby(sp.platform.fillna("NA"))["n_cells"].sum()
    runs_by_plat = sp.groupby(sp.platform.fillna("NA"))["n_runs"].sum()

    # systems: merge tiny buckets
    sys_merge = {
        "Disease/clinical": "Blood/marrow",
        "Splicing/APA": "Blood/marrow",
        "Mouse": "Other",
    }
    st = studies.copy()
    st["sys2"] = st.system.map(lambda s: sys_merge.get(s, s))
    # keep Smart-seq2 separate only if non-trivial; else fold into Other
    if st.loc[st.sys2 == "Smart-seq2", "n_cells"].sum() < 500:
        st.loc[st.sys2 == "Smart-seq2", "sys2"] = "Other"
    sys_order = ["Blood/marrow", "Brain", "Heart", "Cancer", "Differentiation", "Benchmark", "Other"]
    sys_order = [s for s in sys_order if s in set(st.sys2)]
    for s in sorted(set(st.sys2) - set(sys_order)):
        sys_order.append(s)
    colors_sys = {
        "Blood/marrow": TEAL, "Brain": PURPLE, "Heart": CORAL, "Cancer": GOLD,
        "Differentiation": BLUE, "Benchmark": GREY, "Other": "#A0A4AB",
    }

    fig = plt.figure(figsize=(183 * S.MM, 145 * S.MM), facecolor="none")
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.05, 1.0],
        hspace=0.38, wspace=0.28,
        left=0.08, right=0.98, top=0.90, bottom=0.07,
    )

    # ================================================================= a schematic + KPI
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Construction of the scTHREAD resource")
    ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    steps = [
        (0.02, "Public sc\nlong-read\ndatasets", TEAL),
        (0.27, "Uniform\nIsoQuant\nreprocessing", BLUE),
        (0.52, "Four evidence\nlayers\niso·junc·PAS·ASE", GOLD),
        (0.77, "Portal + API\nfree access\nscthread.ai4sc.ac.cn", CORAL),
    ]
    for i, (x0, lab, col) in enumerate(steps):
        card(ax, x0, 0.52, 0.21, 0.40, fc="#FFFFFF", ec=col, lw=1.15, rs=0.04)
        ax.text(x0 + 0.105, 0.72, lab, transform=ax.transAxes, ha="center", va="center",
                fontsize=5.6, fontweight="bold", color=col, linespacing=1.25)
        if i < 3:
            ax.annotate(
                "", xy=(x0 + 0.225, 0.72), xytext=(x0 + 0.205, 0.72),
                xycoords=ax.transAxes, textcoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.0,
                                mutation_scale=9),
            )

    # KPI strip under schematic
    kpis = [
        (f"{cat['n_runs']}", "runs"),
        (f"{cat['n_studies']}", "studies"),
        (f"{cat['n_cells']/1000:.0f}k", "cells"),
        ("4", "layers"),
        (">200k", "isoforms"),
    ]
    for i, (num, lab) in enumerate(kpis):
        x0 = 0.02 + i * 0.195
        card(ax, x0, 0.06, 0.18, 0.36, fc="#F4F8F7", ec=TEAL, lw=0.8, rs=0.03)
        ax.text(x0 + 0.09, 0.30, num, transform=ax.transAxes, ha="center", va="center",
                fontsize=11, fontweight="bold", color=TEAL)
        ax.text(x0 + 0.09, 0.14, lab, transform=ax.transAxes, ha="center", va="center",
                fontsize=5.8, color=SLATE)

    # ================================================================= b species + platform
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Catalog composition: species and platform")
    sp_ord = [s for s in ["human", "mouse", "other"] if s in cells_by_sp.index]
    x = np.arange(len(sp_ord))
    w = 0.36
    ax2 = ax.twinx()
    b1 = ax.bar(x - w / 2, [runs_by_sp.get(s, 0) for s in sp_ord], width=w,
                color=TEAL, edgecolor="white", label="runs", zorder=2)
    b2 = ax2.bar(x + w / 2, [cells_by_sp.get(s, 0) / 1000 for s in sp_ord], width=w,
                 color=CORAL, edgecolor="white", label="cells (×10³)", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([s.capitalize() for s in sp_ord], fontsize=7)
    ax.set_ylabel("Runs", color=TEAL)
    ax2.set_ylabel("isoquant_cells (×10³)", color=CORAL)
    ax.tick_params(axis="y", labelcolor=TEAL)
    ax2.tick_params(axis="y", labelcolor=CORAL)
    plat_txt = "  ·  ".join(
        f"{p}: {int(runs_by_plat.get(p, 0))} runs / {int(cells_by_plat.get(p, 0)):,} cells"
        for p in ["ONT", "PacBio"] if p in runs_by_plat.index
    )
    ax.text(0.02, 0.98, plat_txt, transform=ax.transAxes, ha="left", va="top",
            fontsize=5.2, color=SLATE)
    ax.legend([b1, b2], ["runs", "cells (×10³)"], loc="upper right",
              fontsize=5.5, frameon=False)
    S.style_ax(ax)
    ax.spines["right"].set_visible(True)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    # ================================================================= c biological systems
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "isoquant_cells by biological system")
    rows = []
    for s in sys_order:
        sub = st[st.sys2 == s]
        rows.append((s, int(sub.n_cells.sum()), int(sub.gse.nunique()), int(sub.n_runs.sum())))
    rows = sorted(rows, key=lambda r: r[1])  # ascending for barh
    y = np.arange(len(rows))
    vals = np.array([r[1] for r in rows], dtype=float) / 1000.0
    ax.barh(y, vals, color=[colors_sys.get(r[0], GREY) for r in rows],
            edgecolor="white", height=0.72, zorder=2)
    for yi, (name, ncell, ngse, nrun) in enumerate(rows):
        ax.text(vals[yi] + max(vals) * 0.02, yi,
                f"{ngse} stud. · {nrun} runs · {ncell:,}",
                va="center", fontsize=5.3, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.5)
    ax.set_xlabel("isoquant_cells (×10³)")
    ax.set_xlim(0, max(vals) * 1.55)
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # ================================================================= d dual scope
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Two denominators (do not mix)")
    ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    card(ax, 0.04, 0.52, 0.92, 0.42, fc="#E8F2F1", ec=TEAL, lw=1.15, rs=0.04)
    ax.text(0.08, 0.86, "Catalog  (portal / GA / Fig. 1a–c)",
            transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=TEAL)
    ax.text(0.08, 0.62,
            f"{cat['n_runs']} IsoQuant-complete runs  ·  {cat['n_studies']} studies\n"
            f"{cat['n_cells']:,} cells  (isoquant_cells; write >800k)\n"
            f">200k observed transcript isoforms  ·  4 evidence layers",
            transform=ax.transAxes, fontsize=6.0, color=INK, va="top", linespacing=1.45)

    card(ax, 0.04, 0.06, 0.92, 0.40, fc="#F8EFEA", ec=CORAL, lw=1.15, rs=0.04)
    ax.text(0.08, 0.38, "Labeled analysis subset  (Fig. 2–3 statistics)",
            transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=CORAL)
    ax.text(0.08, 0.16,
            f"{ch['n_studies']} studies  ·  {ch['n_runs']} runs  ·  "
            f"{ch['n_annotated_cells']:,} annotated cells\n"
            f"Requires donor × cell-type labels for DIU / APA / ASE tables.",
            transform=ax.transAxes, fontsize=6.0, color=INK, va="top", linespacing=1.45)

    fig.suptitle(
        "scTHREAD · Content and construction of the single-cell long-read resource",
        x=0.08, y=0.97, ha="left", fontsize=10, fontweight="bold", color=INK,
    )
    fig.text(
        0.08, 0.015,
        f"Source: sample_registry.tsv  (isoquant_status=done) → "
        f"{cat['n_runs']} runs / {cat['n_studies']} studies / {cat['n_cells']:,} isoquant_cells.  "
        f"Does not replace NAR_Fig1 (dense atlas version).",
        fontsize=5.3, color=SLATE,
    )

    # save under distinct name — never overwrite NAR_Fig1
    out = S.save(fig, "NAR_Fig1_narstyle")
    print(f"wrote {out}  catalog={cat['n_runs']} runs / {cat['n_cells']:,} cells")
    return out


if __name__ == "__main__":
    render()
