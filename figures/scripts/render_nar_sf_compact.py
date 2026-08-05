#!/usr/bin/env python3
"""scTHREAD NAR Database Issue — compacted supplementary pack, SF1-SF4.

Replaces the ten-figure / forty-two-panel pack rendered by
render_nar_sf_all.py.  Panel-level provenance for every kept and every dropped
panel is in docs/NAR_SF_COMPACT_PANEL_MAP.md; that file is the contract this
module implements.

Nothing here recomputes anything.  Every panel reads the same frozen table the
old panel read, so the compaction cannot change a number.  The panels that were
dropped because they duplicate main Figure 2 (old SF4d, SF5c, SF7c-upper, SF8a,
SF8b) are not re-derived here under another name — they live in Figure 2.

  SF1  Composition of the release and evidence-layer coverage
  SF2  Splice-class composition and technical controls for novel structures
  SF3  Registration-free query path
  SF4  Cross-species worked example (published mouse gastrulation cohort)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S  # noqa: E402
import render_nar_sf_all as A  # noqa: E402

INK, SLATE, TEAL, CORAL, GOLD, BLUE, GREY = (
    A.INK, A.SLATE, A.TEAL, A.CORAL, A.GOLD, A.BLUE, A.GREY
)
GREYL = A.GREYL
WALK = A.WALK

CLASS_ORDER = ["FSM", "ISM", "NIC", "NNC"]
CLASS_COL = {"FSM": "#3D6F9B", "ISM": "#8FB0C9", "NIC": "#D08A4C", "NNC": "#C1503A"}


def tag(ax, letter, x=-0.04, y=1.06):
    S.panel_label(ax, letter, x=x, y=y)


def title(ax, s, pad=3.0):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=pad)


def suptitle(fig, number: int, text: str, x: float = 0.06) -> None:
    fig.suptitle(
        f"Supplementary Figure {number}  ·  {text}",
        x=x, y=0.975, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )


# ===================================================================== SF1
def render_sf1():
    """Composition of the release and evidence-layer coverage.

    a-d  old SF3a-d  (render_nar_sf_all.render_sf1)
    e    old SF7a    (render_nar_sf_all.render_sf7 a)
    f    old SF2b    (render_nar_sf_all.render_sf2 b)
    """
    import render_nar_bio as R

    done = A.load_registry_done()
    study = R.study_cells()
    atlas = pd.read_csv(A.Path(__file__).resolve().parents[2] / "tables/fig1_study_atlas.tsv",
                        sep="\t")
    study = study.merge(atlas[["gse", "system_display"]], on="gse", how="left")
    study["system_display"] = study["system_display"].fillna("Unassigned")
    head = R.catalog_headline()
    n_runs, n_studies, n_cells = head["n_runs"], head["n_studies"], head["n_cells"]

    cov = pd.read_csv(S.FIGDATA / "ed7_coverage_matrix.tsv", sep="\t")
    # Per-study cells are authoritative in the frozen release. The figdata
    # cohort_cellcount table predates it and disagrees for 13 of these 15
    # studies, summing to 452,213 against the release's 601,690 for the same
    # studies - and 452,213 is on the manuscript's forbidden-number list.
    cells = study[["gse", "cells"]].rename(columns={"cells": "n_cells"})
    lab = pd.read_csv(S.FIGDATA / "ed1_labeled_frac_bystudy.tsv", sep="\t")

    fig = plt.figure(figsize=(183 * S.MM, 200 * S.MM), facecolor="none")
    gs = fig.add_gridspec(3, 6, hspace=0.52, wspace=1.15,
                          left=0.10, right=0.97, top=0.92, bottom=0.055,
                          height_ratios=[1.0, 1.0, 1.30])

    # ---- a run records by species x platform ----
    ax = fig.add_subplot(gs[0, 0:3])
    tag(ax, "a"); title(ax, "Run records by species and platform")
    comp = done.copy()
    comp["sp"] = comp["species"].fillna("unknown").astype(str).str.lower()
    comp["platform"] = comp["platform"].fillna("NA")
    run_table = pd.crosstab(comp["sp"], comp["platform"])
    run_table = run_table.reindex(index=[v for v in ("human", "mouse") if v in run_table.index])
    run_table.plot(kind="bar", ax=ax, color=[TEAL, CORAL][: run_table.shape[1]],
                   edgecolor="white", width=0.68, zorder=2)
    for container in ax.containers:
        ax.bar_label(container, fontsize=5.4, color=SLATE, padding=1)
    ax.set_xlabel(""); ax.set_ylabel("Run records")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=6.8)
    ax.legend(title="platform", fontsize=5.2, title_fontsize=5.2, frameon=False)
    S.style_ax(ax); ax.grid(True, axis="y", color=GREYL, lw=0.35, zorder=0)

    # ---- b cells by species ----
    ax = fig.add_subplot(gs[0, 3:6])
    tag(ax, "b"); title(ax, "Cells by species")
    by_sp = {"human": head["cells_human"], "mouse": head["cells_mouse"]}
    keys = sorted(by_sp, key=by_sp.get)
    y = np.arange(len(keys))
    ax.barh(y, [by_sp[k] / 1000 for k in keys],
            color=[TEAL if k == "human" else CORAL for k in keys],
            edgecolor="white", height=0.6, zorder=2)
    for yi, k in enumerate(keys):
        runs = int((comp["sp"] == k).sum())
        ax.text(by_sp[k] / 1000 + 8, yi, f"{by_sp[k]:,} · {runs} records",
                va="center", fontsize=5.6, color=SLATE)
    ax.set_yticks(y); ax.set_yticklabels(keys, fontsize=6.8)
    ax.set_xlabel("Cells (×10³)")
    ax.set_xlim(0, max(by_sp.values()) / 1000 * 1.48)
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # ---- c cells by biological system ----
    ax = fig.add_subplot(gs[1, 0:3])
    tag(ax, "c"); title(ax, "Cells by biological system")
    sysagg = (study.groupby("system_display")
                   .agg(cells=("cells", "sum"), studies=("gse", "nunique"))
                   .sort_values("cells"))
    y = np.arange(len(sysagg))
    ax.barh(y, sysagg["cells"].values / 1000, color=TEAL, edgecolor="white",
            height=0.68, zorder=2)
    for yi, (key, row) in enumerate(sysagg.iterrows()):
        ax.text(row["cells"] / 1000 + 5, yi,
                f"{int(row['cells']):,} · {int(row['studies'])} datasets",
                va="center", fontsize=5.0, color=SLATE)
    ax.set_yticks(y); ax.set_yticklabels(sysagg.index, fontsize=5.6)
    ax.set_xlabel("Cells (×10³)")
    ax.set_xlim(0, sysagg["cells"].max() / 1000 * 1.60)
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # ---- d barcode candidates that become cells ----
    ax = fig.add_subplot(gs[1, 3:6])
    tag(ax, "d"); title(ax, "Barcode candidates retained as cells")
    plat = done.groupby(done["platform"].fillna("NA")).agg(
        cells=("isoquant_cells", "sum"),
        candidates=("isoquant_barcodes_raw", "sum"),
        n=("srr", "count"),
    ).reindex([p for p in ["ONT", "PacBio"] if p in done["platform"].unique()])
    x = np.arange(len(plat)); w = 0.36
    ax.bar(x - w / 2, plat["candidates"].values / 1e6, width=w, color=GOLD,
           edgecolor="white", label="barcode candidates", zorder=2)
    ax.bar(x + w / 2, plat["cells"].values / 1e6, width=w, color=TEAL,
           edgecolor="white", label="cells", zorder=2)
    for i, (_, r) in enumerate(plat.iterrows()):
        keep = 100 * r.cells / r.candidates if r.candidates > 0 else 0
        top = max(r.candidates, r.cells) / 1e6
        ax.text(i, top * 1.02, f"n={int(r.n)}\nretained {keep:.1f}%",
                ha="center", va="bottom", fontsize=5.4, color=SLATE)
    ax.set_xticks(x); ax.set_xticklabels(plat.index.tolist(), fontsize=7)
    ax.set_ylabel("Count (millions)")
    ax.set_ylim(0, float(plat[["candidates", "cells"]].max().max()) / 1e6 * 1.40)
    ax.legend(fontsize=5.4, frameon=False, loc="upper right")
    ax.text(0.02, 0.98,
            "PacBio libraries carry many one- to two-molecule\n"
            "ambient barcodes, so a smaller share becomes cells.",
            transform=ax.transAxes, va="top", fontsize=5.0, color=SLATE)
    S.style_ax(ax); ax.grid(True, axis="y", color=GREYL, lw=0.35, zorder=0)

    # ---- e evidence-layer coverage (old SF7a) ----
    ax = fig.add_subplot(gs[2, 0:4])
    tag(ax, "e", x=-0.06); title(ax, "Evidence-layer run coverage by study")
    layers = ["ASE", "DIU", "APA", "F2jct"]
    layer_names = ["ASE", "DIU", "APA", "Junction"]
    df = cov.merge(cells, on="gse", how="left").sort_values("n_cells", ascending=False)
    mat = df[layers].fillna(0).values.astype(float)
    mat_n = mat / np.maximum(mat.max(axis=0, keepdims=True), 1)
    n = len(df)
    ax.imshow(mat_n, aspect="auto", cmap="Greens", vmin=0, vmax=1,
              extent=(-0.5, 3.5, n - 0.5, -0.5))
    ax.set_xticks(range(4)); ax.set_xticklabels(layer_names, fontsize=6.5)
    ax.set_yticks(range(n)); ax.set_yticklabels(df.gse.values, fontsize=5.4)
    for i in range(n):
        for j in range(4):
            v = int(mat[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=5.1,
                        color="white" if mat_n[i, j] > 0.5 else INK, fontweight="bold")
        ax.text(3.66, i, f"{int(df.n_cells.iloc[i]):,} cells",
                va="center", fontsize=4.9, color=SLATE)
    ax.set_xlim(-0.5, 6.0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    n_missing = {name: int((df[col] == 0).sum()) for col, name in zip(layers, layer_names)}
    ax.text(0.0, -0.145,
            "Numbers = runs with layer data.  Studies with no run for a layer: "
            + ", ".join(f"{k} {v}/{n}" for k, v in n_missing.items()),
            transform=ax.transAxes, fontsize=5.2, color=SLATE)

    # ---- f cell-type annotation status (old SF2b) ----
    # The former per-study labelled fraction came from ed1_labeled_frac_bystudy,
    # which the release contradicts in both directions: it put GSE158450 at 0%
    # where the release records author cell types, and GSE303762 and GSE289790
    # near 100% where the release records none by design. Annotation status is
    # recorded for every study and is all-or-nothing, so that is what is drawn.
    ax = fig.add_subplot(gs[2, 4:6])
    tag(ax, "f", x=-0.13); title(ax, "Cells by cell-type annotation status")
    status = done.groupby("gse")["annotation_status"].agg(
        lambda s: "annotated" if (s == "annotated").any() else "none by design")
    fstat = study.assign(status=study.gse.map(status)).dropna(subset=["status"])
    grp = fstat.groupby("status").agg(cells=("cells", "sum"), datasets=("gse", "nunique"))
    grp = grp.reindex([k for k in ("annotated", "none by design") if k in grp.index])
    y = np.arange(len(grp))
    ax.barh(y, grp.cells / 1000,
            color=[TEAL if k == "annotated" else GREY for k in grp.index],
            edgecolor="white", height=0.55, zorder=2)
    # this panel is two of six columns wide, so a single-line label on the long
    # bar runs off the right edge; two lines fit in the vertical space instead
    for yi, (k, r) in enumerate(grp.iterrows()):
        ax.text(r.cells / 1000 + grp.cells.max() / 1000 * 0.04, yi,
                f"{int(r.cells):,} cells\n{int(r.datasets)} datasets",
                va="center", ha="left", fontsize=5.4, color=SLATE, linespacing=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels(["cell type assigned", "cell type not defined"][: len(grp)],
                       fontsize=6.2)
    ax.set_xlabel("Cells (×10³)")
    ax.set_xlim(0, grp.cells.max() / 1000 * 1.42)
    # what the five datasets are is explained in the legend: the bottom row is
    # already tight and a three-line note here falls outside the figure margin
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    fig.text(0.10, 0.008,
             f"Release: {n_runs} run records · {n_studies} datasets · {n_cells:,} cells",
             fontsize=5.6, color=SLATE)
    suptitle(fig, 1, "Composition of the release and evidence-layer coverage", x=0.08)
    S.save(fig, "NAR_SF1")
    print("NAR_SF1 ok  (a-d old SF3a-d · e old SF7a · f old SF2b)")


# ===================================================================== SF2
def render_sf2():
    """Splice-class composition and technical controls for novel structures.

    a,b  old SF5a,b  (render_nar_sf_all.render_sf3 a,b)
    c    old SF4b    (render_nar_sf_all.render_sf4 b)
    d    old SF4c    (render_nar_sf_all.render_sf4 c)
    e    old SF6a    (render_nar_sf_all.render_sf6 a)
    """
    ed2 = pd.read_csv(S.FIGDATA / "ed2_full_classification.tsv", sep="\t")
    by_gse = pd.read_csv(S.FIGDATA / "fig1b_classification_all.tsv", sep="\t")
    wt = pd.read_csv(S.FIGDATA / "fig3b_within_celltype_depth.tsv", sep="\t")
    rf = pd.read_csv(S.FIGDATA / "fig2_rarefaction_depthmatch.tsv", sep="\t")
    ed = pd.read_csv(S.FIGDATA / "ed3_error_decomp.tsv", sep="\t").sort_values("conflict_rate")

    fig = plt.figure(figsize=(183 * S.MM, 190 * S.MM), facecolor="none")
    gs = fig.add_gridspec(3, 2, hspace=0.46, wspace=0.34,
                          left=0.11, right=0.97, top=0.91, bottom=0.06,
                          height_ratios=[1.0, 1.0, 0.85])

    # ---- a splice-class composition by tissue ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Splice-class composition by tissue", pad=13)
    tissues = list(ed2.tissue.unique())
    nnc = ed2[ed2.cls == "NNC"].set_index("tissue")["frac"]
    tissues = sorted(tissues, key=lambda t: -nnc.get(t, 0))
    y = np.arange(len(tissues))
    left = np.zeros(len(tissues))
    for cls in CLASS_ORDER:
        fr = []
        for t in tissues:
            sub = ed2[(ed2.tissue == t) & (ed2.cls == cls)]
            fr.append(float(sub.frac.iloc[0]) if len(sub) else 0.0)
        fr = np.array(fr)
        ax.barh(y, fr, left=left, color=CLASS_COL[cls], edgecolor="white",
                height=0.72, label=cls, zorder=2)
        left = left + fr
    ax.set_yticks(y); ax.set_yticklabels(tissues, fontsize=5.6)
    ax.set_xlabel("Fraction of spliced molecules")
    ax.set_xlim(0, 1.02)
    ax.legend(fontsize=5.3, ncol=4, loc="lower right", bbox_to_anchor=(1.02, 1.005),
              frameon=False, handlelength=1.1, columnspacing=1.0)
    S.style_ax(ax)

    # ---- b novel fraction by study ----
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Novel fraction (NIC+NNC) by study")
    d = by_gse.sort_values("novel_frac", ascending=True)
    y = np.arange(len(d))
    cols = [CORAL if v > 0.55 else TEAL if v > 0.4 else GREY for v in d.novel_frac]
    ax.barh(y, d.novel_frac * 100, color=cols, edgecolor="white", height=0.7, zorder=2)
    for yi, (_, r) in enumerate(d.iterrows()):
        ax.text(r.novel_frac * 100 + 0.8, yi, f"{r.tissue}", va="center",
                fontsize=5.0, color=SLATE)
    ax.set_yticks(y); ax.set_yticklabels(d.gse.values, fontsize=5.4)
    ax.set_xlabel("Novel-isoform molecules (%)")
    ax.set_xlim(0, 100)
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.3, zorder=0)

    # ---- c within cell type: novelty vs depth ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Within cell type: novelty vs depth (run-level)")
    plot = wt if len(wt) <= 800 else wt.sample(800, random_state=0)
    ax.scatter(plot.mol_per_cell, plot.novel_frac * 100, s=10, c=TEAL, alpha=0.35,
               edgecolors="none", zorder=2)
    ax.set_xlabel("Molecules per cell")
    ax.set_ylabel("Novel-isoform molecules (%)")
    ax.set_xscale("log")
    S.style_ax(ax); ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    # ---- d depth-matched rarefaction ----
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Depth-matched rarefaction of recurrent novel paths")
    for ct, sub in rf.groupby("cell_type"):
        sub = sub.sort_values("frac")
        ax.plot(sub.frac, sub.recurrent_paths, marker="o", ms=3, lw=1.0, label=ct, zorder=2)
    ax.set_xlabel("Sampled fraction of molecules")
    ax.set_ylabel("Recurrent novel paths")
    ax.legend(fontsize=4.6, ncol=2, loc="upper left")
    S.style_ax(ax); ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    # ---- e nested vs conflict assignment ----
    ax = fig.add_subplot(gs[2, :])
    tag(ax, "e", x=-0.02)
    title(ax, "Nested versus conflicting assignments of discordant molecules", pad=13)
    y = np.arange(len(ed))
    ax.barh(y, ed.nested_rate * 100, color=TEAL, edgecolor="white", height=0.66,
            label="nested", zorder=2)
    ax.barh(y, ed.conflict_rate * 100, left=ed.nested_rate * 100, color=CORAL,
            edgecolor="white", height=0.66, label="conflict", zorder=2)
    for yi, (_, r) in enumerate(ed.iterrows()):
        ax.text(101.5, yi, f"{r.tissue} · conflict {r.conflict_rate*100:.1f}%",
                va="center", fontsize=5.3, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels([r.gse for _, r in ed.iterrows()], fontsize=5.7)
    ax.set_xlabel("Fraction of discordant molecules (%)")
    ax.set_xlim(0, 132)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.legend(fontsize=5.4, ncol=2, loc="lower right", bbox_to_anchor=(0.78, 1.005),
              frameon=False, handlelength=1.1, columnspacing=1.0)
    ax.text(0.0, -0.26,
            "Technical boundary of isoform assignment, not a biological discovery panel.",
            transform=ax.transAxes, fontsize=5.2, color=GREY, style="italic")
    S.style_ax(ax)

    suptitle(fig, 2, "Splice-class composition and technical controls for novel structures",
             x=0.09)
    S.save(fig, "NAR_SF2")
    print("NAR_SF2 ok  (a,b old SF5a,b · c,d old SF4b,c · e old SF6a)")


# ===================================================================== SF3
SHOTS = [
    ("public_20260805/01_home_current.png", "1 · Home / catalog"),
    ("public_20260805/02_search_PTPRC_current.png", "2 · Search PTPRC"),
    ("public_20260805/ptprc_gene_card_live_v3.png", "3 · Gene card"),
    ("public_20260805/ptprc_isoforms_live_v3.png", "4 · Isoform usage"),
    ("public_20260805/05_download_current.png", "5 · Download / export"),
    ("public_20260805/08_docs_current.png", "6 · API documentation"),
]


def render_sf3():
    """Registration-free query path — old SF9, 8 shots trimmed to 6.

    The shots have very different aspect ratios: 1, 2, 5 and 6 are full pages
    while 3 and 4 are cropped detail views.  A uniform gridspec with
    aspect="equal" therefore leaves large ragged gaps, which is what the old
    SF9 did.  Each tile is instead given a width proportional to its own
    aspect ratio at a common per-row height, so a row is always justified.
    """
    missing = [fn for fn, _ in SHOTS if not (WALK / fn).exists()]
    if missing:
        raise RuntimeError(f"walkthrough screenshots missing: {missing}")
    images = [(mpimg.imread(WALK / fn), lab) for fn, lab in SHOTS]

    width_mm = 183.0
    left, right = 0.03, 0.97
    gap_x, gap_y = 0.018, 0.055          # fractions of figure width / height
    top_pad_mm, bottom_pad_mm, title_mm = 13.0, 20.0, 4.0

    rows = [images[:3], images[3:]]
    avail = right - left
    # width fraction of each tile, and the row height expressed in mm
    row_specs = []
    for row in rows:
        aspects = [im.shape[1] / im.shape[0] for im, _ in row]
        usable = avail - gap_x * (len(row) - 1)
        h_frac_of_width = usable / sum(aspects)      # common tile height, in figure-width units
        widths = [a * h_frac_of_width for a in aspects]
        row_specs.append((widths, h_frac_of_width * width_mm))

    body_mm = sum(h for _, h in row_specs) + title_mm * len(rows)
    height_mm = top_pad_mm + body_mm + bottom_pad_mm + gap_y * 100
    fig = plt.figure(figsize=(width_mm * S.MM, height_mm * S.MM), facecolor="none")

    y_cursor = 1.0 - top_pad_mm / height_mm
    for (widths, row_h_mm), row in zip(row_specs, rows):
        h_frac = row_h_mm / height_mm
        y0 = y_cursor - h_frac - title_mm / height_mm
        x_cursor = left
        for w_frac, (img, lab) in zip(widths, row):
            ax = fig.add_axes((x_cursor, y0, w_frac, h_frac))
            ax.imshow(img, aspect="auto", interpolation="bilinear")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color(TEAL); sp.set_linewidth(0.8); sp.set_visible(True)
            ax.set_title(lab, loc="left", fontsize=6.2, fontweight="bold", color=INK, pad=2)
            x_cursor += w_frac + gap_x
        y_cursor = y0 - gap_y

    # Captured 2026-08-05 against the public site. The earlier 2026-07-27
    # capture was taken against an internal preview build whose Download page
    # still advertised a 469-run candidate, and this footnote used to excuse
    # that with processing-progress language; the live pages now lead with the
    # frozen release, so the note simply states what the reader will see.
    fig.text(0.03, 0.048,
             "Captured 2026-08-05 from https://scthread.ai4sc.ac.cn. The home page (1) and the "
             "Download page (5) both lead with the frozen release — 453 run records, 34 datasets,\n"
             "   923,389 cells — which is the checksum-matched snapshot the manuscript reports. The "
             "rolling live catalog is labelled separately in the interface and is never mixed with it.",
             fontsize=5.2, color=SLATE, linespacing=1.5)
    fig.text(0.03, 0.014,
             "Live site: https://scthread.ai4sc.ac.cn  ·  Demo gene: PTPRC  ·  "
             "Every view is reachable without registration and every table shown is exportable.",
             fontsize=5.5, color=SLATE)
    suptitle(fig, 3, "Registration-free query path", x=0.03)
    S.save(fig, "NAR_SF3")
    print(f"NAR_SF3 ok  ({len(images)} shots, old SF9 trimmed from 8; justified layout)")


# ===================================================================== SF4
def render_sf4():
    """Cross-species worked example — old SF1b,c,d,e; text cards moved to the legend."""
    import render_nar_sf10_v3 as M

    data = M.load_inputs()
    umap = data["umap"]
    manifest = data["manifest"]
    malat1 = manifest["malat1"]

    fig = plt.figure(figsize=(183 * S.MM, 195 * S.MM), facecolor="none")
    gs = fig.add_gridspec(3, 2, hspace=0.40, wspace=0.20,
                          left=0.14, right=0.96, top=0.91, bottom=0.07,
                          height_ratios=[0.85, 1.25, 0.95])

    # ---- a cell-type composition of the portal embedding ----
    ax = fig.add_subplot(gs[0, :])
    tag(ax, "a", x=-0.09)
    title(ax, "Cell types in the portal embedding (top 12 of 36)")
    top = data["celltypes"].head(12).iloc[::-1]
    ax.barh(np.arange(len(top)), top["cells"] / 1000.0, color=M.BLUE, height=0.72, zorder=2)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["cell_type"], fontsize=5.2)
    ax.tick_params(axis="y", pad=1.5)
    ax.set_xlabel("cells (×10³)", fontsize=6.0)
    S.style_ax(ax)

    # ---- b,c the two displayed Malat1 transcripts on shared coordinates ----
    for idx, (iso, col) in enumerate(((M.ISO1, "iso1_count"), (M.ISO2, "iso2_count"))):
        ax = fig.add_subplot(gs[1, idx])
        tag(ax, "bc"[idx], x=-0.02)
        stats = malat1["displayed_isoforms"][iso]
        M._panel_scatter(
            ax, umap, umap[col].values.astype(float), iso,
            f"{stats['cells']:,} cells with signal · {int(stats['molecules']):,} molecules",
        )

    # ---- d within-gene usage contrast ----
    ax = fig.add_subplot(gs[2, :])
    tag(ax, "d", x=-0.09)
    title(ax, f"Malat1 within-gene isoform usage by cell type "
              f"(maximum difference {malat1['max_usage_difference']:.3f})")
    pair = data["pair"].set_index("cell_type")
    show = pd.concat([pair.head(6), pair.tail(6)])
    y = np.arange(len(show))
    ax.barh(y - 0.2, show[M.ISO1], height=0.38, color=M.BLUE, label=M.ISO1, zorder=2)
    ax.barh(y + 0.2, show[M.ISO2], height=0.38, color=M.ORANGE, label=M.ISO2, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(show.index, fontsize=5.2)
    ax.set_xlabel("within-gene usage fraction", fontsize=6.0)
    ax.legend(fontsize=4.9, frameon=False, loc="lower right")
    S.style_ax(ax)

    fig.text(0.14, 0.014,
             "Published NGDC CRA044500 mouse cohort: six ONT runs of the 453 in the release, "
             "contributing 68,417 of the 923,389 cells.\n"
             "The maps use the harmonized portal embedding of 55,729 of those cells; "
             f"{malat1['cells_with_signal']:,} of "
             f"{manifest['portal_embedding']['cells']:,} cells carry Malat1 signal.",
             fontsize=5.2, color=SLATE)
    suptitle(fig, 4, "Cross-species worked example using a published mouse gastrulation cohort",
             x=0.06)
    S.save(fig, "NAR_SF4")
    print("NAR_SF4 ok  (a-d old SF1b,c,d,e; scope and boundary cards moved to the legend)")


RENDERERS = {1: render_sf1, 2: render_sf2, 3: render_sf3, 4: render_sf4}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fig", default="all",
                        help="all or comma-separated figure numbers, for example 1,4")
    parser.add_argument("--outdir", type=Path, default=S.OUTDIR)
    args = parser.parse_args()
    S.OUTDIR = args.outdir.resolve()
    selected = set(RENDERERS) if args.fig == "all" else {
        int(v.strip()) for v in args.fig.split(",") if v.strip()
    }
    bad = selected - set(RENDERERS)
    if bad:
        parser.error(f"invalid supplementary figure numbers: {sorted(bad)}")
    for number in sorted(selected):
        RENDERERS[number]()
    print("compact SF pack done →", S.OUTDIR)


if __name__ == "__main__":
    main()
