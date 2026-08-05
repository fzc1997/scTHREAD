#!/usr/bin/env python3
"""scTHREAD NAR Database Issue — Supplementary Figures SF1–SF8 (minimum pack).

Nature-style via nar_style.py. All numbers from frozen figdata / registry.
No invented values. Not discovery claims — supporting / QC / dual-scope.
"""
from __future__ import annotations

import os

import argparse
import json
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

INK, SLATE, TEAL, CORAL, GOLD, BLUE, GREY, CREAM, SOFT = (
    S.INK, S.SLATE, S.TEAL, S.CORAL, S.GOLD, S.BLUE, S.GREY, S.CREAM, S.SOFT
)
GREYL = "#D9DCE0"
PURPLE = "#6B5B95"
REG = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
CURRENT_REG = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv"
)
WALK = Path(__file__).resolve().parents[2] / "figures/website_walkthrough"
WALK_CORRECTED = WALK / "corrected_20260727"
F2 = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD") + "/results/paper1/f2_grammar/figdata")
P0 = Path(__file__).resolve().parents[2] / "tables/p0_biological_unit_rerun"


def tag(ax, letter, x=-0.04, y=1.06):
    S.panel_label(ax, letter, x=x, y=y)


def title(ax, s, pad=3.0):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=pad)


def card(ax, x, y, w, h, fc=CREAM, ec=GREYL, lw=0.7, rs=0.04):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.01,rounding_size={rs}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=ax.transAxes, clip_on=False,
    ))


# Supplementary figures must be numbered in order of first mention in the text.
# Content is authored under a stable key; this map is the only place the printed
# number lives, so re-ordering the narrative is a one-line change here and can
# never leave a title and a filename disagreeing.
SF_NUMBER = {1: 3, 2: 2, 3: 5, 4: 4, 5: 8, 6: 6, 7: 7, 8: 9, 9: 10, 10: 1}


def sf_title(content_key: int, title: str) -> str:
    return f"Supplementary Figure {SF_NUMBER[content_key]}  \u00b7  {title}"


def sf_name(content_key: int) -> str:
    return f"NAR_SF{SF_NUMBER[content_key]}"


def load_registry_done():
    """Return the released manuscript snapshot.

    Reads the frozen release manifest through render_nar_bio so that the
    modality audit and the stale-row correction applied at release time are
    honoured; re-deriving from sample_registry.tsv reintroduced both defects.
    """
    import render_nar_bio as R

    done = R.registry_done().copy()
    for c in ["annotation_cells", "isoquant_cells", "isoquant_barcodes_raw"]:
        if c in done.columns:
            done[c] = pd.to_numeric(done[c], errors="coerce")
    # Cells are authoritative per study, not per run: several studies take an
    # author-supplied or STARsolo total with no per-run decomposition, so a sum
    # over these rows overstates the release. Check runs here, cells there.
    cells = int(R.study_cells()["cells"].sum())
    if len(done) != 453 or cells != 923_389:
        raise RuntimeError(
            f"Released snapshot changed: {len(done)} runs / {cells} cells"
        )
    if done["species"].fillna("").isin(["", "?"]).any():
        raise RuntimeError("Released snapshot retains unresolved species")
    return done


# ===================================================================== SF1
def render_sf1():
    """Catalog composition: what the release contains, by species, platform and system.

    Cells come from the per-study authority, not a sum over run rows: several
    studies take an author-supplied or STARsolo total that has no per-run
    decomposition, so summing rows overstates the release.
    """
    import render_nar_bio as R

    done = load_registry_done()
    study = R.study_cells()
    atlas = pd.read_csv(
        Path(__file__).resolve().parents[2] / "tables/fig1_study_atlas.tsv", sep="\t")
    study = study.merge(atlas[["gse", "system_display"]], on="gse", how="left")
    study["system_display"] = study["system_display"].fillna("Unassigned")

    head = R.catalog_headline()
    n_runs, n_studies, n_cells = head["n_runs"], head["n_studies"], head["n_cells"]

    fig = plt.figure(figsize=(183 * S.MM, 150 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.32,
                          left=0.10, right=0.98, top=0.90, bottom=0.09)

    # ---- a run records by species x platform ----
    ax = fig.add_subplot(gs[0, 0])
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
    ax = fig.add_subplot(gs[0, 1])
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
    ax.set_xlabel("Cells (\u00d710\u00b3)")
    ax.set_xlim(0, max(by_sp.values()) / 1000 * 1.42)
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # ---- c cells by biological system ----
    ax = fig.add_subplot(gs[1, 0])
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
                va="center", fontsize=5.2, color=SLATE)
    ax.set_yticks(y); ax.set_yticklabels(sysagg.index, fontsize=5.8)
    ax.set_xlabel("Cells (\u00d710\u00b3)")
    ax.set_xlim(0, sysagg["cells"].max() / 1000 * 1.52)
    S.style_ax(ax); ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # ---- d barcode candidates that become cells, by platform ----
    ax = fig.add_subplot(gs[1, 1])
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
    ax.set_ylim(0, float(plat[["candidates", "cells"]].max().max()) / 1e6 * 1.34)
    ax.legend(fontsize=5.4, frameon=False, loc="upper right")
    ax.text(0.02, 0.98,
            "PacBio libraries carry many one- to two-molecule\n"
            "ambient barcodes, so a smaller share becomes cells.",
            transform=ax.transAxes, va="top", fontsize=5.0, color=SLATE)
    S.style_ax(ax); ax.grid(True, axis="y", color=GREYL, lw=0.35, zorder=0)

    fig.text(0.10, 0.035,
             f"Release: {n_runs} run records \u00b7 {n_studies} datasets \u00b7 {n_cells:,} cells",
             fontsize=5.6, color=SLATE)
    fig.suptitle(
        sf_title(1, "Composition of the scTHREAD release"),
        x=0.08, y=0.965, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(1))
    print("NAR_SF1 ok")


# ===================================================================== SF2
def render_sf2():
    """Cohort QC."""
    qc = pd.read_csv(S.FIGDATA / "ed1_cohort_qc.tsv", sep="\t")
    lab = pd.read_csv(S.FIGDATA / "ed1_labeled_frac_bystudy.tsv", sep="\t")
    cells = pd.read_csv(S.FIGDATA / "cohort_cellcount.tsv", sep="\t")
    df = qc.merge(lab, on="gse", how="left").merge(cells, on="gse", how="left")
    df = df.sort_values("n_cells", ascending=True)

    fig = plt.figure(figsize=(183 * S.MM, 155 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.32, left=0.12, right=0.98, top=0.90, bottom=0.08)

    # a cells by study
    ax = fig.add_subplot(gs[0, :])
    tag(ax, "a", x=-0.02); title(ax, "Uniform cohort: annotated cells and label rate per study")
    y = np.arange(len(df))
    ax.barh(y, df.n_cells / 1000, color=TEAL, edgecolor="white", height=0.72, zorder=2)
    for yi, (_, r) in enumerate(df.iterrows()):
        lf = r.labeled_frac if pd.notna(r.get("labeled_frac")) else np.nan
        lab_s = f"  label {lf*100:.0f}%" if pd.notna(lf) else ""
        novel = f"  novel {r.novel_pct:.0f}%" if pd.notna(r.get("novel_pct")) else ""
        ax.text(r.n_cells / 1000 + 1, yi,
                f"{r.tissue} · {int(r.n_runs)} runs{lab_s}{novel}",
                va="center", fontsize=5.2, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels(df.gse.values, fontsize=5.6)
    ax.set_xlabel("Annotated cells (×10³)")
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.35, zorder=0)

    # b labeled frac
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "b"); title(ax, "Fraction of cells with cell-type labels")
    d2 = df.dropna(subset=["labeled_frac"]).sort_values("labeled_frac")
    y = np.arange(len(d2))
    cols = [TEAL if v >= 0.9 else GOLD if v >= 0.5 else CORAL for v in d2.labeled_frac]
    ax.barh(y, d2.labeled_frac * 100, color=cols, edgecolor="white", height=0.7, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(d2.gse.values, fontsize=5.5)
    ax.set_xlabel("Labeled cells (%)")
    ax.set_xlim(0, 105)
    ax.axvline(90, color=GREYL, ls=":", lw=0.7)
    S.style_ax(ax)

    # c admitted reads
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "c"); title(ax, "Admitted long-read depth vs cells")
    ax.scatter(df.admitted_reads_M, df.n_cells / 1000, s=36, c=TEAL,
               edgecolors="white", linewidths=0.5, zorder=3)
    for _, r in df.iterrows():
        if r.n_cells > 40000 or r.admitted_reads_M > 400:
            ax.annotate(r.gse, (r.admitted_reads_M, r.n_cells / 1000),
                        fontsize=4.8, color=SLATE, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Admitted reads (millions)")
    ax.set_ylabel("Annotated cells (×10³)")
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    fig.suptitle(
        sf_title(2, "Uniform-cohort QC and cell-type label coverage"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(2))
    print("NAR_SF2 ok")


# ===================================================================== SF4
def render_sf4():
    """Depth / rarefaction controls."""
    nv = pd.read_csv(S.FIGDATA / "fig3_novel_vs_depth.tsv", sep="\t")
    wt = pd.read_csv(S.FIGDATA / "fig3b_within_celltype_depth.tsv", sep="\t")
    rf = pd.read_csv(S.FIGDATA / "fig2_rarefaction_depthmatch.tsv", sep="\t")
    jct = pd.read_csv(S.FIGDATA / "fig4_cross_study_junction.tsv", sep="\t")

    fig = plt.figure(figsize=(183 * S.MM, 150 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30, left=0.09, right=0.98, top=0.90, bottom=0.08)

    # a novel vs mol_per_cell
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Novel-isoform fraction vs depth (molecules / cell)")
    for gse, sub in nv.groupby("gse"):
        ax.scatter(sub.mol_per_cell, sub.novel_frac * 100, s=28, alpha=0.85,
                   edgecolors="white", linewidths=0.4, label=sub.tissue.iloc[0], zorder=3)
    ax.set_xlabel("Molecules per cell")
    ax.set_ylabel("Novel-isoform molecules (%)")
    ax.legend(fontsize=5.0, loc="best")
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    # b within CT
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Within cell type: novelty vs depth (run-level)")
    # subsample if huge
    plot = wt
    if len(plot) > 800:
        plot = plot.sample(800, random_state=0)
    ax.scatter(plot.mol_per_cell, plot.novel_frac * 100, s=10, c=TEAL, alpha=0.35,
               edgecolors="none", zorder=2)
    ax.set_xlabel("Molecules per cell")
    ax.set_ylabel("Novel-isoform molecules (%)")
    ax.set_xscale("log")
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    # c rarefaction
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Depth-matched rarefaction of recurrent novel paths")
    for ct, sub in rf.groupby("cell_type"):
        sub = sub.sort_values("frac")
        ax.plot(sub.frac, sub.recurrent_paths, marker="o", ms=3, lw=1.0, label=ct, zorder=2)
    ax.set_xlabel("Sampled fraction of molecules")
    ax.set_ylabel("Recurrent novel paths")
    ax.legend(fontsize=4.8, ncol=2, loc="upper left")
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    # d junction canonicity reminder
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Cross-study junction canonicity (biology axis)")
    ax.plot(jct.n_studies, jct.gtag_rate * 100, "o-", color=TEAL, ms=6, lw=1.4,
            markeredgecolor="white", markeredgewidth=0.6, zorder=3)
    ax.fill_between(jct.n_studies, jct.gtag_rate * 100, alpha=0.12, color=TEAL)
    r0, r1 = jct.gtag_rate.iloc[0] * 100, jct.gtag_rate.iloc[-1] * 100
    ax.annotate(f"{r0:.1f}% → {r1:.1f}% GT–AG\n"
                f"(x = independent studies, not depth)",
                xy=(jct.n_studies.iloc[-1], r1),
                xytext=(0.45, 0.25), textcoords="axes fraction",
                fontsize=6.0, color=TEAL, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=TEAL, lw=0.8))
    ax.set_xlabel("Independent studies reporting junction")
    ax.set_ylabel("GT–AG canonical rate (%)")
    ax.set_ylim(40, 102)
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.3, zorder=0)

    fig.suptitle(
        sf_title(4, "Depth and rarefaction controls for novelty claims"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(4))
    print("NAR_SF4 ok")


# ===================================================================== SF5
def render_sf5(observed_suffix: str = "", stem: str | None = None):
    stem = stem or sf_name(5)
    """Biological-unit negative controls + corrected inventory."""
    summary = pd.read_csv(
        P0 / "validation_summary.tsv", sep="\t", keep_default_na=False
    )
    observed = {
        axis: pd.read_csv(
            P0 / f"{axis.lower()}_observed{observed_suffix}.tsv", sep="\t"
        )
        for axis in ("ASE", "DIU", "APA")
    }
    null = summary[summary["dataset"].eq("null")].copy()
    null["analysis"] = null["analysis"].str.upper()
    null_rates = {
        axis: null.loc[
            null["analysis"].eq(axis), "raw_p_lt_0_05_fraction"
        ].astype(float).tolist()
        for axis in ("ASE", "DIU", "APA")
    }
    if any(len(values) != 3 for values in null_rates.values()):
        raise RuntimeError("SF5 requires three complete null seeds per analysis")

    fig = plt.figure(figsize=(183 * S.MM, 140 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30, left=0.09, right=0.98, top=0.90, bottom=0.08)

    # a null rates
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Biological-unit label-permutation nulls (3 seeds)")
    order = ["ASE", "DIU", "APA"]
    x = np.arange(len(order))
    vals = np.array([np.mean(null_rates[axis]) * 100 for axis in order])
    mins = np.array([np.min(null_rates[axis]) * 100 for axis in order])
    maxs = np.array([np.max(null_rates[axis]) * 100 for axis in order])
    cols = [TEAL, CORAL, GOLD]
    ax.bar(
        x, vals, color=cols, edgecolor="white", width=0.65, zorder=2,
        yerr=np.vstack([vals - mins, maxs - vals]), capsize=2.5,
        error_kw={"elinewidth": 0.7, "ecolor": INK},
    )
    ax.axhline(5, color=INK, ls="--", lw=0.9, zorder=1, label="nominal 5%")
    for i, (v, lo, hi) in enumerate(zip(vals, mins, maxs)):
        ax.text(
            i, hi + 0.25, f"mean {v:.1f}%\n{lo:.1f}–{hi:.1f}%",
            ha="center", va="bottom", fontsize=5.2, color=SLATE,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=7)
    ax.set_ylabel("Fraction p < 0.05 under null (%)")
    ax.set_ylim(0, max(8, vals.max() + 2))
    ax.legend(fontsize=5.5, loc="upper right")
    S.style_ax(ax)

    # b tested vs sig
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Corrected maps: genes tested vs significant")
    axes_n = [
        (axis, len(observed[axis]), int(observed[axis].sig.sum()), color)
        for axis, color in (("ASE", TEAL), ("DIU", CORAL), ("APA", GOLD))
    ]
    x = np.arange(3)
    w = 0.36
    ax.bar(x - w / 2, [a[1] / 1000 for a in axes_n], width=w, color=GREYL, edgecolor="white",
           label="tested", zorder=2)
    ax.bar(x + w / 2, [a[2] / 1000 for a in axes_n], width=w,
           color=[a[3] for a in axes_n], edgecolor="white", label="significant (q)", zorder=2)
    for i, (name, n, s, _) in enumerate(axes_n):
        ax.text(i, max(n, s) / 1000 + 0.3, f"{s:,}/{n:,}", ha="center", fontsize=5.5, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([a[0] for a in axes_n])
    ax.set_ylabel("Genes (×10³)")
    ax.legend(fontsize=5.5)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    # c ASE null emphasis
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "ASE cell-type interactions: no screen discoveries")
    ax.axis("off")
    card(ax, 0.06, 0.15, 0.88, 0.72, fc="#E8F2F1", ec=TEAL, lw=1.1)
    ax.text(0.5, 0.72, "0", transform=ax.transAxes, ha="center", fontsize=28,
            fontweight="bold", color=TEAL)
    ase = observed["ASE"]
    ax.text(0.5, 0.52, f"significant cell-type ASE interactions\n"
            f"(of {len(ase):,} eligible genes; q and effect gate)",
            transform=ax.transAxes, ha="center", fontsize=6.5, color=INK)
    ax.text(0.5, 0.28,
            f"raw p<0.05 ≈ {(ase.pval < 0.05).mean()*100:.1f}%  ·  "
            f"null mean ≈ {np.mean(null_rates['ASE'])*100:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=5.8, color=SLATE)
    ax.text(0.5, 0.12, "No ASE discovery passed this screen; this is not evidence of biological absence.",
            transform=ax.transAxes, ha="center", fontsize=5.4, color=GREY, style="italic")

    # d honesty card
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Interpretation bounds (supporting, not discovery)")
    ax.axis("off")
    lines = [
        ("DIU significant", f"{int(observed['DIU'].sig.sum()):,} / {len(observed['DIU']):,} genes", CORAL),
        ("APA significant", f"{int(observed['APA'].sig.sum()):,} / {len(observed['APA']):,} genes", GOLD),
        ("ASE CT significant", f"0 / {len(ase):,} genes", TEAL),
        ("Negative controls", "3 full-workflow seeds per axis; 19 sources; 200 inner permutations", GREY),
    ]
    for i, (k, v, c) in enumerate(lines):
        y0 = 0.78 - i * 0.20
        card(ax, 0.05, y0 - 0.08, 0.90, 0.16, fc=CREAM, ec=GREYL, lw=0.6)
        ax.text(0.10, y0 + 0.02, k, transform=ax.transAxes, fontsize=6.2, fontweight="bold", color=c)
        ax.text(0.90, y0 + 0.02, v, transform=ax.transAxes, fontsize=6.2, ha="right", color=INK)

    fig.suptitle(
        sf_title(5, "Negative-control diagnostics for ASE / DIU / APA tables"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, stem)
    print(f"{stem} ok")


# ===================================================================== SF6
def render_sf6():
    """Error decomp."""
    ed = pd.read_csv(S.FIGDATA / "ed3_error_decomp.tsv", sep="\t")
    ed = ed.sort_values("conflict_rate", ascending=True)

    fig = plt.figure(figsize=(183 * S.MM, 120 * S.MM), facecolor="none")
    gs = fig.add_gridspec(1, 2, wspace=0.32, left=0.12, right=0.98, top=0.88, bottom=0.12)

    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Nested vs conflict assignments by study")
    y = np.arange(len(ed))
    ax.barh(y, ed.nested_rate * 100, color=TEAL, edgecolor="white", height=0.7,
            label="nested", zorder=2)
    ax.barh(y, ed.conflict_rate * 100, left=ed.nested_rate * 100, color=CORAL,
            edgecolor="white", height=0.7, label="conflict", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.gse}  {r.tissue}" for _, r in ed.iterrows()], fontsize=5.8)
    ax.set_xlabel("Fraction of discordant molecules (%)")
    ax.set_xlim(0, 100)
    ax.legend(fontsize=5.8, loc="lower right")
    S.style_ax(ax)

    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Conflict rate (technical boundary)")
    y = np.arange(len(ed))
    ax.barh(y, ed.conflict_rate * 100, color=CORAL, edgecolor="white", height=0.7, zorder=2)
    for yi, (_, r) in enumerate(ed.iterrows()):
        ax.text(r.conflict_rate * 100 + 0.8, yi, f"{r.conflict_rate*100:.1f}%",
                va="center", fontsize=5.5, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels(ed.gse.values, fontsize=5.8)
    ax.set_xlabel("Conflict rate (%)")
    ax.text(0.98, 0.05,
            "Method-supporting view of assignment ambiguity;\n"
            "not a biological discovery panel.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.4, color=GREY, style="italic")
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.3, zorder=0)

    fig.suptitle(
        sf_title(6, "Technical assignment decomposition (nested vs conflict)"),
        x=0.08, y=0.97, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(6))
    print("NAR_SF6 ok")


# ===================================================================== SF3
def render_sf3():
    """Full isoform classification landscape (FSM/ISM/NIC/NNC)."""
    ed2 = pd.read_csv(S.FIGDATA / "ed2_full_classification.tsv", sep="\t")
    by_gse = pd.read_csv(S.FIGDATA / "fig1b_classification_all.tsv", sep="\t")
    CLASS_ORDER = ["FSM", "ISM", "NIC", "NNC"]
    CLASS_COL = {"FSM": "#3D6F9B", "ISM": "#8FB0C9", "NIC": "#D08A4C", "NNC": "#C1503A"}

    fig = plt.figure(figsize=(183 * S.MM, 155 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30, left=0.10, right=0.98, top=0.90, bottom=0.08)

    # a tissue stacked
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Splice-class composition by tissue")
    tissues = list(ed2.tissue.unique())
    # order by NNC frac descending
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
    ax.set_yticks(y)
    ax.set_yticklabels(tissues, fontsize=5.6)
    ax.set_xlabel("Fraction of spliced molecules")
    ax.set_xlim(0, 1.02)
    ax.legend(fontsize=5.5, ncol=4, loc="lower right", bbox_to_anchor=(1.0, 1.01))
    S.style_ax(ax)

    # b study-level novel frac
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Novel fraction (NIC+NNC) by study")
    d = by_gse.sort_values("novel_frac", ascending=True)
    y = np.arange(len(d))
    cols = [CORAL if v > 0.55 else TEAL if v > 0.4 else GREY for v in d.novel_frac]
    ax.barh(y, d.novel_frac * 100, color=cols, edgecolor="white", height=0.7, zorder=2)
    for yi, (_, r) in enumerate(d.iterrows()):
        ax.text(r.novel_frac * 100 + 0.8, yi, f"{r.tissue}", va="center", fontsize=5.0, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels(d.gse.values, fontsize=5.4)
    ax.set_xlabel("Novel-isoform molecules (%)")
    ax.set_xlim(0, 100)
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.3, zorder=0)

    # c system heatmap of class
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Class profile by biological system (study means)")
    # map study to system via by_gse
    mat = by_gse.groupby("system")[["FSM", "ISM", "NIC", "NNC"]].mean()
    # order systems
    sys_ord = [s for s in ["Blood/marrow", "Brain", "Heart", "Cancer", "Differentiation", "Benchmark", "Mouse"]
               if s in mat.index]
    sys_ord += [s for s in mat.index if s not in sys_ord]
    mat = mat.reindex(sys_ord)[CLASS_ORDER]
    im = ax.imshow(mat.values * 100, aspect="auto", cmap="YlOrRd", vmin=0, vmax=60)
    ax.set_xticks(range(4))
    ax.set_xticklabels(CLASS_ORDER, fontsize=6.5)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=5.8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j] * 100
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5.5,
                    color="white" if v > 35 else INK, fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("% of spliced molecules", fontsize=5.2)
    cbar.ax.tick_params(labelsize=4.8)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(GREYL); sp.set_linewidth(0.5)

    # d global pie-like bars: overall molecule-weighted
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Four-class composition, renormalized over classified spliced reads")
    tot = ed2.groupby("cls")["n_reads"].sum().reindex(CLASS_ORDER)
    frac = tot / tot.sum()
    x = np.arange(len(CLASS_ORDER))
    ax.bar(x, frac.values * 100, color=[CLASS_COL[c] for c in CLASS_ORDER],
           edgecolor="white", width=0.7, zorder=2)
    for i, (c, f, n) in enumerate(zip(CLASS_ORDER, frac.values, tot.values)):
        ax.text(i, f * 100 + 1.2, f"{f*100:.1f}%\n{n/1e9:.2f}B",
                ha="center", va="bottom", fontsize=5.5, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels(CLASS_ORDER, fontsize=7)
    ax.set_ylabel("% of spliced molecules")
    ax.set_ylim(0, max(frac.values) * 100 * 1.35)
    ax.text(0.98, 0.05, "FSM full-splice-match · ISM incomplete\n"
            "NIC novel-in-catalog · NNC novel-not-in-catalog",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    fig.suptitle(
        sf_title(3, "Full isoform classification landscape (FSM / ISM / NIC / NNC)"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(3))
    print("NAR_SF3 ok")


# ===================================================================== SF7
def render_sf7():
    """Layer coverage + gene-loss honesty boundary."""
    cov = pd.read_csv(S.FIGDATA / "ed7_coverage_matrix.tsv", sep="\t")
    cells = pd.read_csv(S.FIGDATA / "cohort_cellcount.tsv", sep="\t")
    d14 = pd.read_csv(S.FIGDATA / "fig5_d14_summary.tsv", sep="\t").iloc[0]
    debt = pd.read_csv(S.FIGDATA / "fig1c_debt_saturation.tsv", sep="\t")
    # gene rescue: summarize only — avoid loading all for plot density
    rescue = pd.read_csv(S.FIGDATA / "fig5_gene_rescue.tsv", sep="\t",
                         usecols=["gene_id", "n_reads", "overlapped", "rescue_frac"])

    df = cov.merge(cells, on="gse", how="left").sort_values("n_cells", ascending=False)
    layers = ["ASE", "DIU", "APA", "F2jct"]
    layer_names = ["ASE", "DIU", "APA", "Junction"]

    fig = plt.figure(figsize=(183 * S.MM, 155 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30, left=0.10, right=0.98, top=0.90, bottom=0.08)

    # a coverage heatmap
    ax = fig.add_subplot(gs[0, :])
    tag(ax, "a", x=-0.02); title(ax, "Evidence-layer run coverage by study (uniform cohort)")
    mat = df[layers].fillna(0).values.astype(float)
    mat_n = mat / np.maximum(mat.max(axis=0, keepdims=True), 1)
    n = len(df)
    im = ax.imshow(mat_n, aspect="auto", cmap="Greens", vmin=0, vmax=1,
                   extent=(-0.5, 3.5, n - 0.5, -0.5))
    ax.set_xticks(range(4))
    ax.set_xticklabels(layer_names, fontsize=6.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(df.gse.values, fontsize=5.5)
    for i in range(n):
        for j in range(4):
            v = int(mat[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=5.2,
                        color="white" if mat_n[i, j] > 0.5 else INK, fontweight="bold")
        ax.text(3.7, i, f"{int(df.n_cells.iloc[i]):,} cells",
                va="center", fontsize=5.0, color=SLATE)
    ax.set_xlim(-0.5, 6.2)
    ax.set_xlabel("Numbers = runs with layer data")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)

    # b missing layers
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "b"); title(ax, "Studies lacking each layer")
    miss = [(name, int((df[col] == 0).sum())) for col, name in zip(layers, layer_names)]
    x = np.arange(len(miss))
    ax.bar(x, [m[1] for m in miss], color=[TEAL, CORAL, GOLD, BLUE],
           edgecolor="white", width=0.65, zorder=2)
    for i, (name, v) in enumerate(miss):
        ax.text(i, v + 0.15, f"{v}/{n}", ha="center", fontsize=6.0, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([m[0] for m in miss])
    ax.set_ylabel("Studies with zero runs")
    ax.set_ylim(0, n + 1)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    # c gene-loss honesty + debt
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "c"); title(ax, "Reference recovery & short-read gene-loss bound")
    ax.axis("off")
    # debt final
    debt_last = float(debt.cum_frac.iloc[-1])
    card(ax, 0.04, 0.55, 0.92, 0.40, fc="#E8F2F1", ec=TEAL, lw=0.9)
    ax.text(0.08, 0.88, "FSM reference recovery (main Fig. 2e support)",
            transform=ax.transAxes, fontsize=6.2, fontweight="bold", color=TEAL)
    ax.text(0.08, 0.68,
            f"Observed FSM catalogue fraction: {debt_last*100:.1f}%\n"
            f"Remaining annotation debt: {(1-debt_last)*100:.1f}%\n"
            f"after {int(debt.n_runs.iloc[-1])} cumulative runs",
            transform=ax.transAxes, fontsize=5.8, color=INK, va="top", linespacing=1.4)

    card(ax, 0.04, 0.05, 0.92, 0.45, fc="#F8EFEA", ec=CORAL, lw=0.9)
    n_tot = len(rescue)
    n_overlap = int(rescue.overlapped.sum()) if rescue.overlapped.dtype == bool else int(rescue.overlapped.astype(bool).sum())
    ax.text(0.08, 0.42, "Short-read gene-loss boundary (not a discovery claim)",
            transform=ax.transAxes, fontsize=6.2, fontweight="bold", color=CORAL)
    ax.text(0.08, 0.22,
            f"Genes losing ≥20% short-read coverage: {int(d14.n_lose20):,}\n"
            f"Genes losing ≥50%: {int(d14.n_lose50):,}  ·  overlapped: {int(d14.n_overlapped):,}\n"
            f"Global gene accuracy: {d14.global_gene_acc:.3f}\n"
            f"Immune enrichment: {d14.immune_verdict}",
            transform=ax.transAxes, fontsize=5.5, color=INK, va="top", linespacing=1.35)

    fig.suptitle(
        sf_title(7, "Evidence-layer coverage and reference / gene-loss bounds"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(7))
    print("NAR_SF7 ok")


# ===================================================================== SF10
def render_sf10_catalog_legacy():
    """Legacy species × platform split; its content is now incorporated into SF1."""
    done = load_registry_done()
    done["species"] = done["species"].fillna("unknown").astype(str).str.lower()
    # normalize species
    def sp_norm(s):
        s = s.lower()
        if "human" in s or s in ("homo sapiens", "hs", "hsa"):
            return "human"
        if "mouse" in s or "mus" in s or s in ("mm", "mmu"):
            return "mouse"
        if s in ("", "nan", "?", "unknown", "none"):
            return "unknown"
        return s
    done["sp"] = done["species"].map(sp_norm)
    done["platform"] = done["platform"].fillna("NA")

    fig = plt.figure(figsize=(183 * S.MM, 130 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.30, left=0.10, right=0.98, top=0.90, bottom=0.09)

    # a runs by species × platform
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Frozen-release runs by species × platform")
    ct = pd.crosstab(done["sp"], done["platform"])
    # keep main
    ct = ct.reindex(index=[i for i in ["human", "mouse", "unknown"] if i in ct.index])
    ct.plot(kind="bar", ax=ax, color=[TEAL, CORAL, GOLD, BLUE, GREY][: ct.shape[1]],
            edgecolor="white", width=0.75, zorder=2)
    ax.set_xlabel("")
    ax.set_ylabel("Number of runs")
    ax.legend(title="platform", fontsize=5.2, title_fontsize=5.2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=6.5)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    # b cells by species
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "isoquant_cells by species")
    g = done.groupby("sp")["isoquant_cells"].sum().sort_values(ascending=True)
    y = np.arange(len(g))
    cols = [TEAL if k == "human" else CORAL if k == "mouse" else GREY for k in g.index]
    ax.barh(y, g.values / 1000, color=cols, edgecolor="white", height=0.65, zorder=2)
    for yi, (k, v) in enumerate(g.items()):
        ax.text(v / 1000 + 5, yi, f"{int(v):,}  ({int((done.sp==k).sum())} runs)",
                va="center", fontsize=5.5, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels(g.index, fontsize=6.5)
    ax.set_xlabel("isoquant_cells (×10³)")
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.3, zorder=0)

    # c cells by platform
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "isoquant_cells by platform")
    g = done.groupby("platform")["isoquant_cells"].sum().sort_values(ascending=True)
    y = np.arange(len(g))
    ax.barh(y, g.values / 1000, color=BLUE, edgecolor="white", height=0.65, zorder=2)
    for yi, (k, v) in enumerate(g.items()):
        bc = done.loc[done.platform == k, "isoquant_barcodes_raw"].sum()
        ax.text(v / 1000 + 5, yi,
                f"{int(v):,} cells  ·  raw BC {int(bc):,}",
                va="center", fontsize=5.3, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels(g.index, fontsize=6.5)
    ax.set_xlabel("isoquant_cells (×10³)")
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.3, zorder=0)

    # d method × species
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Cell-calling method × species (run counts)")
    done["method"] = done["isoquant_cells_method"].fillna("NA")
    ct2 = pd.crosstab(done["method"], done["sp"])
    ct2 = ct2.reindex(columns=[c for c in ["human", "mouse", "unknown"] if c in ct2.columns])
    im = ax.imshow(ct2.values, aspect="auto", cmap="Blues")
    ax.set_xticks(range(ct2.shape[1]))
    ax.set_xticklabels(ct2.columns, fontsize=6.5)
    ax.set_yticks(range(ct2.shape[0]))
    ax.set_yticklabels(ct2.index, fontsize=5.2)
    for i in range(ct2.shape[0]):
        for j in range(ct2.shape[1]):
            v = int(ct2.values[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=5.5,
                        color="white" if v > ct2.values.max() * 0.5 else INK, fontweight="bold")
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(GREYL); sp.set_linewidth(0.5)
    ax.tick_params(length=0)

    fig.suptitle(
        sf_title(10, "Catalog composition by species and platform"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, sf_name(10))
    print(f"{sf_name(10)} ok")


# ===================================================================== SF9
def render_sf9(
    observed_suffix: str = "",
    stem: str | None = None,
    ptprc_snapshot: Path | None = None,
):
    """Second gene case (MS4A1) — not cherry-pick PTPRC only."""
    stem = stem or sf_name(9)
    tab = Path(__file__).resolve().parents[2] / "tables"
    ov_path = tab / "MS4A1_api_overview.json"
    use_path = tab / "MS4A1_isoform_usage_by_ct.tsv"
    ptprc_path = ptprc_snapshot or (
        tab / "PTPRC_api_overview_corrected_20260727.json"
    )
    if not ov_path.exists() or not use_path.exists():
        print(f"{sf_name(9)} skip: missing MS4A1 tables")
        return

    ov = json.loads(ov_path.read_text())
    usage = pd.read_csv(use_path, sep="\t")
    pt = json.loads(ptprc_path.read_text()) if ptprc_path.exists() else None
    gid = ov["gene"]["gid"]
    for axis, feature_name in (("diu", "n_iso"), ("apa", "n_pas")):
        corrected = pd.read_csv(
            P0 / f"{axis}_observed{observed_suffix}.tsv", sep="\t"
        )
        row = corrected.loc[corrected["gene"].eq(gid)]
        if len(row) != 1:
            raise RuntimeError(f"Expected one corrected {axis.upper()} row for {gid}")
        row = row.iloc[0]
        ov["analyses"][axis] = {
            "gene": gid,
            "pval": float(row["pval"]),
            "effect": float(row["effect_equal_donor"]),
            feature_name: int(row["n_features"]),
            "qval": float(row["qval"]),
            "sig": bool(row["sig"]),
        }
    ov["analyses"]["ase"] = None

    a = ov["analyses"]
    cov = ov["coverage"]
    g = ov["gene"]

    fig = plt.figure(figsize=(183 * S.MM, 150 * S.MM), facecolor="none")
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28, left=0.09, right=0.98, top=0.90, bottom=0.07)

    # a story cards
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Second worked gene: MS4A1 (CD20) multi-layer status")
    ax.axis("off")
    card(ax, 0.04, 0.55, 0.92, 0.40, fc="#E8F2F1", ec=TEAL, lw=0.9)
    ax.text(0.08, 0.88, f"MS4A1 / CD20  ·  {g['gid']}", transform=ax.transAxes,
            fontsize=7, fontweight="bold", color=TEAL)
    ax.text(0.08, 0.68,
            "B-cell surface marker; textbook lineage gene.\n"
            "Independent of PTPRC to show multi-layer query is general.",
            transform=ax.transAxes, fontsize=5.6, color=INK, va="top", linespacing=1.35)
    # three badges
    badges = [
        ("DIU", a["diu"]["sig"], a["diu"]["qval"], a["diu"]["effect"], f"n_iso={a['diu']['n_iso']}"),
        ("APA", a["apa"]["sig"], a["apa"]["qval"], a["apa"]["effect"], f"n_pas={a['apa']['n_pas']}"),
        ("ASE", None, None, None, "not eligible"),
    ]
    for i, (name, sig, q, eff, note) in enumerate(badges):
        x0 = 0.04 + i * 0.32
        col = TEAL if sig is True else GREY
        card(ax, x0, 0.08, 0.30, 0.40, fc=CREAM, ec=col, lw=1.0)
        ax.text(x0 + 0.15, 0.40, name, transform=ax.transAxes, ha="center",
                fontsize=7, fontweight="bold", color=col)
        status = "sig" if sig is True else ("n.s." if sig is False else "not tested")
        ax.text(x0 + 0.15, 0.28, status, transform=ax.transAxes,
                ha="center", fontsize=6.5, color=col)
        q_s = f"q={q:.4f}" if q is not None else "—"
        e_s = f"eff={eff:.2f}" if eff is not None else ""
        ax.text(x0 + 0.15, 0.14, f"{q_s}\n{e_s}\n{note}", transform=ax.transAxes,
                ha="center", fontsize=5.0, color=SLATE, va="top", linespacing=1.25)

    # b PTPRC vs MS4A1 comparison
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Same framework on two textbook genes")
    genes = []
    if pt is not None:
        pa = pt["analyses"]
        genes.append(("PTPRC", pa["diu"]["effect"], pa["apa"]["effect"],
                      1 if pa["diu"]["sig"] else 0, 1 if pa["apa"]["sig"] else 0,
                      0 if not pa.get("ase") or not pa["ase"]["sig"] else 1))
    genes.append(("MS4A1", a["diu"]["effect"], a["apa"]["effect"],
                  1 if a["diu"]["sig"] else 0, 1 if a["apa"]["sig"] else 0,
                  0))
    x = np.arange(len(genes))
    w = 0.28
    ax.bar(x - w, [g[1] for g in genes], width=w, color=CORAL, edgecolor="white", label="DIU effect", zorder=2)
    ax.bar(x, [g[2] for g in genes], width=w, color=GOLD, edgecolor="white", label="APA effect", zorder=2)
    ax.bar(x + w, [0.05 if g[5] else 0.0 for g in genes], width=w, color=TEAL, edgecolor="white",
           label="ASE CT sig (binary)", zorder=2)
    for i, g in enumerate(genes):
        ax.text(i - w, g[1] + 0.02, "sig" if g[3] else "ns", ha="center", fontsize=5.5,
                color=CORAL, fontweight="bold")
        ax.text(i, g[2] + 0.02, "sig" if g[4] else "ns", ha="center", fontsize=5.5,
                color=GOLD, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in genes], fontsize=7)
    ax.set_ylabel("Effect size (or ASE flag)")
    ax.legend(fontsize=5.2, loc="upper right")
    ax.set_ylim(0, max(0.7, max(g[1] for g in genes) + 0.15))
    ax.text(0.02, 0.98, "sig = passes the joint q and effect-size gate",
            transform=ax.transAxes, va="top", fontsize=5.2, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.3, zorder=0)

    # c isoform usage heatmap
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "MS4A1 isoform usage by cell type (API aggregates)")
    iso_cols = [c for c in usage.columns if c.startswith("ENST")]
    # focus main immune CTs with enough mol
    keep_ct = [c for c in ["B cell", "Monocyte", "Monocyte/Myeloid", "NK", "CD8 T",
                           "Dendritic cell", "Erythroid", "HSPC"] if c in set(usage.cell_type)]
    u = usage[usage.cell_type.isin(keep_ct)].set_index("cell_type").reindex(keep_ct)
    mat = u[iso_cols].values.astype(float)
    # short labels
    short = [c.replace("ENST00000", "…") for c in iso_cols]
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=0.7)
    ax.set_xticks(range(len(short)))
    ax.set_xticklabels(short, rotation=40, ha="right", fontsize=5.2)
    ax.set_yticks(range(len(keep_ct)))
    display_ct = [
        "Mono/Myeloid" if value == "Monocyte/Myeloid" else value
        for value in keep_ct
    ]
    ax.set_yticklabels(display_ct, fontsize=5.8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if v >= 0.08:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.8,
                        color="white" if v > 0.4 else INK)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("usage fraction", fontsize=5.2)
    cbar.ax.tick_params(labelsize=4.8)
    ax.text(0.0, -0.28, "B cells dominate MS4A1 molecules; usage still multi-isoform.",
            transform=ax.transAxes, fontsize=5.2, color=SLATE)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(GREYL); sp.set_linewidth(0.5)

    # d coverage + export path
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Portal / API coverage for MS4A1")
    ax.axis("off")
    rows = [
        ("Isoform features / molecules",
         f"{cov['isoforms']['features']} / {cov['isoforms']['molecules']:,}"),
        ("PAS features / molecules",
         f"{cov['pas']['features']} / {cov['pas']['molecules']:,}"),
        ("Studies (isoform / PAS)",
         f"{cov['isoforms']['studies']} / {cov['pas']['studies']}"),
        ("Cell types (isoform)", f"{cov['isoforms']['cell_types']}"),
        ("Browse URL", "scthread.ai4sc.ac.cn/browse?query=MS4A1"),
        ("API", "GET /api/gene/MS4A1/overview|isoforms"),
    ]
    for i, (k, v) in enumerate(rows):
        y0 = 0.88 - i * 0.14
        card(ax, 0.04, y0 - 0.08, 0.92, 0.12, fc=CREAM, ec=GREYL, lw=0.55)
        ax.text(0.08, y0 - 0.01, k, transform=ax.transAxes, fontsize=5.8, color=SLATE, va="center")
        ax.text(0.92, y0 - 0.01, v, transform=ax.transAxes, fontsize=5.6, color=INK,
                ha="right", va="center", fontweight="bold")

    fig.suptitle(
        sf_title(9, "Second multi-layer gene case (MS4A1 / CD20)"),
        x=0.08, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    S.save(fig, stem)
    print(f"{stem} ok")


# ===================================================================== SF8
def render_sf8():
    """Portal walkthrough gallery."""
    # prefer clearer shots when available
    shots = [
        ("corrected_20260727/01_home_current.png", "1 · Home / catalog"),
        ("corrected_20260727/02_search_PTPRC_current.png", "2 · Search PTPRC"),
        ("corrected_20260727/ptprc_gene_card_live_v3.png", "3 · Corrected gene card"),
        ("corrected_20260727/ptprc_isoforms_live_v3.png", "4 · Corrected isoforms"),
        ("corrected_20260727/05_download_current.png", "5 · Download"),
        ("corrected_20260727/06_analyze_current.png", "6 · Analyze"),
        ("corrected_20260727/07_about_current.png", "7 · About"),
        ("corrected_20260727/08_docs_current.png", "8 · Docs / API"),
    ]
    present = [(fn, lab) for fn, lab in shots if (WALK / fn).exists()]
    if not present:
        print("NAR_SF8 skip: no screenshots in", WALK)
        return

    n = len(present)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(183 * S.MM, 42 * S.MM * nrows + 20 * S.MM), facecolor="none")
    gs = fig.add_gridspec(nrows, ncols, hspace=0.28, wspace=0.12,
                          left=0.03, right=0.97, top=0.90, bottom=0.04)

    for i, (fn, lab) in enumerate(present):
        r, c = divmod(i, ncols)
        ax = fig.add_subplot(gs[r, c])
        img = mpimg.imread(WALK / fn)
        ax.imshow(img, aspect="equal", interpolation="bilinear")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(TEAL if i < 4 else GREY)
            sp.set_linewidth(0.9)
            sp.set_visible(True)
        ax.set_title(lab, loc="left", fontsize=6.0, fontweight="bold", color=INK, pad=2)

    fig.suptitle(
        sf_title(8, "Portal walkthrough (registration-free query path)"),
        x=0.03, y=0.98, ha="left", fontsize=9.5, fontweight="bold", color=INK,
    )
    fig.text(0.03, 0.01,
             "Live site: https://scthread.ai4sc.ac.cn  ·  Demo gene: PTPRC  ·  "
             "Purpose: show multi-layer evidence is queryable without local reanalysis.",
             fontsize=5.5, color=SLATE)
    S.save(fig, sf_name(8))
    print("NAR_SF8 ok")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fig",
        default="all",
        help="all or comma-separated supplementary figure numbers (for example 1,5,8,9)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=S.OUTDIR,
        help="Output directory; defaults to the submission figures directory.",
    )
    parser.add_argument(
        "--observed-suffix",
        default="_9999",
        help=(
            "Suffix for the observed ASE/DIU/APA tables used by SF5 and SF9. "
            "Defaults to _9999: the manuscript reports the 9,999-permutation "
            "screen (docs/NAR_OBSERVED_9999_SENSITIVITY_20260727.md), so the "
            "figures must be rebuilt from the same tables. Pass '' only to "
            "reproduce the superseded 1,000-permutation render."
        ),
    )
    parser.add_argument(
        "--sf5-stem",
        default=None,
        help="Override the SF5 output stem",
    )
    parser.add_argument(
        "--sf9-stem",
        default=None,
        help="Override the SF9 output stem",
    )
    parser.add_argument(
        "--ptprc-snapshot",
        type=Path,
        default=None,
        help="Optional PTPRC overview snapshot used by SF9",
    )
    args = parser.parse_args()
    S.OUTDIR = args.outdir.resolve()
    selected = set(range(1, 11)) if args.fig == "all" else {
        int(value.strip()) for value in args.fig.split(",") if value.strip()
    }
    invalid = selected - set(range(1, 11))
    if invalid:
        parser.error(f"invalid supplementary figure numbers: {sorted(invalid)}")
    renderers = {
        1: render_sf1, 2: render_sf2, 3: render_sf3, 4: render_sf4,
        6: render_sf6, 7: render_sf7, 8: render_sf8,
    }
    for number in sorted(selected & set(renderers)):
        renderers[number]()
    if 5 in selected:
        render_sf5(
            observed_suffix=args.observed_suffix,
            stem=args.sf5_stem,
        )
    if 9 in selected:
        render_sf9(
            observed_suffix=args.observed_suffix,
            stem=args.sf9_stem,
            ptprc_snapshot=args.ptprc_snapshot,
        )
    if 10 in selected:
        # v3 renders from the tables the portal serves; v2 asserted the superseded
        # 25,621-cell portal embedding and cannot rebuild.
        import render_nar_sf10_v3 as mouse_scont
        mouse_scont.render(
            S.OUTDIR / sf_name(10),
            dpi=450,
            suptitle=(
                sf_title(10, "Cross-species utility using a previously "
                             "published mouse gastrulation dataset")
            ),
            role=(
                "supplementary cross-species utility case using published NGDC "
                "CRA044500 mouse data"
            ),
        )
    print("All NAR SF done →", S.OUTDIR)


if __name__ == "__main__":
    main()
