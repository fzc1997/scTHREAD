#!/usr/bin/env python3
"""Supplementary Figure 10 v3 — cross-species utility on the current portal embedding.

v2 asserted a 25,621-cell / 24-cell-type portal cell map and became
unreproducible when the portal was rebuilt on 30 Jul 2026 onto the full six-run
embedding. v3 renders directly from the files the portal serves
(tables/sf10_v3/, built by scripts/build_sf10_v3_sources.py), so the figure and
the live site cannot drift apart.

Two denominators are kept separate and are labelled on the figure:
  * portal cell map      55,729 cells / 36 cell types / 6 runs (Forward + Reciprocal)
  * ANCHOR DTU analysis  25,621 cells / 24 cell types, forward cross only
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S  # noqa: E402
import render_nar_fig3_mouse_scont as base  # noqa: E402

PROJECT = Path(__file__).resolve().parents[2]
SRC = PROJECT / "tables/sf10_v3"
ROLLING_ATLAS = PROJECT / "tables/fig1_rolling_view_atlas.tsv"

ISO1 = "ENSMUST00000245150"
ISO2 = "ENSMUST00000172812"
WIDTH_MM, HEIGHT_MM = 183.0, 200.0

INK, SLATE, GREY = base.INK, base.SLATE, base.GREY
BLUE, ORANGE, PURPLE = base.BLUE, base.ORANGE, base.PURPLE
LIGHT_GREY = base.LIGHT_GREY
CMAP = LinearSegmentedColormap.from_list(
    "sf10_expr", ["#EDEFF2", "#A7C4C0", "#2F6E6B", "#12233A"]
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_inputs() -> dict[str, object]:
    manifest = json.loads((SRC / "sf10_v3_manifest.json").read_text())
    embedding = manifest["portal_embedding"]
    _require(embedding["cells"] == 55_729, "Portal embedding cell count changed")
    _require(embedding["cell_types"] == 36, "Portal cell-type count changed")
    _require(embedding["runs"] == 6, "Portal run count changed")
    _require(all(manifest["audit_json_agrees"].values()),
             "Portal audit JSON disagrees with the served points table")

    umap = pd.read_parquet(SRC / "malat1_umap.parquet")
    _require(len(umap) == embedding["cells"], "UMAP table does not match the manifest")
    composition = pd.read_csv(SRC / "embedding_composition.tsv", sep="\t")
    celltypes = pd.read_csv(SRC / "celltype_composition.tsv", sep="\t")
    pair = pd.read_csv(SRC / "malat1_iso_pair_by_celltype.tsv", sep="\t")

    # Under the single-cell long-read scope rule this cohort is part of the
    # release, not a rolling-view extra, so its counts come from the release
    # tables. Cells are authoritative per study.
    import render_nar_bio as R
    own_runs = R.registry_done()
    own_runs = own_runs.loc[own_runs["gse"].eq("OWN_ASE_scONT")]
    own_cells = R.study_cells()
    own_cells = own_cells.loc[own_cells["gse"].eq("OWN_ASE_scONT")]
    _require(len(own_cells) == 1, "OWN_ASE_scONT missing from the release")
    _require(own_runs["srr"].nunique() == 6, "Expected six OWN_ASE_scONT runs")
    _require(int(own_cells.iloc[0]["cells"]) == 68_417, "Registered cell count changed")

    # the ANCHOR DTU analysis keeps its own, smaller denominator
    stats_raw = pd.read_csv(base.SAMPLE_STATS, sep="\t")
    stats = dict(zip(stats_raw["Category"], stats_raw["Value"]))
    _require(stats["Total Cells"] == "25,621", "ANCHOR DTU subset cell count changed")
    _require(stats["Cell Types"] == "24", "ANCHOR DTU cell-type count changed")
    lineage = pd.read_csv(base.LINEAGE_DTU, sep="\t")

    return {
        "manifest": manifest,
        "umap": umap,
        "composition": composition,
        "celltypes": celltypes,
        "pair": pair,
        "registered_cells": int(own_cells.iloc[0]["cells"]),
        "lineage": lineage,
    }


def _panel_scatter(ax, umap: pd.DataFrame, values: np.ndarray, title: str, note: str):
    x, y = umap["umap1"].values, umap["umap2"].values
    off = values <= 0
    ax.scatter(x[off], y[off], s=0.30, color="#E4E7EB", linewidths=0, rasterized=True)
    on = ~off
    signal = np.log1p(values[on])
    vmax = float(np.percentile(signal, 99)) if signal.size else 1.0
    order = np.argsort(signal)
    sc = ax.scatter(
        x[on][order], y[on][order], c=signal[order], s=0.45, cmap=CMAP,
        vmin=0.0, vmax=max(vmax, 1e-6), linewidths=0, rasterized=True,
    )
    bar = ax.figure.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    bar.set_label("log(1+molecules)", fontsize=4.8)
    bar.ax.tick_params(labelsize=4.4)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=3)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(LIGHT_GREY); sp.set_linewidth(0.5)
    ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=5.0, color=SLATE)


def render(stem: Path, dpi: int = 450, suptitle: str | None = None,
           role: str | None = None) -> None:
    del dpi  # nar_style.save fixes the export contract at 450 DPI
    data = load_inputs()
    umap: pd.DataFrame = data["umap"]
    manifest = data["manifest"]
    malat1 = manifest["malat1"]

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.18,
                          left=0.06, right=0.98, top=0.90, bottom=0.06)

    # (a) cohort scope: two denominators, stated side by side
    ax = fig.add_subplot(gs[0, 0]); ax.axis("off")
    S.panel_label(ax, "a", x=-0.02, y=1.06)
    ax.set_title("Dataset scope (rolling web view, not the manuscript snapshot)",
                 loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=3)
    comp = data["composition"]
    lines = [
        "CRA044500 · OWN_ASE_scONT · 6 ONT runs",
        f"{data['registered_cells']:,} registered cells",
        "",
        f"Portal cell map: {manifest['portal_embedding']['cells']:,} cells · "
        f"{manifest['portal_embedding']['cell_types']} cell types",
    ]
    for row in comp.itertuples(index=False):
        lines.append(f"  {row.cross:<11s}{row.stage}  {row.cells:>6,} cells  "
                     f"{row.cell_types} types")
    lines += [
        "",
        "ANCHOR DTU analysis (forward cross only):",
        "25,621 cells · 24 cell types",
        "The two denominators are never mixed.",
    ]
    ax.text(0.0, 0.94, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=5.1, color=INK, family="monospace", linespacing=1.5)

    # (b) cell-type composition of the portal embedding
    ax = fig.add_subplot(gs[0, 1])
    S.panel_label(ax, "b", x=-0.02, y=1.06)
    top = data["celltypes"].head(12).iloc[::-1]
    ax.barh(np.arange(len(top)), top["cells"] / 1000.0, color=BLUE, height=0.72)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top["cell_type"], fontsize=5.0)
    ax.tick_params(axis="y", pad=1.5)
    ax.set_xlabel("cells (×10³)", fontsize=6.0)
    ax.set_title("Cell types in the portal embedding (top 12 of 36)",
                 loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=3)
    S.style_ax(ax)

    # (c, d) the two displayed Malat1 isoforms on shared coordinates
    for idx, (iso, col) in enumerate(((ISO1, "iso1_count"), (ISO2, "iso2_count"))):
        ax = fig.add_subplot(gs[1, idx])
        S.panel_label(ax, "cd"[idx], x=-0.02, y=1.06)
        stats = malat1["displayed_isoforms"][iso]
        _panel_scatter(
            ax, umap, umap[col].values.astype(float),
            f"{iso}",
            f"{stats['cells']:,} cells with signal · {int(stats['molecules']):,} molecules",
        )

    # (e) within-gene usage contrast for the two isoforms
    ax = fig.add_subplot(gs[2, 0])
    ax.set_position(ax.get_position().translated(0.055, 0.0).shrunk(0.86, 1.0))
    S.panel_label(ax, "e", x=-0.20, y=1.06)
    pair = data["pair"].set_index("cell_type")
    show = pd.concat([pair.head(6), pair.tail(6)])
    y = np.arange(len(show))
    ax.barh(y - 0.2, show[ISO1], height=0.38, color=BLUE, label=ISO1)
    ax.barh(y + 0.2, show[ISO2], height=0.38, color=ORANGE, label=ISO2)
    ax.set_yticks(y); ax.set_yticklabels(show.index, fontsize=5.2)
    ax.set_xlabel("within-gene usage fraction", fontsize=6.0)
    ax.legend(fontsize=4.8, frameon=False, loc="lower right")
    ax.set_title(
        f"Malat1 isoform usage by cell type (max difference {malat1['max_usage_difference']:.3f})",
        loc="left", fontweight="bold", fontsize=7.2, color=INK, pad=3)
    S.style_ax(ax)

    # (f) evidence boundary
    ax = fig.add_subplot(gs[2, 1]); ax.axis("off")
    S.panel_label(ax, "f", x=-0.02, y=1.06)
    ax.set_title("Evidence boundary", loc="left", fontweight="bold",
                 fontsize=7.2, color=INK, pad=3)
    boundary = (
        "· Descriptive usage proportions from the served portal tables.\n"
        "· Cell-bootstrap pseudo-replicates are not independent embryos, so no\n"
        "  biological-replicate P-value or FDR is reported for this dataset.\n"
        f"· {malat1['cells_with_signal']:,} of "
        f"{manifest['portal_embedding']['cells']:,} cells carry Malat1 signal;\n"
        "  usage fractions use the full gene denominator.\n"
        "· Forward and Reciprocal crosses are shown together but retain their\n"
        "  run labels; they are not pooled as replicates.\n"
        "· This dataset is served by the rolling web view and is outside the\n"
        "  frozen manuscript snapshot."
    )
    ax.text(0.0, 0.92, boundary, transform=ax.transAxes, va="top",
            fontsize=5.4, color=SLATE, linespacing=1.5)

    if suptitle:
        fig.suptitle(suptitle, x=0.06, y=0.965, ha="left", fontsize=9.2,
                     fontweight="bold", color=INK)
    if role:
        fig.text(0.06, 0.015, role, fontsize=5.0, color=SLATE)
    S.save(fig, Path(stem).name)
    print(f"{stem.name} v3 ok  portal={manifest['portal_embedding']['cells']:,} cells")


def main() -> None:
    render(
        S.OUTDIR / "NAR_SF10",
        dpi=450,
        suptitle=("Supplementary Figure 10 · Cross-species utility using a previously "
                  "published mouse gastrulation dataset"),
        role=("supplementary cross-species utility case using published NGDC CRA044500 "
              "mouse data, rendered from the tables served by the portal"),
    )


if __name__ == "__main__":
    main()
