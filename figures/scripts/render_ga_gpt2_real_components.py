#!/usr/bin/env python3
"""Render real-data mini-plots and a composition wireframe for the GPT2 GA.

All quantitative components are derived from frozen release files. The wireframe is
only a layout reference for GPT Image 2 and contains no manuscript-ready artwork.
"""

from __future__ import annotations

import os

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
SC_ROOT = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD"))
REGISTRY_SNAPSHOT = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
UMAP_POINTS = SC_ROOT / "results/paper1/web/atlas_umap_points.parquet"
UMAP_AUDIT = SC_ROOT / "results/paper1/web/atlas_umap_audit.json"
PTPRC_TABLE = ROOT / "tables/Table_S_PTPRC_isoform_usage.tsv"
AXIS_TABLE = ROOT / "tables/Table_S_three_axis_summary.tsv"
DEFAULT_OUTDIR = ROOT / "figures/ga_gpt2_components_v2"
ARIAL = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))

TEAL = "#0B7C86"
BLUE = "#2F6FB2"
PURPLE = "#5B4B9A"
ORANGE = "#D88B1B"
CORAL = "#E34A33"
GREEN = "#4C956C"
INK = "#202124"
MID = "#6D6E73"
LIGHT = "#D7D9DF"
WHITE = "#FFFFFF"
PALE_BLUE = "#EFF6FB"
PALE_TEAL = "#EEF8F8"
PALE_GOLD = "#FFF8E9"
PALE_LAV = "#F4F1FA"

SYSTEM_MAP = {
    "GSE307660": "Blood/immune",
    "GSE276974": "Blood/immune",
    "GSE292324": "Blood/immune",
    "GSE178175": "Neural/sensory",
    "GSE274249": "Neural/sensory",
    "GSE283629": "Neural/sensory",
    "GSE314176": "Neural/sensory",
    "GSE255520": "Neural/sensory",
    "GSE295352": "Neural/sensory",
    "GSE130708": "Neural/sensory",
    "GSE76026": "Neural/sensory",
    "GSE114157": "Neural/sensory",
    "GSE288222": "Heart/vascular",
    "GSE309071": "Heart/vascular",
    "GSE289790": "Cancer",
    "GSE301658": "Cancer",
    "GSE303762": "Cancer",
    "GSE212945": "Cancer",
    "GSE224045": "Cancer",
    "GSE295932": "Cancer",
    "GSE295353": "Endocrine",
    "GSE283658": "Development/embryo",
    "GSE140890": "Development/embryo",
    "GSE185554": "Development/embryo",
    "GSE274527": "Development/embryo",
    "GSE289428": "Development/embryo",
    "GSE248118": "Reproductive",
    "GSE214231": "Other tissues",
    "GSE250381": "Other tissues",
    "GSE252344": "Other tissues",
    "GSE252416": "Other tissues",
}

SYSTEM_COLORS = {
    "Blood/immune": TEAL,
    "Neural/sensory": PURPLE,
    "Cancer": "#D55E00",
    "Endocrine": BLUE,
    "Heart/vascular": "#7A9E3A",
    "Development/embryo": ORANGE,
    "Reproductive": "#CC78BC",
    "Other tissues": "#8B7355",
}

UMAP_COLORS = {
    "Lymphoid": BLUE,
    "Myeloid": TEAL,
    "Erythroid": CORAL,
    "Progenitor": ORANGE,
    "Neural": PURPLE,
    "Cardiovascular": "#7A9E3A",
    "Stromal": "#B06A3C",
    "Other": "#A0A4AB",
}


def setup_style() -> None:
    if ARIAL.exists():
        font_manager.fontManager.addfont(str(ARIAL))
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 450,
            "figure.dpi": 150,
            "savefig.transparent": True,
            "axes.grid": False,
        }
    )


def save_component(fig: plt.Figure, stem: Path, dpi: int = 450) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    svg = stem.with_suffix(".svg")
    png = stem.with_suffix(".png")
    fig.savefig(svg, transparent=True, bbox_inches="tight", pad_inches=0.01)
    fig.savefig(png, transparent=True, bbox_inches="tight", pad_inches=0.01, dpi=dpi)
    plt.close(fig)
    return [svg, png]


def classify_system(row: pd.Series) -> str:
    if row["gse"] in SYSTEM_MAP:
        return SYSTEM_MAP[row["gse"]]
    group = str(row.get("biology_group", ""))
    description = str(row.get("description", "")).lower()
    if group == "组织-分化轨迹":
        return "Development/embryo"
    if group == "脑-神经-翻译":
        return "Neural/sensory"
    if group == "方法-Benchmark":
        if any(key in description for key in ("ovarian", "ovary", "卵巢")):
            return "Reproductive"
        if any(key in description for key in ("tumor", "cancer", "癌", "瘤")):
            return "Cancer"
        return "Other tissues"
    if group == "剪接-APA-调控":
        return "Other tissues"
    if group == "疾病-临床":
        if any(key in description for key in ("marrow", "ccus", "骨髓", "骨髓瘤")):
            return "Blood/immune"
        if any(key in description for key in ("cancer", "glioma", "癌", "瘤")):
            return "Cancer"
        if any(key in description for key in ("heart", "心")):
            return "Heart/vascular"
    return "Other tissues"


def load_catalog_systems(outdir: Path) -> pd.DataFrame:
    registry = pd.read_csv(REGISTRY_SNAPSHOT, sep="\t", dtype=str)
    registry["isoquant_cells"] = pd.to_numeric(
        registry["isoquant_cells"], errors="coerce"
    ).fillna(0)
    done = registry[
        registry["isoquant_status"].fillna("").str.lower().eq("done")
    ].copy()
    if len(done) != 469 or int(done["isoquant_cells"].sum()) != 845_781:
        raise RuntimeError(
            "Frozen registry snapshot no longer matches 469 / 845,781."
        )
    done["system"] = done.apply(classify_system, axis=1)
    study_audit = (
        done.groupby(["gse", "system"], as_index=False, dropna=False)
        .agg(
            cells=("isoquant_cells", "sum"),
            samples=("srr", "nunique"),
            species=("species", lambda x: "|".join(sorted(set(x.dropna())))),
            description=(
                "description",
                lambda x: "|".join(sorted(v for v in set(x.dropna()) if v)),
            ),
        )
        .sort_values(["system", "cells"], ascending=[True, False])
    )
    study_audit.to_csv(
        outdir / "catalog_study_biology_classification_20260726.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "display_label": "Human ovary",
                "registry_id": "benagen_ovary_Young2",
                "species": "human",
                "assay": "scONT",
                "stage_or_tissue": "snap-frozen ovary",
                "included_in_469_release_denominator": "no",
            },
            {
                "display_label": "Mouse gastrula",
                "registry_id": "OWN_ASE_scONT",
                "species": "mouse",
                "assay": "scONT",
                "stage_or_tissue": "E6.5/E7.5/E8.5 embryo",
                "included_in_469_release_denominator": "no",
            },
        ]
    ).to_csv(
        outdir / "featured_project_cohorts_20260726.tsv",
        sep="\t",
        index=False,
    )
    summary = (
        done.groupby("system", as_index=False)
        .agg(
            cells=("isoquant_cells", "sum"),
            samples=("srr", "nunique"),
            studies=("gse", "nunique"),
        )
        .sort_values("cells", ascending=False)
    )
    return summary


def plot_catalog_systems(summary: pd.DataFrame, outdir: Path) -> list[Path]:
    summary.to_csv(outdir / "catalog_system_composition_20260726.tsv", sep="\t", index=False)
    display = summary.sort_values("cells", ascending=True).copy()
    fig, ax = plt.subplots(figsize=(43 / 25.4, 28 / 25.4), facecolor="none")
    y = np.arange(len(display))
    values = display["cells"].to_numpy() / 1000
    ax.barh(
        y,
        values,
        color=[SYSTEM_COLORS.get(v, MID) for v in display["system"]],
        height=0.68,
        edgecolor=WHITE,
        linewidth=0.35,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(display["system"], fontsize=4.5)
    ax.set_xlabel("Cells (×10³)", fontsize=4.5)
    ax.tick_params(axis="x", labelsize=4, width=0.35, length=2)
    ax.tick_params(axis="y", width=0, length=0, pad=1.5)
    xmax = max(values) * 1.30
    ax.set_xlim(0, xmax)
    for yi, (_, row) in enumerate(display.iterrows()):
        ax.text(
            values[yi] + xmax * 0.02,
            yi,
            f"{int(row.cells):,}",
            va="center",
            ha="left",
            fontsize=4.0,
            color=MID,
        )
    ax.set_title(
        "Catalog cells by system",
        loc="left",
        fontsize=6.2,
        fontweight="bold",
        pad=2,
        color=INK,
    )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(LIGHT)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.grid(axis="x", color="#ECEDEF", linewidth=0.35)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.30, right=0.98, top=0.84, bottom=0.20)
    return save_component(fig, outdir / "catalog_system_composition")


def plot_ptprc_switch(outdir: Path) -> list[Path]:
    table = pd.read_csv(PTPRC_TABLE, sep="\t")
    target = ["ENST00000367364", "ENST00000697630"]
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
    sub = table[
        table["transcript_id"].isin(target) & table["ct"].isin(cell_order)
    ].copy()
    matrix = (
        sub.pivot(index="transcript_id", columns="ct", values="frac")
        .reindex(index=target, columns=cell_order)
        .fillna(0)
    )
    source = (
        matrix.rename_axis("transcript_id")
        .reset_index()
        .melt(id_vars="transcript_id", var_name="cell_type", value_name="fraction")
    )
    source.to_csv(outdir / "ptprc_two_isoform_switch.tsv", sep="\t", index=False)

    cmap = LinearSegmentedColormap.from_list(
        "ptprc", ["#F7FAFC", "#BFDCEB", "#4E91BC", "#175B8D"]
    )
    fig, ax = plt.subplots(figsize=(47 / 25.4, 18 / 25.4), facecolor="none")
    im = ax.imshow(matrix.to_numpy(), cmap=cmap, vmin=0, vmax=0.60, aspect="auto")
    short_ct = ["B", "DC", "Ery", "Mono", "NK", "Plasma", "Prog", "T"]
    ax.set_xticks(range(len(short_ct)))
    ax.set_xticklabels(short_ct, fontsize=4.2)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["ENST…67364", "ENST…97630"], fontsize=4.2)
    ax.tick_params(length=0, pad=1.2)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            ax.text(
                j,
                i,
                f"{value * 100:.0f}",
                ha="center",
                va="center",
                fontsize=3.7,
                color=WHITE if value > 0.34 else INK,
                fontweight="bold" if value > 0.30 else "normal",
            )
    ax.set_title(
        "PTPRC isoform fraction (%)",
        loc="left",
        fontsize=6.1,
        fontweight="bold",
        color=INK,
        pad=2,
    )
    for spine in ax.spines.values():
        spine.set_color(WHITE)
        spine.set_linewidth(0.5)
    fig.subplots_adjust(left=0.21, right=0.995, top=0.78, bottom=0.23)
    outputs = save_component(fig, outdir / "ptprc_isoform_switch")
    return outputs


def broad_cell_class(label: str) -> str:
    text = label.lower()
    if any(k in text for k in ("t cell", "nk", "b cell", "plasma", "lymphoid")):
        return "Lymphoid"
    if any(k in text for k in ("monocyte", "myeloid", "dendritic")):
        return "Myeloid"
    if "erythroid" in text:
        return "Erythroid"
    if any(k in text for k in ("progenitor", "cycling", "pluripotent")):
        return "Progenitor"
    if any(
        k in text
        for k in (
            "neuron",
            "neural",
            "astro",
            "oligodend",
            "opc",
            "retinal",
            "amacrine",
            "bipolar",
            "muller",
            "rod",
            "cajal",
            "hipp_",
            "npc-like",
            "mes-like",
        )
    ):
        return "Neural"
    if any(
        k in text
        for k in (
            "cardiomyocyte",
            "endothelial",
            "pericyte",
            "smooth muscle",
            "endocardial",
            "vascular",
            "lymphatic",
        )
    ):
        return "Cardiovascular"
    if "fibroblast" in text:
        return "Stromal"
    return "Other"


def plot_real_umap(outdir: Path) -> list[Path]:
    points = pd.read_parquet(
        UMAP_POINTS, columns=["cell_id", "cell_type", "umap1", "umap2"]
    )
    sampled = pd.concat(
        [
            frame.sample(n=min(140, len(frame)), random_state=42)
            for _, frame in points.groupby("cell_type", observed=True)
        ],
        ignore_index=True,
    )
    sampled["broad_class"] = sampled["cell_type"].map(broad_cell_class)
    sampled.to_csv(outdir / "atlas_umap_stratified_sample.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(43 / 25.4, 31 / 25.4), facecolor="none")
    for label in UMAP_COLORS:
        sub = sampled[sampled["broad_class"] == label]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap1"],
            sub["umap2"],
            s=2.8,
            c=UMAP_COLORS[label],
            alpha=0.76,
            linewidths=0,
            rasterized=False,
            label=label,
        )
    ax.set_title(
        "Harmonized cell-type atlas",
        loc="left",
        fontsize=6.2,
        fontweight="bold",
        color=INK,
        pad=1.5,
    )
    ax.text(
        0.0,
        -0.04,
        "74,906 portal cells · exploratory UMAP",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=4.2,
        color=MID,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.08)
    return save_component(fig, outdir / "atlas_umap_real")


def plot_axis_inventory(outdir: Path) -> list[Path]:
    table = pd.read_csv(AXIS_TABLE, sep="\t")
    order = ["DIU", "APA", "ASE"]
    table = table.set_index("axis").reindex(order).reset_index()
    table["sig_fraction"] = (
        table["n_sig_fdr05_effect_gate"] / table["n_genes_tested"]
    )
    table.to_csv(outdir / "three_axis_inventory_source.tsv", sep="\t", index=False)

    colors = {"DIU": BLUE, "APA": ORANGE, "ASE": TEAL}
    fig, ax = plt.subplots(figsize=(42 / 25.4, 20 / 25.4), facecolor="none")
    y = np.arange(len(table))
    ax.barh(
        y,
        table["sig_fraction"] * 100,
        color=[colors[a] for a in table["axis"]],
        height=0.56,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(table["axis"], fontsize=5.0, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 42)
    ax.set_xticks([0, 20, 40])
    ax.set_xticklabels(["0", "20", "40%"], fontsize=4.0)
    ax.tick_params(length=0, pad=1.2)
    for yi, row in table.iterrows():
        ax.text(
            max(row["sig_fraction"] * 100 + 1.0, 1.0),
            yi,
            f"{int(row['n_sig_fdr05_effect_gate']):,} / {int(row['n_genes_tested']):,}",
            va="center",
            ha="left",
            fontsize=4.2,
            color=INK,
        )
    ax.set_title(
        "Precomputed cell-type maps",
        loc="left",
        fontsize=6.2,
        fontweight="bold",
        color=INK,
        pad=1.5,
    )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(LIGHT)
    ax.spines["bottom"].set_linewidth(0.4)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.80, bottom=0.22)
    return save_component(fig, outdir / "three_axis_inventory")


def rounded(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face: str = WHITE,
    edge: str = PURPLE,
    lw: float = 1.0,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.2,rounding_size=2.2",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
        )
    )


def plot_wireframe(outdir: Path) -> list[Path]:
    width, height = 183.0, 78.0
    fig = plt.figure(figsize=(width / 25.4, height / 25.4), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), width, height, facecolor=WHITE, edgecolor="none"))

    ax.text(
        width / 2,
        74.2,
        "TITLE + URL",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )
    rounded(ax, 3, 37, 42, 31, face=PALE_TEAL)
    rounded(ax, 47, 37, 31, 31, face="#FFF1E8")
    rounded(ax, 80, 37, 50, 31, face=PALE_LAV)
    rounded(ax, 132, 8, 48, 60, face=PALE_BLUE)
    rounded(ax, 3, 8, 127, 27, face=WHITE)

    ax.text(24, 64.5, "CATALOG COVERAGE", ha="center", fontsize=7, fontweight="bold")
    ax.text(62.5, 64.5, "DATABASE CYLINDER", ha="center", fontsize=7, fontweight="bold")
    ax.text(105, 64.5, "4 RNA LAYER CARDS", ha="center", fontsize=7, fontweight="bold")
    ax.text(156, 64.5, "QUERY + ONLINE ANALYSIS", ha="center", fontsize=7, fontweight="bold")

    for x, y, w, h, label in [
        (5.5, 39.5, 37, 14, "REAL SYSTEM COMPOSITION"),
        (78.5, 11.0, 24, 17, "REAL PTPRC HEATMAP"),
        (135.5, 27.0, 41, 23, "REAL UMAP"),
        (104.5, 11.0, 23, 17, "REAL DIU / APA / ASE"),
    ]:
        ax.add_patch(
            Rectangle(
                (x, y),
                w,
                h,
                facecolor=WHITE,
                edgecolor="#80868B",
                linewidth=0.9,
                linestyle="--",
            )
        )
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=5.5, color=MID)

    for x, label in zip(
        [15, 39, 64, 90, 116],
        ["BROWSE", "SEARCH", "DOWNLOAD", "GENE DETAIL", "PRECOMPUTED"],
    ):
        ax.text(x, 31.5, label, ha="center", va="center", fontsize=6.0, fontweight="bold")

    ax.text(
        width / 2,
        3.7,
        "FOOTER",
        ha="center",
        va="center",
        fontsize=6.5,
        fontweight="bold",
        color=TEAL,
    )
    stem = outdir / "ga_gpt2_layout_wireframe"
    svg = stem.with_suffix(".svg")
    png = stem.with_suffix(".png")
    fig.savefig(svg, transparent=False, bbox_inches=None)
    fig.savefig(png, transparent=False, bbox_inches=None, dpi=240)
    plt.close(fig)
    return [svg, png]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    setup_style()
    outputs: list[Path] = []
    systems = load_catalog_systems(outdir)
    outputs.extend(plot_catalog_systems(systems, outdir))
    outputs.extend(plot_ptprc_switch(outdir))
    outputs.extend(plot_real_umap(outdir))
    outputs.extend(plot_axis_inventory(outdir))
    outputs.extend(plot_wireframe(outdir))
    for output in outputs:
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
