#!/usr/bin/env python3
"""Render an OWN_ASE_scONT mouse-gastrulation alternative for NAR Figure 3.

Figure contract
---------------
Core conclusion
    The reciprocal-cross mouse gastrulation scONT cohort provides a
    cell-lineage-resolved use case spanning developmental isoform usage and
    gene-level allelic-ratio evidence.
Evidence chain
    a, project-cohort provenance and denominator;
    b, cell-type composition of the DTU analysis subset;
    c, stage-associated isoform-usage screens in three represented lineages;
    d, a traceable two-isoform Tnrc6c developmental switch;
    e, five E6.5 genes from a cell-type allelic-ratio screen;
    f, analysis-layer inventory and export/source traceability.
Guardrails
    The 68,417-cell project registry count and the 25,621-cell DTU analysis
    subset are distinct denominators. Isoform trajectories are descriptive
    proportions, not proof of mechanism. Allele-ratio heterogeneity is not
    equivalent to parent-of-origin imprinting. Cell-bootstrap pseudo-replicates
    are not independent embryos, and the allele panel is not isoform-level
    AS-DTU.

The renderer uses only existing project tables and writes a compact source-data
bundle beside the manuscript tables. It does not overwrite the preferred human
PTPRC Figure 3.
"""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyBboxPatch
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
TABLES = PROJECT / "tables"
ASE_ROOT = Path(
    os.environ.get("SCTHREAD_ISOFORM_ROOT", "/gpfs/home/fuzc/project/ASE/cc/config/ISOform") + "/"
    "06_isoform_identification/dtu_analysis"
)

STUDY_ATLAS = TABLES / "fig1_study_atlas.tsv"
# OWN_ASE_scONT is a rolling-web-view dataset, deliberately outside the frozen
# manuscript snapshot, so it is described by the rolling-view atlas.
ROLLING_ATLAS = TABLES / "fig1_rolling_view_atlas.tsv"
SAMPLE_STATS = ASE_ROOT / "tables/Table1_sample_statistics.tsv"
CELL_COMPOSITION = ASE_ROOT / "tables/ST1_celltype_annotation.tsv"
LINEAGE_DTU = ASE_ROOT / "results/developmental_DTU_multi_celltype_summary.tsv"
TEMPORAL_CASES = ASE_ROOT / "results/nature_cases_developmental.tsv"
AS_DTU = ASE_ROOT / "results/E65_LR_allele_specific_DTU_v2.tsv"
MOUSE_GENES = Path(
    os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD") + "/web/data/mouse_gene_intervals_grcm39.tsv"
)

WIDTH_MM = 183.0
HEIGHT_MM = 170.0

INK = "#1C2130"
SLATE = "#4A5060"
BLUE = "#0173B2"
ORANGE = "#DE8F05"
GREEN = "#029E73"
PURPLE = "#9C4E97"
VERMILLION = "#CC3311"
SKY = "#56B4E9"
GREY = "#8A8F98"
LIGHT_GREY = "#D8DDE2"
PALE_BLUE = "#E8F2F8"
PALE_ORANGE = "#FBF0DE"
PALE_GREEN = "#E6F4EF"
PALE_PURPLE = "#F3EAF3"
WHITE = "#FFFFFF"

STAGES = ["E6.5", "E7.5", "E8.5"]
STAGE_COLORS = {"E6.5": BLUE, "E7.5": ORANGE, "E8.5": GREEN}
DISPLAY_LINEAGES = [
    "Neuroectoderm",
    "Visceral_Endoderm_Yolk_Sac",
    "Embryonic_Epiblast",
    "Primitive_Streak_Epiblast",
    "Presomitic_Mesoderm",
    "Streak_Mesoderm",
    "Early_Erythrocytes",
]
LINEAGE_LABELS = {
    "Neuroectoderm": "Neuroectoderm",
    "Visceral_Endoderm_Yolk_Sac": "Visceral endoderm/\nyolk sac",
    "Embryonic_Epiblast": "Embryonic epiblast",
    "Primitive_Streak_Epiblast": "Primitive-streak\nepiblast",
    "Presomitic_Mesoderm": "Presomitic mesoderm",
    "Streak_Mesoderm": "Streak mesoderm",
    "Early_Erythrocytes": "Early erythrocytes",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _configure_style() -> str:
    arial = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))
    if arial.is_file():
        fm.fontManager.addfont(str(arial))
        family = fm.FontProperties(fname=str(arial)).get_name()
    else:
        family = "DejaVu Sans"
    plt.rcParams.update(
        {
            "font.family": family,
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 450,
            "savefig.transparent": True,
            "figure.facecolor": "none",
            "axes.facecolor": "none",
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )
    return family


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _load_inputs() -> dict[str, object]:
    for path in (
        STUDY_ATLAS,
        SAMPLE_STATS,
        CELL_COMPOSITION,
        LINEAGE_DTU,
        TEMPORAL_CASES,
        AS_DTU,
        MOUSE_GENES,
    ):
        _require(path.is_file(), f"Missing input: {path}")

    atlas = pd.read_csv(STUDY_ATLAS, sep="\t")
    own = atlas.loc[atlas["gse"].eq("OWN_ASE_scONT")].copy()
    if own.empty and ROLLING_ATLAS.is_file():
        rolling = pd.read_csv(ROLLING_ATLAS, sep="\t")
        own = rolling.loc[rolling["gse"].eq("OWN_ASE_scONT")].copy()
    _require(len(own) == 1, "OWN_ASE_scONT must occur once in the study atlas")
    own_row = own.iloc[0]
    _require(int(own_row["n_runs"]) == 6, "Expected six OWN_ASE_scONT runs")
    _require(int(own_row["n_cells"]) == 68_417, "Unexpected project cell count")

    stats_raw = pd.read_csv(SAMPLE_STATS, sep="\t")
    stats = dict(zip(stats_raw["Category"], stats_raw["Value"]))
    _require(stats["Total Cells"] == "25,621", "Unexpected DTU subset cell count")
    _require(stats["Cell Types"] == "24", "Unexpected cell-type count")
    _require(stats["Genes with Isoforms"] == "47,972", "Unexpected gene count")
    _require(stats["Total Isoforms"] == "182,689", "Unexpected isoform count")

    composition = pd.read_csv(CELL_COMPOSITION, sep="\t")
    _require(
        int(composition["Total"].sum()) == 25_621,
        "Cell composition does not sum to the analysis denominator",
    )
    _require(set(STAGES).issubset(composition.columns), "Missing stage columns")

    lineage = pd.read_csv(LINEAGE_DTU, sep="\t")
    expected_lineages = {
        "Neuroectoderm",
        "Presomitic_Mesoderm",
        "Visceral_Endoderm_Yolk_Sac",
    }
    _require(
        set(lineage["celltype"]) == expected_lineages,
        "Unexpected lineage-DTU table content",
    )

    temporal = pd.read_csv(TEMPORAL_CASES, sep="\t")
    tnrc6c = temporal.loc[temporal["gene_name"].eq("Tnrc6c")].copy()
    _require(len(tnrc6c) == 2, "Expected two curated Tnrc6c isoforms")
    _require(
        set(tnrc6c["trend"]) == {"increasing", "decreasing"},
        "Tnrc6c must contain reciprocal trends",
    )

    as_dtu = pd.read_csv(AS_DTU, sep="\t")
    _require(len(as_dtu) == 5, "Expected five significant AS-DTU genes")
    _require((as_dtu["padj"] < 0.05).all(), "AS-DTU source contains non-significant rows")
    genes = pd.read_csv(MOUSE_GENES, sep="\t", usecols=["gid", "gname"])
    gene_map = (
        genes.loc[genes["gid"].isin(as_dtu["gene_id"])]
        .drop_duplicates("gid")
        .set_index("gid")["gname"]
    )
    as_dtu["gene_name"] = as_dtu["gene_id"].map(gene_map)
    _require(as_dtu["gene_name"].notna().all(), "Missing mouse gene symbols")

    return {
        "own": own_row,
        "stats": stats,
        "composition": composition,
        "lineage": lineage,
        "tnrc6c": tnrc6c,
        "as_dtu": as_dtu,
    }


def _write_source_bundle(data: dict[str, object]) -> list[Path]:
    TABLES.mkdir(parents=True, exist_ok=True)
    composition = data["composition"].copy()
    lineage = data["lineage"].copy()
    tnrc6c = data["tnrc6c"].copy()
    as_dtu = data["as_dtu"].copy()

    outputs = [
        TABLES / "Fig3_mouse_scONT_cell_composition.tsv",
        TABLES / "Fig3_mouse_scONT_lineage_DTU.tsv",
        TABLES / "Fig3_mouse_scONT_Tnrc6c_switch.tsv",
        TABLES / "Fig3_mouse_scONT_AS_DTU.tsv",
    ]
    composition.to_csv(outputs[0], sep="\t", index=False)
    lineage.to_csv(outputs[1], sep="\t", index=False)
    tnrc6c.to_csv(outputs[2], sep="\t", index=False)
    as_dtu.to_csv(outputs[3], sep="\t", index=False)

    manifest = {
        "figure": "NAR_Fig3_mouse_scONT_v1",
        "role": "alternative mouse OWN_ASE_scONT use case; does not replace preferred PTPRC Figure 3",
        "project_cohort": {
            "study_id": "OWN_ASE_scONT",
            "runs": 6,
            "registered_cells": 68417,
            "stages": STAGES,
            "platform": "ONT",
            "species": "mouse",
        },
        "analysis_subset": {
            "cells": 25621,
            "cell_types": 24,
            "genes_with_isoforms": 47972,
            "isoforms": 182689,
        },
        "guardrails": [
            "68,417 registered project cells and 25,621 DTU-analysis cells are distinct denominators",
            "isoform trajectories are descriptive proportions and do not establish mechanism",
            "developmental DRIMSeq uses cell-bootstrap pseudo-replicates, not independent embryos",
            "cell-type allele-ratio heterogeneity is a gene-level max-versus-min screen, not isoform-level AS-DTU",
            "cell-type allele-ratio heterogeneity is not equivalent to imprinting",
        ],
        "inputs": {
            str(path): {"sha256": _sha256(path)}
            for path in (
                STUDY_ATLAS,
                SAMPLE_STATS,
                CELL_COMPOSITION,
                LINEAGE_DTU,
                TEMPORAL_CASES,
                AS_DTU,
                MOUSE_GENES,
            )
        },
        "source_tables": [str(path) for path in outputs],
    }
    manifest_path = TABLES / "Fig3_mouse_scONT_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    outputs.append(manifest_path)
    return outputs


def _panel_label(ax: plt.Axes, letter: str, x: float = -0.12) -> None:
    ax.text(
        x,
        1.13,
        letter,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        style="normal",
        ha="right",
        va="top",
        clip_on=False,
    )


def _panel_title(ax: plt.Axes, title: str) -> None:
    ax.set_title(title, loc="left", fontsize=7, fontweight="bold", y=1.105, pad=0)


def _rounded(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    facecolor: str,
    edgecolor: str = LIGHT_GREY,
    linewidth: float = 0.6,
    radius: float = 0.025,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )
    ax.add_patch(patch)
    return patch


def _draw_mouse(ax: plt.Axes, center: tuple[float, float], scale: float) -> None:
    x, y = center
    ax.add_patch(Circle((x, y), 0.072 * scale, facecolor=SLATE, edgecolor=SLATE))
    ax.add_patch(Circle((x - 0.048 * scale, y + 0.058 * scale), 0.033 * scale, facecolor=SLATE, edgecolor=SLATE))
    ax.add_patch(Circle((x + 0.048 * scale, y + 0.058 * scale), 0.033 * scale, facecolor=SLATE, edgecolor=SLATE))
    ax.add_patch(Circle((x - 0.025 * scale, y + 0.012 * scale), 0.008 * scale, facecolor=WHITE, edgecolor="none"))
    ax.add_patch(Circle((x + 0.025 * scale, y + 0.012 * scale), 0.008 * scale, facecolor=WHITE, edgecolor="none"))
    ax.add_patch(Circle((x, y - 0.025 * scale), 0.007 * scale, facecolor="#E7A0A0", edgecolor="none"))
    ax.plot(
        [x - 0.010 * scale, x - 0.075 * scale],
        [y - 0.018 * scale, y - 0.035 * scale],
        color=SLATE,
        linewidth=0.55,
    )
    ax.plot(
        [x + 0.010 * scale, x + 0.075 * scale],
        [y - 0.018 * scale, y - 0.035 * scale],
        color=SLATE,
        linewidth=0.55,
    )


def _draw_cohort(ax: plt.Axes, data: dict[str, object]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _panel_label(ax, "a", x=-0.035)
    _panel_title(ax, "OWN_ASE_scONT mouse gastrulation cohort")
    _rounded(ax, (0.0, 0.04), 1.0, 0.82, PALE_BLUE, edgecolor="#BDD3E1")
    _draw_mouse(ax, (0.12, 0.52), 0.90)
    ax.text(
        0.27,
        0.70,
        "Reciprocal-cross mouse embryo · scONT",
        fontsize=7,
        fontweight="bold",
        ha="left",
        va="center",
    )
    ax.text(
        0.27,
        0.56,
        "6 runs · 68,417 registered project cells",
        fontsize=6.2,
        color=SLATE,
        ha="left",
        va="center",
    )
    x_positions = [0.52, 0.70, 0.88]
    for x, stage in zip(x_positions, STAGES):
        ax.add_patch(Circle((x, 0.35), 0.065, facecolor=STAGE_COLORS[stage], edgecolor=WHITE, linewidth=1.0))
        ax.text(x, 0.35, stage, fontsize=5.7, fontweight="bold", color=WHITE, ha="center", va="center")
    ax.plot([0.455, 0.945], [0.35, 0.35], color="#A8B0BA", linewidth=0.8, zorder=0)
    ax.text(
        0.23,
        0.19,
        "Project registry denominator",
        fontsize=5.3,
        fontweight="bold",
        color=BLUE,
        ha="left",
        va="center",
    )
    ax.text(
        0.48,
        0.19,
        "≠",
        fontsize=7,
        fontweight="bold",
        color=VERMILLION,
        ha="center",
        va="center",
    )
    ax.text(
        0.52,
        0.19,
        "25,621-cell DTU analysis subset",
        fontsize=5.3,
        fontweight="bold",
        color=GREEN,
        ha="left",
        va="center",
    )


def _draw_composition(ax: plt.Axes, data: dict[str, object]) -> None:
    composition = data["composition"].set_index("Cell_Type")
    selected = composition.reindex(DISPLAY_LINEAGES)[STAGES].copy()
    other = composition[STAGES].sum() - selected.sum()
    selected.loc["Other"] = other
    shares = selected.div(selected.sum(axis=0), axis=1) * 100.0

    colors = [BLUE, SKY, PURPLE, "#B07AA1", ORANGE, "#D4A72C", VERMILLION, LIGHT_GREY]
    bottom = np.zeros(len(STAGES))
    for lineage, color in zip(shares.index, colors):
        values = shares.loc[lineage].to_numpy()
        ax.bar(
            np.arange(len(STAGES)),
            values,
            bottom=bottom,
            width=0.62,
            color=color,
            edgecolor=WHITE,
            linewidth=0.35,
            label=LINEAGE_LABELS.get(lineage, lineage),
        )
        bottom += values
    ax.set_ylim(0, 100)
    ax.set_xticks(np.arange(len(STAGES)), STAGES)
    ax.set_ylabel("Cells (%)")
    ax.set_yticks([0, 25, 50, 75, 100])
    _panel_label(ax, "b")
    _panel_title(ax, "Cellular context of the DTU subset")
    ax.text(
        0.0,
        1.015,
        "25,621 cells · 24 cell types",
        transform=ax.transAxes,
        fontsize=5.2,
        color=SLATE,
        ha="left",
        va="bottom",
    )
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.03, 1.01),
        fontsize=4.25,
        ncol=1,
        handlelength=1.0,
        labelspacing=0.45,
        borderaxespad=0,
    )


def _draw_lineage_dtu(ax: plt.Axes, data: dict[str, object]) -> None:
    lineage = data["lineage"].set_index("celltype").loc[
        ["Neuroectoderm", "Presomitic_Mesoderm", "Visceral_Endoderm_Yolk_Sac"]
    ]
    labels = ["Neuroectoderm", "Presomitic\nmesoderm", "Visceral endoderm/\nyolk sac"]
    values = lineage["sig_genes"].to_numpy()
    isoforms = lineage["sig_isoforms"].to_numpy()
    colors = [BLUE, ORANGE, GREEN]
    bars = ax.bar(np.arange(3), values, width=0.54, color=colors, edgecolor=INK, linewidth=0.35)
    for bar, genes, iso in zip(bars, values, isoforms):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            genes + 28,
            f"{genes:,} genes\n{iso:,} isoform–comparison hits",
            fontsize=4.7,
            ha="center",
            va="bottom",
        )
    ax.set_xticks(np.arange(3), labels)
    ax.set_ylabel("Genes with ≥1 significant isoform")
    ax.set_ylim(0, 1020)
    ax.set_yticks([0, 250, 500, 750, 1000])
    _panel_label(ax, "c")
    _panel_title(ax, "Stage-associated isoform shifts across lineages")
    ax.text(
        0.0,
        1.015,
        "Cell-bootstrap pseudo-replicates · feature-level BH q < 0.05",
        transform=ax.transAxes,
        fontsize=5.2,
        color=SLATE,
        ha="left",
        va="bottom",
    )


def _draw_tnrc6c(ax: plt.Axes, data: dict[str, object]) -> None:
    tnrc6c = data["tnrc6c"].sort_values("trend")
    colors = {"increasing": BLUE, "decreasing": ORANGE}
    for _, row in tnrc6c.iterrows():
        values = [row["prop_E65"], row["prop_E75"], row["prop_E85"]]
        short_id = str(row["isoform_id"])[-6:]
        ax.plot(
            np.arange(3),
            values,
            marker="o",
            markersize=4,
            linewidth=1.4,
            color=colors[row["trend"]],
            label=f"…{short_id} ({row['trend']})",
        )
    ax.set_xlim(-0.15, 2.15)
    ax.set_ylim(-0.04, 1.04)
    ax.set_xticks(np.arange(3), STAGES)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Isoform usage fraction")
    _panel_label(ax, "d")
    _panel_title(ax, "Tnrc6c reciprocal isoform switch")
    ax.text(
        0.0,
        1.015,
        "Primitive-streak epiblast · descriptive proportions",
        transform=ax.transAxes,
        fontsize=5.2,
        color=SLATE,
        ha="left",
        va="bottom",
    )
    ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.52), fontsize=5.0)
    ax.text(
        0.99,
        0.84,
        "Δ usage = ±0.996",
        transform=ax.transAxes,
        fontsize=5.2,
        fontweight="bold",
        color=PURPLE,
        ha="right",
        va="center",
    )


def _draw_as_dtu(ax: plt.Axes, data: dict[str, object]) -> None:
    as_dtu = data["as_dtu"].sort_values("ratio_range")
    y = np.arange(len(as_dtu))
    colors = np.where(as_dtu["ratio_range"].to_numpy() >= 0.4, VERMILLION, PURPLE)
    ax.hlines(y, 0, as_dtu["ratio_range"], color=LIGHT_GREY, linewidth=1.5)
    ax.scatter(
        as_dtu["ratio_range"],
        y,
        s=28,
        c=colors,
        edgecolor=INK,
        linewidth=0.35,
        zorder=3,
    )
    for yy, (_, row) in zip(y, as_dtu.iterrows()):
        ax.text(
            row["ratio_range"] + 0.012,
            yy,
            f"q={row['padj']:.1e}",
            fontsize=4.5,
            color=SLATE,
            ha="left",
            va="center",
        )
    ax.set_yticks(y, as_dtu["gene_name"])
    ax.set_xlim(0, 0.50)
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xlabel("Allelic-ratio range across cell types")
    _panel_label(ax, "e")
    _panel_title(ax, "Allelic-ratio heterogeneity at E6.5")
    ax.text(
        0.0,
        1.01,
        "Gene-level max–min screen · BH q < 0.05 · not isoform-level AS-DTU",
        transform=ax.transAxes,
        fontsize=5.2,
        color=SLATE,
        ha="left",
        va="bottom",
    )


def _draw_access(ax: plt.Axes, data: dict[str, object]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _panel_label(ax, "f")
    _panel_title(ax, "Traceable evidence layers and export")
    _rounded(ax, (0.0, 0.02), 1.0, 0.84, "#222936", edgecolor="#222936")
    layers = [
        ("CELL ATLAS", "25,621 cells · 24 types", BLUE),
        ("ISOFORMS", "182,689 transcripts", ORANGE),
        ("STAGE SHIFTS", "feature-level screens", GREEN),
        ("ALLELIC RATIO", "gene-level E6.5 screen", PURPLE),
    ]
    for y, (label, value, color) in zip([0.70, 0.54, 0.38, 0.22], layers):
        ax.add_patch(Circle((0.09, y), 0.028, facecolor=color, edgecolor=WHITE, linewidth=0.5))
        ax.text(0.15, y + 0.025, label, fontsize=4.5, fontweight="bold", color="#AEB8C6", ha="left", va="center")
        ax.text(0.15, y - 0.027, value, fontsize=5.1, color=WHITE, ha="left", va="center")
    _rounded(ax, (0.59, 0.18), 0.34, 0.56, "#2C3545", edgecolor="#59677B")
    ax.text(0.76, 0.65, "SOURCE DATA", fontsize=4.8, fontweight="bold", color="#B7D8D3", ha="center")
    ax.text(0.76, 0.54, "4 TSV + manifest", fontsize=5.4, fontweight="bold", color=WHITE, ha="center")
    ax.text(0.76, 0.42, "input SHA-256", fontsize=4.6, color="#AEB8C6", ha="center")
    ax.text(0.76, 0.32, "editable SVG/PDF", fontsize=4.6, color="#AEB8C6", ha="center")
    ax.text(0.76, 0.22, "450-dpi preview", fontsize=4.6, color="#AEB8C6", ha="center")
    ax.text(
        0.05,
        0.085,
        "Alternative mouse use case\nPreferred PTPRC Figure 3 remains unchanged",
        fontsize=4.05,
        fontweight="bold",
        color="#8FD1C6",
        ha="left",
        va="center",
    )


def render(stem: Path, dpi: int) -> tuple[list[Path], list[Path]]:
    family = _configure_style()
    data = _load_inputs()
    source_outputs = _write_source_bundle(data)

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="none")
    ax_a = fig.add_axes([0.065, 0.765, 0.91, 0.175])
    ax_b = fig.add_axes([0.070, 0.425, 0.245, 0.245])
    ax_c = fig.add_axes([0.565, 0.425, 0.375, 0.245])
    ax_d = fig.add_axes([0.070, 0.085, 0.245, 0.235])
    ax_e = fig.add_axes([0.395, 0.085, 0.245, 0.235])
    ax_f = fig.add_axes([0.720, 0.085, 0.245, 0.235])

    _draw_cohort(ax_a, data)
    _draw_composition(ax_b, data)
    _draw_lineage_dtu(ax_c, data)
    _draw_tnrc6c(ax_d, data)
    _draw_as_dtu(ax_e, data)
    _draw_access(ax_f, data)

    for ax in (ax_b, ax_c, ax_d, ax_e):
        ax.tick_params(length=2.0, width=0.5, pad=1.5)
        for spine in ax.spines.values():
            spine.set_color(INK)
            spine.set_linewidth(0.5)

    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [
        stem.with_suffix(".pdf"),
        stem.with_suffix(".svg"),
        stem.with_suffix(".png"),
    ]
    fig.savefig(outputs[0], format="pdf", transparent=True, bbox_inches=None, dpi=dpi)
    fig.savefig(outputs[1], format="svg", transparent=True, bbox_inches=None, dpi=dpi)
    fig.savefig(outputs[2], format="png", transparent=True, bbox_inches=None, dpi=dpi)
    plt.close(fig)

    print(
        json.dumps(
            {
                "font_family": family,
                "canvas_mm": [WIDTH_MM, HEIGHT_MM],
                "outputs": [str(path) for path in outputs],
                "source_outputs": [str(path) for path in source_outputs],
            },
            indent=2,
        )
    )
    return outputs, source_outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_mouse_scONT_v1",
    )
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()
    render(args.stem.resolve(), args.dpi)


if __name__ == "__main__":
    main()
