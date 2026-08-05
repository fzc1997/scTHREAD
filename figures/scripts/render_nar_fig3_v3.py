#!/usr/bin/env python3
"""Render the publication-grade scTHREAD NAR Figure 3 v3.

The figure combines one live portal crop with vector summaries rebuilt from the
frozen PTPRC API snapshot and the underlying DIU, APA, ASE and isoform-usage
tables.  The renderer never queries the network; live parity is checked by the
companion QA script.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

import nar_style as S


PROJECT = Path(__file__).resolve().parents[2]
SC_ROOT = PROJECT.parent
FIGURES = PROJECT / "figures"
ASSETS = FIGURES / "assets"
TABLES = PROJECT / "tables"
FIGDATA = SC_ROOT / "results" / "paper1" / "figdata"
F2DATA = SC_ROOT / "results" / "paper1" / "f2_grammar" / "figdata"

GENE = "PTPRC"
GID = "ENSG00000081237"
WIDTH_MM = 183.0
HEIGHT_MM = 165.0

INK = "#1C2130"
SLATE = "#55606F"
TEAL = "#2F6E6B"
BLUE = "#3D6F9B"
GOLD = "#B37A20"
CORAL = "#C1503A"
GREY = "#8A8F98"
PALE = "#F2F5F5"
PALE_BLUE = "#EDF3F7"
WHITE = "#FFFFFF"
GRID = "#D8DEE2"

CELL_TYPE_ORDER = [
    "Progenitor",
    "B cell",
    "T cell",
    "NK",
    "Monocyte",
    "Dendritic cell",
    "Plasma cell",
    "Erythroid",
]


@dataclass(frozen=True)
class Inputs:
    diu: pd.Series
    apa: pd.Series
    ase: pd.Series
    usage: pd.DataFrame
    api: dict
    capture: dict
    portal_image: Path


def _single_row(frame: pd.DataFrame, mask: pd.Series, label: str) -> pd.Series:
    rows = frame.loc[mask]
    if len(rows) != 1:
        raise ValueError(f"Expected one {label} row for {GENE}; found {len(rows)}")
    return rows.iloc[0]


def _as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def _close(left: object, right: object, label: str, atol: float = 1e-12) -> None:
    if not np.isclose(float(left), float(right), rtol=0.0, atol=atol):
        raise ValueError(f"{label} differs: {left!r} versus {right!r}")


def load_inputs() -> Inputs:
    diu_table = pd.read_csv(F2DATA / "diu_celltype.tsv", sep="\t")
    apa_table = pd.read_csv(F2DATA / "apa_celltype.tsv", sep="\t")
    ase_table = pd.read_csv(F2DATA / "ase_interaction.tsv", sep="\t")
    usage = pd.read_csv(FIGDATA / "ptprc_isoform_usage.tsv", sep="\t")
    api = json.loads((TABLES / "PTPRC_api_overview.json").read_text())
    capture = json.loads((ASSETS / "ptprc_live_capture_v3.json").read_text())

    diu = _single_row(diu_table, diu_table["gene"].eq(GID), "DIU")
    apa = _single_row(apa_table, apa_table["gene"].eq(GID), "APA")
    ase = _single_row(
        ase_table,
        ase_table["gene"].astype(str).str.upper().eq(GENE),
        "ASE",
    )
    inputs = Inputs(
        diu=diu,
        apa=apa,
        ase=ase,
        usage=usage,
        api=api,
        capture=capture,
        portal_image=ASSETS / "ptprc_gene_card_live_v3.png",
    )
    validate_inputs(inputs)
    return inputs


def validate_inputs(inputs: Inputs) -> None:
    if not inputs.portal_image.is_file():
        raise FileNotFoundError(inputs.portal_image)
    gene = inputs.api["gene"]
    if gene["gid"] != GID or gene["gname"] != GENE:
        raise ValueError("Frozen API snapshot is not the expected PTPRC record")

    for name, row in (("diu", inputs.diu), ("apa", inputs.apa), ("ase", inputs.ase)):
        api_row = inputs.api["analyses"][name]
        for field in ("pval", "effect", "qval"):
            _close(row[field], api_row[field], f"{name}.{field}")
        if _as_bool(row["sig"]) != bool(api_row["sig"]):
            raise ValueError(f"{name}.sig differs between source table and API snapshot")
    if int(inputs.diu["n_iso"]) != int(inputs.api["analyses"]["diu"]["n_iso"]):
        raise ValueError("DIU isoform count differs between table and API snapshot")
    if int(inputs.apa["n_pas"]) != int(inputs.api["analyses"]["apa"]["n_pas"]):
        raise ValueError("APA site count differs between table and API snapshot")

    usage = inputs.usage
    expected_columns = {"ct", "transcript_id", "count", "frac"}
    if not expected_columns.issubset(usage.columns):
        raise ValueError(f"PTPRC usage table lacks {expected_columns - set(usage.columns)}")
    if usage["transcript_id"].nunique() != 23 or usage["ct"].nunique() != 8:
        raise ValueError("PTPRC usage table no longer contains 23 isoforms × 8 cell types")
    observed_ct = set(usage["ct"])
    if observed_ct != set(CELL_TYPE_ORDER):
        raise ValueError(f"Unexpected PTPRC cell types: {sorted(observed_ct)}")
    totals = usage.groupby("ct")["count"].transform("sum")
    max_error = np.abs(usage["frac"] - usage["count"] / totals).max()
    if max_error > 1e-12:
        raise ValueError(f"Usage fractions do not match molecule counts (max error {max_error})")
    row_sums = usage.groupby("ct")["frac"].sum()
    if not np.allclose(row_sums.to_numpy(), 1.0, atol=1e-12, rtol=0.0):
        raise ValueError("PTPRC usage fractions are not normalized over all 23 isoforms")
    if int(round(usage["count"].sum())) != 98_445:
        raise ValueError("Eight-cell-type PTPRC molecule total changed from 98,445")

    location = inputs.capture.get("browser_location", "")
    evidence = inputs.capture.get("evidence_text", "")
    required_capture_strings = (
        GENE,
        GID,
        "317 junctions",
        "8,994 junctions",
        "2,341,574 molecules",
    )
    capture_text = f"{location}\n{evidence}"
    missing = [item for item in required_capture_strings if item not in capture_text]
    if missing:
        raise ValueError(f"Live portal capture metadata lacks: {missing}")


def _panel_heading(
    fig: plt.Figure,
    letter: str,
    title: str,
    *,
    x: float,
    y: float,
    title_size: float = 7.2,
) -> None:
    fig.text(
        x - 0.030,
        y,
        letter,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
        fontfamily=S._FAM,
        color=INK,
    )
    fig.text(
        x,
        y,
        title,
        ha="left",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        fontfamily=S._FAM,
        color=INK,
    )


def _rounded(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    *,
    facecolor: str,
    edgecolor: str = "none",
    linewidth: float = 0.6,
    radius: float = 0.025,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def _analysis_card(
    ax: plt.Axes,
    x: float,
    *,
    acronym: str,
    label: str,
    qval: float,
    effect: float,
    detail: str,
    significant: bool,
    color: str,
) -> None:
    width = 0.143
    _rounded(
        ax,
        (x, 0.05),
        width,
        0.86,
        facecolor=WHITE,
        edgecolor=GRID,
        linewidth=0.55,
        radius=0.025,
    )
    ax.add_patch(
        Rectangle(
            (x, 0.83),
            width,
            0.08,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            clip_on=False,
        )
    )
    ax.text(
        x + 0.012,
        0.74,
        acronym,
        transform=ax.transAxes,
        fontsize=6.5,
        fontweight="bold",
        color=color,
        ha="left",
        va="center",
    )
    ax.text(
        x + width - 0.012,
        0.74,
        "SIGNIFICANT" if significant else "NOT SIGNIFICANT",
        transform=ax.transAxes,
        fontsize=4.3,
        fontweight="bold",
        color=color if significant else GREY,
        ha="right",
        va="center",
    )
    ax.text(
        x + 0.012,
        0.57,
        label,
        transform=ax.transAxes,
        fontsize=5.0,
        color=SLATE,
        ha="left",
        va="center",
    )
    ax.text(
        x + 0.012,
        0.38,
        f"q = {qval:.3g}",
        transform=ax.transAxes,
        fontsize=6.1,
        fontweight="bold",
        color=INK,
        ha="left",
        va="center",
    )
    ax.text(
        x + 0.012,
        0.24,
        f"effect = {effect:.3f}",
        transform=ax.transAxes,
        fontsize=5.1,
        color=INK,
        ha="left",
        va="center",
    )
    ax.text(
        x + 0.012,
        0.10,
        detail,
        transform=ax.transAxes,
        fontsize=4.8,
        color=SLATE,
        ha="left",
        va="center",
    )


def _draw_portal_panel(fig: plt.Figure, inputs: Inputs) -> None:
    _panel_heading(
        fig,
        "a",
        "A live PTPRC query exposes filtered and database-wide evidence in one gene card",
        x=0.055,
        y=0.975,
    )
    ax = fig.add_axes([0.055, 0.680, 0.930, 0.262])
    image = mpimg.imread(inputs.portal_image)
    ax.imshow(image, interpolation="lanczos", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#AEB8BF")
        spine.set_linewidth(0.6)
    fig.text(
        0.055,
        0.660,
        f"{GENE} · {GID} · 317 junctions = current browser subset after the "
        "≥10-molecule filter; "
        "8,994 junctions = database-wide PTPRC coverage.",
        fontsize=5.2,
        color=SLATE,
        ha="left",
        va="center",
    )


def _draw_evidence_panel(fig: plt.Figure, inputs: Inputs) -> None:
    _panel_heading(
        fig,
        "b",
        "Significant isoform and poly(A) usage are supported by broad multi-layer coverage",
        x=0.055,
        y=0.630,
    )
    ax = fig.add_axes([0.055, 0.455, 0.930, 0.145])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _analysis_card(
        ax,
        0.000,
        acronym="DIU",
        label="Isoform usage",
        qval=float(inputs.diu["qval"]),
        effect=float(inputs.diu["effect"]),
        detail=f"{int(inputs.diu['n_iso'])} isoforms",
        significant=_as_bool(inputs.diu["sig"]),
        color=TEAL,
    )
    _analysis_card(
        ax,
        0.156,
        acronym="APA",
        label="Poly(A) usage",
        qval=float(inputs.apa["qval"]),
        effect=float(inputs.apa["effect"]),
        detail=f"{int(inputs.apa['n_pas'])} tested sites",
        significant=_as_bool(inputs.apa["sig"]),
        color=GOLD,
    )
    _analysis_card(
        ax,
        0.312,
        acronym="ASE",
        label="Cell-type interaction",
        qval=float(inputs.ase["qval"]),
        effect=float(inputs.ase["effect"]),
        detail=f"{inputs.ase['hi_ct']} vs {inputs.ase['lo_ct']}",
        significant=_as_bool(inputs.ase["sig"]),
        color=GREY,
    )

    coverage = inputs.api["coverage"]
    x0, x1 = 0.500, 0.995
    _rounded(
        ax,
        (x0, 0.05),
        x1 - x0,
        0.86,
        facecolor=PALE,
        edgecolor=GRID,
        linewidth=0.55,
        radius=0.022,
    )
    columns = [
        ("LAYER", x0 + 0.020, "left"),
        ("FEATURES", x0 + 0.220, "right"),
        ("MOLECULES", x0 + 0.350, "right"),
        ("STUDIES / CELL TYPES", x0 + 0.478, "right"),
    ]
    for text, x, align in columns:
        ax.text(
            x,
            0.78,
            text,
            transform=ax.transAxes,
            fontsize=4.4,
            fontweight="bold",
            color=SLATE,
            ha=align,
            va="center",
        )
    rows = [
        (
            "Isoforms",
            coverage["isoforms"]["features"],
            coverage["isoforms"]["molecules"],
            coverage["isoforms"]["studies"],
            coverage["isoforms"]["cell_types"],
        ),
        (
            "Poly(A) sites",
            coverage["pas"]["features"],
            coverage["pas"]["molecules"],
            coverage["pas"]["studies"],
            coverage["pas"]["cell_types"],
        ),
        (
            "Junctions",
            coverage["junctions"]["features"],
            coverage["junctions"]["molecules"],
            coverage["junctions"]["studies"],
            "—",
        ),
    ]
    row_y = [0.58, 0.37, 0.16]
    for idx, ((label, features, molecules, studies, cell_types), y) in enumerate(
        zip(rows, row_y)
    ):
        if idx:
            ax.plot(
                [x0 + 0.015, x1 - 0.010],
                [y + 0.105, y + 0.105],
                transform=ax.transAxes,
                color=GRID,
                linewidth=0.45,
            )
        values = [
            (label, x0 + 0.020, "left", True),
            (f"{int(features):,}", x0 + 0.220, "right", False),
            (f"{int(molecules):,}", x0 + 0.350, "right", False),
            (
                f"{int(studies):,} / "
                + (f"{int(cell_types):,}" if cell_types != "—" else "—"),
                x0 + 0.478,
                "right",
                False,
            ),
        ]
        for value, x, align, bold in values:
            ax.text(
                x,
                y,
                value,
                transform=ax.transAxes,
                fontsize=5.35,
                fontweight="bold" if bold else "normal",
                color=INK,
                ha=align,
                va="center",
            )


def _usage_matrix(usage: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    top = (
        usage.groupby("transcript_id", sort=False)["count"]
        .sum()
        .nlargest(5)
        .index.tolist()
    )
    subset = usage.loc[usage["transcript_id"].isin(top)].copy()
    matrix = (
        subset.pivot_table(
            index="ct",
            columns="transcript_id",
            values="frac",
            aggfunc="sum",
        )
        .reindex(CELL_TYPE_ORDER)
        .fillna(0.0)
    )
    matrix = matrix[matrix.loc["Monocyte"].sort_values(ascending=False).index]
    totals = usage.groupby("ct")["count"].sum().reindex(CELL_TYPE_ORDER)
    top_five = subset.groupby("ct")["count"].sum().reindex(CELL_TYPE_ORDER)
    return matrix, totals, top_five / totals


def _draw_usage_panel(fig: plt.Figure, inputs: Inputs) -> None:
    _panel_heading(
        fig,
        "c",
        "PTPRC isoform usage contrasts monocytes with progenitor/B/NK cells",
        x=0.055,
        y=0.425,
        title_size=6.8,
    )
    matrix, totals, coverage = _usage_matrix(inputs.usage)

    ax = fig.add_axes([0.120, 0.105, 0.315, 0.255])
    image = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1),
        np.arange(matrix.shape[0] + 1),
        matrix.to_numpy(),
        cmap="viridis_r",
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        linewidth=0.35,
        edgecolors=WHITE,
        rasterized=False,
    )
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    short_labels = ["…" + transcript[-6:] for transcript in matrix.columns]
    ax.set_xticks(np.arange(len(short_labels)) + 0.5)
    ax.set_xticklabels(short_labels, rotation=27, ha="left", rotation_mode="anchor")
    ax.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        length=0,
        pad=1.5,
        labelsize=5.2,
    )
    for idx, tick in enumerate(ax.get_xticklabels()):
        if idx == 0:
            tick.set_color(TEAL)
            tick.set_fontweight("bold")
        elif idx == 1:
            tick.set_color(BLUE)
            tick.set_fontweight("bold")
    ax.set_yticks(np.arange(len(CELL_TYPE_ORDER)) + 0.5)
    ax.set_yticklabels(CELL_TYPE_ORDER, fontsize=5.7)
    ax.tick_params(axis="y", length=0, pad=2.0)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = float(matrix.iat[row, col])
            ax.text(
                col + 0.5,
                row + 0.5,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=4.9,
                fontweight="bold" if value >= 0.35 else "normal",
                color=WHITE if value >= 0.32 else INK,
            )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#BBC3C8")
        spine.set_linewidth(0.5)

    cax = fig.add_axes([0.120, 0.069, 0.315, 0.008])
    color_edges = np.linspace(0.0, 0.65, 66)
    color_values = (color_edges[:-1] + color_edges[1:]) / 2.0
    cax.pcolormesh(
        color_edges,
        [0.0, 1.0],
        color_values[np.newaxis, :],
        cmap="viridis_r",
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        rasterized=False,
    )
    cax.set_xlim(0.0, 0.65)
    cax.set_ylim(0.0, 1.0)
    cax.set_yticks([])
    cax.set_xticks([0.0, 0.3, 0.6])
    cax.tick_params(axis="x", labelsize=4.5, length=1.8, width=0.4, pad=0.8)
    for spine in cax.spines.values():
        spine.set_linewidth(0.4)
        spine.set_color(INK)
    fig.text(
        0.2775,
        0.050,
        "Usage fraction across all 23 isoforms",
        fontsize=4.8,
        color=SLATE,
        ha="center",
        va="center",
    )

    bar = fig.add_axes([0.468, 0.105, 0.135, 0.255])
    y = np.arange(len(CELL_TYPE_ORDER))
    bar.barh(
        y,
        totals.to_numpy(),
        left=1.0,
        color=PALE_BLUE,
        edgecolor=BLUE,
        linewidth=0.45,
        height=0.60,
    )
    bar.set_xscale("log")
    bar.set_xlim(100, 100_000)
    bar.set_ylim(len(CELL_TYPE_ORDER) - 0.5, -0.5)
    bar.set_yticks([])
    bar.set_xticks([])
    for spine in bar.spines.values():
        spine.set_visible(False)
    bar.text(
        0.0,
        1.075,
        "All-isoform molecules",
        transform=bar.transAxes,
        fontsize=5.2,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
    )
    bar.text(
        0.0,
        1.015,
        f"8-cell matrix · Σ = {int(totals.sum()):,} · log scale",
        transform=bar.transAxes,
        fontsize=4.5,
        color=SLATE,
        ha="left",
        va="bottom",
    )
    for row, value in enumerate(totals.astype(int)):
        bar.text(
            1.02,
            row,
            f"{value:,}",
            transform=bar.get_yaxis_transform(),
            fontsize=4.8,
            color=INK,
            ha="left",
            va="center",
            clip_on=False,
        )
    fig.text(
        0.468,
        0.050,
        f"Top five = {coverage.min() * 100:.0f}–{coverage.max() * 100:.0f}% per row "
        "(values not renormalized)",
        fontsize=4.65,
        color=SLATE,
        ha="left",
        va="center",
    )


def _draw_api_panel(fig: plt.Figure, inputs: Inputs) -> None:
    _panel_heading(
        fig,
        "d",
        "API fields match portal and source tables",
        x=0.680,
        y=0.425,
        title_size=6.9,
    )
    ax = fig.add_axes([0.680, 0.055, 0.305, 0.335])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    _rounded(
        ax,
        (0.0, 0.0),
        1.0,
        1.0,
        facecolor=INK,
        edgecolor=INK,
        linewidth=0.0,
        radius=0.035,
    )
    ax.text(
        0.045,
        0.935,
        "GET /api/gene/PTPRC/overview?species=human",
        fontsize=4.8,
        fontweight="bold",
        color="#B7D8D3",
        ha="left",
        va="center",
    )
    ax.plot([0.04, 0.96], [0.885, 0.885], color="#465063", linewidth=0.5)

    ax.text(0.045, 0.835, "ANALYSIS OBJECT", fontsize=4.15, color="#AAB2C0", va="center")
    ax.text(0.570, 0.835, "QVAL", fontsize=4.15, color="#AAB2C0", ha="right", va="center")
    ax.text(0.790, 0.835, "EFFECT", fontsize=4.15, color="#AAB2C0", ha="right", va="center")
    ax.text(0.945, 0.835, "SIG", fontsize=4.15, color="#AAB2C0", ha="right", va="center")
    analysis_rows = [
        ("analyses.diu", inputs.api["analyses"]["diu"]),
        ("analyses.apa", inputs.api["analyses"]["apa"]),
        ("analyses.ase", inputs.api["analyses"]["ase"]),
    ]
    for y, (name, values) in zip([0.755, 0.675, 0.595], analysis_rows):
        ax.text(0.045, y, name, fontsize=4.9, color=WHITE, ha="left", va="center")
        ax.text(
            0.570,
            y,
            f"{float(values['qval']):.6f}",
            fontsize=4.75,
            color=WHITE,
            ha="right",
            va="center",
        )
        ax.text(
            0.790,
            y,
            f"{float(values['effect']):.6f}",
            fontsize=4.75,
            color=WHITE,
            ha="right",
            va="center",
        )
        ax.text(
            0.945,
            y,
            str(bool(values["sig"])).lower(),
            fontsize=4.7,
            fontweight="bold",
            color="#70B7AC" if values["sig"] else "#C9CED6",
            ha="right",
            va="center",
        )

    ax.plot([0.04, 0.96], [0.535, 0.535], color="#465063", linewidth=0.5)
    ax.text(0.045, 0.485, "COVERAGE OBJECT", fontsize=4.15, color="#AAB2C0", va="center")
    ax.text(0.690, 0.485, "FEATURES", fontsize=4.15, color="#AAB2C0", ha="right", va="center")
    ax.text(0.945, 0.485, "MOLECULES", fontsize=4.15, color="#AAB2C0", ha="right", va="center")
    coverage_rows = [
        ("coverage.isoforms", inputs.api["coverage"]["isoforms"]),
        ("coverage.pas", inputs.api["coverage"]["pas"]),
        ("coverage.junctions", inputs.api["coverage"]["junctions"]),
    ]
    for y, (name, values) in zip([0.405, 0.325, 0.245], coverage_rows):
        ax.text(0.045, y, name, fontsize=4.9, color=WHITE, ha="left", va="center")
        ax.text(
            0.690,
            y,
            f"{int(values['features']):,}",
            fontsize=4.8,
            color=WHITE,
            ha="right",
            va="center",
        )
        ax.text(
            0.945,
            y,
            f"{int(values['molecules']):,}",
            fontsize=4.8,
            color=WHITE,
            ha="right",
            va="center",
        )

    _rounded(
        ax,
        (0.040, 0.065),
        0.920,
        0.105,
        facecolor="#263B43",
        edgecolor="#55767B",
        linewidth=0.5,
        radius=0.018,
    )
    ax.text(
        0.065,
        0.118,
        "NUMERIC PARITY VALIDATED",
        fontsize=4.55,
        fontweight="bold",
        color="#8FD1C6",
        ha="left",
        va="center",
    )
    ax.text(
        0.935,
        0.118,
        "tables = snapshot = live API",
        fontsize=4.45,
        color=WHITE,
        ha="right",
        va="center",
    )


def render(inputs: Inputs, stem: Path, dpi: int) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(
        {
            "font.family": S._FAM,
            "font.size": 6.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
            "axes.unicode_minus": False,
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
        }
    ):
        fig = plt.figure(figsize=(WIDTH_MM * S.MM, HEIGHT_MM * S.MM), facecolor="none")
        _draw_portal_panel(fig, inputs)
        _draw_evidence_panel(fig, inputs)
        _draw_usage_panel(fig, inputs)
        _draw_api_panel(fig, inputs)

        outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
        for output in outputs:
            fig.savefig(
                output,
                dpi=dpi,
                transparent=True,
                bbox_inches=None,
                pad_inches=0,
            )
        plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v3",
        help="Output path without an extension.",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all frozen inputs without rendering.",
    )
    args = parser.parse_args()

    inputs = load_inputs()
    print(
        "INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "cell_types": inputs.usage["ct"].nunique(),
            "isoforms": inputs.usage["transcript_id"].nunique(),
            "matrix_molecules": int(round(inputs.usage["count"].sum())),
        },
    )
    if args.validate_only:
        return
    outputs = render(inputs, args.stem.resolve(), args.dpi)
    for output in outputs:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
