#!/usr/bin/env python3
"""Render the six-panel scTHREAD NAR Figure 3 v6.

Figure contract
---------------
Core conclusion
    Matched UMAPs localize two contrasting PTPRC isoforms, while a separate
    usage-fraction panel provides the quantitative cell-type DTU evidence.
Evidence chain
    a, query and analysis status;
    b, shared cell-type UMAP reference;
    c-d, matched isoform-expression UMAPs;
    e, all-23-isoform-denominator usage fractions and DTU statistics;
    f, junction support plus a live selected-row detail and export/API access.
Guardrails
    Expression UMAPs localize the contrast but are not themselves DTU tests.
    Their native color scales differ (0-5 versus 0-3 molecules), so spatial
    distribution rather than absolute color intensity is compared. UMAP is
    exploratory and does not imply a trajectory.

The renderer is offline and deterministic. It reuses validated v4/v5 source
loaders and portal crops without altering the underlying data.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import nar_style as S
import render_nar_fig3_v3 as V3
import render_nar_fig3_v4 as V4
import render_nar_fig3_v5 as V5


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
JUNCTION_SOURCE = (
    FIGURES / "website_walkthrough" / "ptprc_views_v6_junction_vector"
)
JUNCTION_TABLE = JUNCTION_SOURCE / "junction_first100_live.tsv"
JUNCTION_METADATA = JUNCTION_SOURCE / "ptprc_all_views_metadata.json"
JUNCTION_REFERENCE = JUNCTION_SOURCE / "01_junctions.png"

WIDTH_MM = 183.0
HEIGHT_MM = 180.0

INK = V3.INK
SLATE = V3.SLATE
TEAL = V3.TEAL
BLUE = V3.BLUE
WHITE = V3.WHITE
PALE_BLUE = V3.PALE_BLUE
PURPLE = "#8B6F9E"
GOLD = V3.GOLD

GENE_START = 198_638_457
GENE_END = 198_757_476


def _load_junction_table() -> pd.DataFrame:
    frame = pd.read_csv(JUNCTION_TABLE, sep="\t")
    expected = {
        "junction",
        "span",
        "molecules",
        "reads",
        "runs",
        "studies",
    }
    if set(frame.columns) != expected:
        raise ValueError(
            f"Unexpected junction columns: {sorted(frame.columns)}"
        )
    parsed = frame["junction"].str.extract(
        r"^chr(?P<chrom>[^:]+):(?P<strand>[+-]):"
        r"(?P<start>[0-9]+)-(?P<end>[0-9]+)$"
    )
    if parsed.isna().any(axis=None):
        raise ValueError("Could not parse every exported junction identifier")
    frame = frame.join(parsed)
    for column in ("start", "end"):
        frame[column] = frame[column].astype(int)
    for column in ("span", "molecules", "reads", "runs", "studies"):
        frame[f"{column}_n"] = (
            frame[column]
            .astype(str)
            .str.replace(",", "", regex=False)
            .str.replace(" bp", "", regex=False)
            .astype(int)
        )
    if len(frame) != 100 or frame["junction"].nunique() != 100:
        raise ValueError("Expected 100 unique live portal junction rows")
    if not np.array_equal(
        frame["span_n"].to_numpy(),
        (frame["end"] - frame["start"]).to_numpy(),
    ):
        raise ValueError("Exported junction spans disagree with coordinates")
    if (
        (frame["end"] <= GENE_START).any()
        or (frame["start"] >= GENE_END).any()
        or set(frame["chrom"]) != {"1"}
        or not set(frame["strand"]).issubset({"+", "-"})
    ):
        raise ValueError("An exported junction does not overlap the PTPRC view")
    return frame


def _select_display_junction(junctions: pd.DataFrame) -> pd.Series:
    """Choose a strong, visibly long arc solely as an interface example."""
    candidates = junctions.loc[junctions["span_n"] >= 10_000].sort_values(
        ["molecules_n", "span_n", "junction"],
        ascending=[False, False, True],
    )
    if candidates.empty:
        raise ValueError("No junction is long enough for the display example")
    return candidates.iloc[0]


def _draw_cell_type_umap(fig: plt.Figure) -> None:
    V4._panel_heading(
        fig,
        "b",
        "Shared cell-type UMAP",
        x=0.055,
        y=0.705,
        title_size=6.35,
    )
    fig.text(
        0.055,
        0.674,
        "Harmonized cell-type annotation",
        fontsize=4.55,
        fontweight="bold",
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax = V4._image_axes(
        fig,
        [0.055, 0.445, 0.285, 0.215],
        V5.PRIMARY_CELL_MAP,
        x_slice=(0, 1140),
        y_slice=(290, 1150),
    )
    V5._label_relevant_cell_types(ax)
    fig.text(
        0.055,
        0.430,
        "71,913 sampled cells · shared coordinates",
        fontsize=4.2,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.412,
        "Median-coordinate labels · no trajectory inference",
        fontsize=4.15,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_isoform_umap(
    fig: plt.Figure,
    inputs: V3.Inputs,
    *,
    letter: str,
    x: float,
    image: Path,
    isoform: str,
    molecules: str,
    scale: str,
    usage_parts: tuple[tuple[str, str], ...],
    color: str,
) -> None:
    V4._panel_heading(
        fig,
        letter,
        f"{isoform} localization",
        x=x,
        y=0.705,
        title_size=5.95,
    )
    fig.text(
        x,
        0.674,
        "Isoform expression",
        fontsize=4.45,
        fontweight="bold",
        color=color,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    V4._image_axes(
        fig,
        [x, 0.445, 0.285, 0.215],
        image,
        x_slice=(1140, 2280),
        y_slice=(290, 1150),
    )
    fig.text(
        x,
        0.430,
        f"{molecules} molecules · native {scale} long-read molecules",
        fontsize=4.1,
        color=color,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    usage_text = " · ".join(
        f"{label} "
        f"{V5._usage_fraction(inputs, isoform, cell_type):.2f}"
        for label, cell_type in usage_parts
    )
    fig.text(
        x,
        0.412,
        f"Usage: {usage_text}",
        fontsize=4.1,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_umap_guardrail(fig: plt.Figure) -> None:
    fig.text(
        0.375,
        0.392,
        "Expression UMAPs localize the contrast; panel e quantifies DTU. "
        "Matched coordinates; compare spatial distribution, not native color intensity.",
        fontsize=4.15,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_dtu_quantification(
    fig: plt.Figure,
    inputs: V3.Inputs,
    *,
    replication_text: str | None = None,
) -> None:
    V4._panel_heading(
        fig,
        "e",
        (
            "Usage fractions and study-stratified DTU"
            if replication_text
            else "Usage fractions confirm cell-type DTU"
        ),
        x=0.055,
        y=0.355,
        title_size=6.25,
    )
    matrix = (
        inputs.usage.pivot(
            index="ct",
            columns="transcript_id",
            values="frac",
        )
        .reindex(V3.CELL_TYPE_ORDER)
        [[V5.PRIMARY_ISOFORM, V5.SECONDARY_ISOFORM]]
    )
    totals = (
        inputs.usage.groupby("ct")["count"]
        .sum()
        .reindex(V3.CELL_TYPE_ORDER)
        .round()
        .astype(int)
    )
    fig.text(
        0.055,
        0.324,
        f"CELL-TYPE DTU · q = {float(inputs.diu['qval']):.5f} "
        f"· effect = {float(inputs.diu['effect']):.3f} · 23 isoforms",
        fontsize=4.25,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.305,
        (
            replication_text
            if replication_text
            else (
                f"Fractions use the all-23-isoform denominator · "
                f"Σ = {int(totals.sum()):,} molecules"
            )
        ),
        fontsize=4.2,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    ax = fig.add_axes([0.135, 0.108, 0.145, 0.175])
    mesh = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1),
        np.arange(matrix.shape[0] + 1),
        matrix.to_numpy(),
        cmap="viridis_r",
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        linewidth=0.36,
        edgecolors=WHITE,
        rasterized=False,
    )
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    ax.set_xticks([0.5, 1.5])
    ax.set_xticklabels(
        ["…367364", "…697630"],
        fontsize=4.45,
        fontweight="bold",
    )
    ax.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        length=0,
        pad=1.5,
    )
    for tick, color in zip(ax.get_xticklabels(), [TEAL, BLUE]):
        tick.set_color(color)
    ax.set_yticks(np.arange(len(V3.CELL_TYPE_ORDER)) + 0.5)
    ax.set_yticklabels(V3.CELL_TYPE_ORDER, fontsize=4.35)
    ax.tick_params(axis="y", length=0, pad=1.8)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = float(matrix.iat[row, col])
            ax.text(
                col + 0.5,
                row + 0.5,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=4.1,
                fontweight="bold" if value >= 0.35 else "normal",
                color=WHITE if value >= 0.32 else INK,
            )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#BBC3C8")
        spine.set_linewidth(0.5)

    bar = fig.add_axes([0.315, 0.108, 0.075, 0.175])
    y = np.arange(len(V3.CELL_TYPE_ORDER))
    bar.barh(
        y,
        totals.to_numpy(),
        left=1.0,
        color=PALE_BLUE,
        edgecolor=BLUE,
        linewidth=0.42,
        height=0.58,
    )
    bar.set_xscale("log")
    bar.set_xlim(100, 100_000)
    bar.set_ylim(len(V3.CELL_TYPE_ORDER) - 0.5, -0.5)
    bar.set_xticks([])
    bar.set_yticks([])
    for spine in bar.spines.values():
        spine.set_visible(False)
    bar.text(
        0.0,
        1.08,
        "All 23 molecules",
        transform=bar.transAxes,
        fontsize=4.1,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
    )
    for row, value in enumerate(totals):
        bar.text(
            1.04,
            row,
            f"{value:,}",
            transform=bar.get_yaxis_transform(),
            fontsize=3.9,
            color=INK,
            ha="left",
            va="center",
            clip_on=False,
        )

    cax = fig.add_axes([0.135, 0.080, 0.145, 0.006])
    edges = np.linspace(0.0, 0.65, 66)
    values = (edges[:-1] + edges[1:]) / 2.0
    cax.pcolormesh(
        edges,
        [0.0, 1.0],
        values[np.newaxis, :],
        cmap=mesh.cmap,
        vmin=0.0,
        vmax=0.65,
        shading="flat",
        rasterized=False,
    )
    cax.set_xlim(0.0, 0.65)
    cax.set_ylim(0.0, 1.0)
    cax.set_yticks([])
    cax.set_xticks([0.0, 0.3, 0.6])
    cax.tick_params(axis="x", labelsize=3.9, length=1.4, width=0.4, pad=0.5)
    for spine in cax.spines.values():
        spine.set_linewidth(0.4)
        spine.set_color(INK)
    fig.text(
        0.2075,
        0.061,
        "Usage fraction",
        fontsize=4.05,
        color=SLATE,
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.315,
        0.061,
        (
            "Fractions use all 23 isoforms; q/effect summarize the gene-level test"
            if replication_text
            else "Fractions are descriptive; q/effect summarize the gene-level test"
        ),
        fontsize=3.75,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_access_card(fig: plt.Figure, selected: pd.Series) -> None:
    ax = fig.add_axes([0.825, 0.105, 0.160, 0.190])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    V4._rounded(
        ax,
        (0.0, 0.0),
        1.0,
        1.0,
        facecolor=INK,
        edgecolor=INK,
        linewidth=0,
        radius=0.035,
    )
    ax.text(
        0.07,
        0.915,
        "SELECTED ARC · DISPLAY EXAMPLE",
        fontsize=3.1,
        fontweight="bold",
        color="#E4B15B",
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.07,
        0.825,
        f"chr{selected['chrom']}:{selected['strand']}",
        fontsize=3.7,
        fontweight="bold",
        color=WHITE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.07,
        0.755,
        f"{int(selected['start']):,}–{int(selected['end']):,}",
        fontsize=3.55,
        color=WHITE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.plot([0.055, 0.945], [0.685, 0.685], color="#465063", linewidth=0.5)
    rows = [
        (0.07, 0.585, "MOLECULES", int(selected["molecules_n"])),
        (0.53, 0.585, "READS", int(selected["reads_n"])),
        (0.07, 0.395, "RUNS", int(selected["runs_n"])),
        (0.53, 0.395, "STUDIES", int(selected["studies_n"])),
    ]
    for x, y, label, value in rows:
        ax.text(
            x,
            y,
            label,
            fontsize=2.85,
            fontweight="bold",
            color="#AAB2C0",
            ha="left",
            va="center",
            fontfamily=S._FAM,
        )
        ax.text(
            x,
            y - 0.072,
            f"{value:,}",
            fontsize=3.8,
            fontweight="bold",
            color=WHITE,
            ha="left",
            va="center",
            fontfamily=S._FAM,
        )
    ax.plot([0.055, 0.945], [0.255, 0.255], color="#465063", linewidth=0.5)
    V4._rounded(
        ax,
        (0.055, 0.105),
        0.42,
        0.115,
        facecolor="#2D5D58",
        edgecolor="#6EA49D",
        linewidth=0.5,
        radius=0.025,
    )
    ax.text(
        0.265,
        0.162,
        "EXPORT CSV",
        fontsize=3.55,
        fontweight="bold",
        color=WHITE,
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.53,
        0.184,
        "GET /api/gene/PTPRC/",
        fontsize=2.8,
        fontweight="bold",
        color="#B7D8D3",
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.53,
        0.132,
        "overview?species=human",
        fontsize=2.75,
        color="#B7D8D3",
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax.text(
        0.50,
        0.050,
        "row: live table · overview: frozen snapshot",
        fontsize=2.65,
        color="#8FD1C6",
        ha="center",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_vector_junction_track(
    fig: plt.Figure,
    junctions: pd.DataFrame,
    selected: pd.Series,
) -> None:
    fig.text(
        0.525,
        0.319,
        "LIVE TABLE · 100 rows shown / 317 filtered",
        fontsize=4.15,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.525,
        0.300,
        "+ strand",
        fontsize=3.95,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.575,
        0.300,
        "− strand",
        fontsize=3.95,
        fontweight="bold",
        color=PURPLE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.620,
        0.300,
        "selected example",
        fontsize=3.65,
        fontweight="bold",
        color=GOLD,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.681,
        0.300,
        "width / opacity scale with log10(molecules)",
        fontsize=3.45,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    ax = fig.add_axes([0.525, 0.123, 0.285, 0.158])
    molecule_log = np.log10(junctions["molecules_n"].to_numpy(dtype=float))
    support = (molecule_log - molecule_log.min()) / (
        molecule_log.max() - molecule_log.min()
    )
    spans = junctions["span_n"].to_numpy(dtype=float)
    span_scale = np.sqrt(spans / spans.max())
    draw_order = np.argsort(support)
    t = np.linspace(0.0, 1.0, 96)
    for row_index in draw_order:
        row = junctions.iloc[int(row_index)]
        x = float(row["start"]) + float(row["span_n"]) * t
        height = 0.10 + 0.76 * float(span_scale[int(row_index)])
        y = height * 4.0 * t * (1.0 - t)
        ax.plot(
            x,
            y,
            color=TEAL if row["strand"] == "+" else PURPLE,
            linewidth=0.28 + 1.32 * float(support[int(row_index)]),
            alpha=0.18 + 0.68 * float(support[int(row_index)]),
            solid_capstyle="round",
        )

    selected_support = float(
        support[
            junctions.index.get_loc(
                junctions.index[
                    junctions["junction"] == selected["junction"]
                ][0]
            )
        ]
    )
    selected_span_scale = float(
        np.sqrt(float(selected["span_n"]) / spans.max())
    )
    selected_x = (
        float(selected["start"]) + float(selected["span_n"]) * t
    )
    selected_y = (
        0.10 + 0.76 * selected_span_scale
    ) * 4.0 * t * (1.0 - t)
    ax.plot(
        selected_x,
        selected_y,
        color=GOLD,
        linewidth=1.55 + 0.55 * selected_support,
        alpha=0.98,
        solid_capstyle="round",
        zorder=30,
    )
    ax.scatter(
        [float(selected["start"]), float(selected["end"])],
        [0.0, 0.0],
        s=5.0,
        color=GOLD,
        edgecolors=INK,
        linewidths=0.25,
        zorder=31,
        clip_on=True,
    )

    ax.axhline(0.0, color=INK, linewidth=0.62, zorder=20)
    ax.set_xlim(GENE_START, GENE_END)
    ax.set_ylim(-0.03, 0.95)
    ticks = np.array([GENE_START, (GENE_START + GENE_END) // 2, GENE_END])
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{value / 1_000_000:.3f}" for value in ticks],
        fontsize=4.0,
        color=SLATE,
    )
    ax.get_xticklabels()[0].set_ha("left")
    ax.get_xticklabels()[-1].set_ha("right")
    ax.set_xlabel(
        "chr1 genomic coordinate (Mb)",
        fontsize=4.15,
        color=SLATE,
        labelpad=2.0,
    )
    ax.set_yticks([])
    ax.tick_params(
        axis="x",
        length=2.0,
        width=0.45,
        colors=SLATE,
        pad=1.5,
    )
    for name, spine in ax.spines.items():
        spine.set_visible(name == "bottom")
        spine.set_color(INK)
        spine.set_linewidth(0.55)


def _draw_junction_access(
    fig: plt.Figure,
    junctions: pd.DataFrame,
) -> None:
    selected = _select_display_junction(junctions)
    V4._panel_heading(
        fig,
        "f",
        "Junction evidence remains inspectable and exportable",
        x=0.525,
        y=0.355,
        title_size=6.05,
    )
    _draw_vector_junction_track(fig, junctions, selected)
    fig.text(
        0.525,
        0.080,
        "Gold arc is a display example; selecting any arc exposes "
        "molecules, reads, runs and studies",
        fontsize=3.85,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    _draw_access_card(fig, selected)


def render(
    inputs: V3.Inputs,
    junctions: pd.DataFrame,
    stem: Path,
    dpi: int,
) -> list[Path]:
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
        fig = plt.figure(
            figsize=(WIDTH_MM * S.MM, HEIGHT_MM * S.MM),
            facecolor="none",
        )
        V4._draw_query_panel(fig, inputs)
        _draw_cell_type_umap(fig)
        _draw_isoform_umap(
            fig,
            inputs,
            letter="c",
            x=0.375,
            image=V5.PRIMARY_CELL_MAP,
            isoform=V5.PRIMARY_ISOFORM,
            molecules="57,185",
            scale="0–5",
            usage_parts=(
                ("Monocyte", "Monocyte"),
                ("Plasma", "Plasma cell"),
                ("Progenitor", "Progenitor"),
            ),
            color=TEAL,
        )
        _draw_isoform_umap(
            fig,
            inputs,
            letter="d",
            x=0.695,
            image=V5.SECONDARY_CELL_MAP,
            isoform=V5.SECONDARY_ISOFORM,
            molecules="25,697",
            scale="0–3",
            usage_parts=(
                ("Progenitor", "Progenitor"),
                ("B cell", "B cell"),
                ("Monocyte", "Monocyte"),
            ),
            color=BLUE,
        )
        _draw_umap_guardrail(fig)
        _draw_dtu_quantification(fig, inputs)
        _draw_junction_access(fig, junctions)

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
        default=FIGURES / "NAR_Fig3_v6",
        help="Output path without an extension.",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate all frozen inputs without rendering.",
    )
    args = parser.parse_args()

    for source in (
        V5.PRIMARY_CELL_MAP,
        V5.SECONDARY_CELL_MAP,
        V5.CENTROID_TABLE,
        JUNCTION_TABLE,
        JUNCTION_METADATA,
        JUNCTION_REFERENCE,
    ):
        if not source.is_file():
            raise FileNotFoundError(source)
    inputs = V3.load_inputs()
    junctions = _load_junction_table()
    print(
        "INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "cell_types": inputs.usage["ct"].nunique(),
            "isoforms": inputs.usage["transcript_id"].nunique(),
            "matrix_molecules": int(round(inputs.usage["count"].sum())),
            "panels": 6,
            "portal_images": 4,
            "vector_junctions": len(junctions),
            "canvas_mm": [WIDTH_MM, HEIGHT_MM],
        },
    )
    if args.validate_only:
        return
    outputs = render(inputs, junctions, args.stem.resolve(), args.dpi)
    for output in outputs:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
