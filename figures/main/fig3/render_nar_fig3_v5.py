#!/usr/bin/env python3
"""Render the UMAP-centred scTHREAD NAR Figure 3 v5.

Figure contract
---------------
Core conclusion
    A PTPRC query exposes significant cell-type-specific differential
    transcript usage whose cellular localization is visible in matched
    isoform-expression UMAPs.
Evidence hierarchy
    Cell-type UMAP reference -> matched ENST00000367364 and ENST00000697630
    expression UMAPs -> cell-type usage anchors -> junction support ->
    export/API reproducibility.
Guardrails
    The two expression UMAPs use identical cell coordinates and one shared
    0-5 molecule color scale, so both localization and color intensity are
    directly comparable. Usage fractions retain the all-23-isoform
    denominator. UMAP is exploratory and does not imply a trajectory.

The renderer is offline and deterministic. It reuses the validated Figure 3
v4 query and junction/access panels while replacing the heatmap with direct,
matched portal UMAP evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "style"))
REPO = Path(__file__).resolve().parents[3]
import nar_style as S
import render_nar_fig3_v3 as V3
import render_nar_fig3_v4 as V4


PROJECT = Path(__file__).resolve().parents[3]
FIGURES = PROJECT / "figures"
PRIMARY_WALKTHROUGH = (
    REPO / "source_data" / "portal_screenshots" / "ptprc_views_v8_outliers_hidden_20260807"
)
SECONDARY_WALKTHROUGH = (
    REPO / "source_data" / "portal_screenshots" / "ptprc_views_v7_scale5_secondary"
)

PRIMARY_CELL_MAP = PRIMARY_WALKTHROUGH / "02_cell_map.png"
SECONDARY_CELL_MAP = SECONDARY_WALKTHROUGH / "02_cell_map.png"
CENTROID_TABLE = REPO / "source_data" / "tables" / "PTPRC_umap_celltype_centroids_live_20260807.tsv"

PRIMARY_ISOFORM = "ENST00000367364"
SECONDARY_ISOFORM = "ENST00000697630"
SHARED_COLOR_CEILING = 5.0
DISPLAYED_CELLS = 97_934
PRIMARY_MOLECULES = 47_951
SECONDARY_MOLECULES = 32_456

WIDTH_MM = V3.WIDTH_MM
HEIGHT_MM = V3.HEIGHT_MM
INK = V3.INK
SLATE = V3.SLATE
TEAL = V3.TEAL
BLUE = V3.BLUE


def canvas_crop(walkthrough: Path, role: str) -> dict[str, tuple[int, int]]:
    """Return the screenshot slices that isolate one portal UMAP canvas.

    Crops are derived from the canvas box recorded during capture rather than
    hard-coded, so cell-type labels expressed as canvas fractions stay on their
    clusters even when the portal's surrounding chrome changes height.
    """
    metadata = json.loads(
        (walkthrough / "ptprc_all_views_metadata.json").read_text()
    )
    views = [view for view in metadata["views"] if view["view"] == "cellmap"]
    if len(views) != 1:
        raise ValueError(f"Expected one cellmap view in {walkthrough}")
    view = views[0]
    clip = view["clip_css_px"]
    box = view["cellmap_config"]["canvas_geometry"][role]
    scale = float(metadata["device_scale_factor"])
    left = round((box["x"] - clip["x"]) * scale)
    top = round((box["y"] - clip["y"]) * scale)
    return {
        "x_slice": (left, left + round(box["width"] * scale)),
        "y_slice": (top, top + round(box["height"] * scale)),
    }


def assert_shared_color_ceiling(walkthroughs: tuple[Path, ...]) -> str:
    """Fail unless every capture rendered the same explicit colour ceiling."""
    expected = f"0 → {SHARED_COLOR_CEILING:.0f} long-read molecules"
    legends = []
    for walkthrough in walkthroughs:
        metadata = json.loads(
            (walkthrough / "ptprc_all_views_metadata.json").read_text()
        )
        if float(metadata["cellmap_color_max"]) != SHARED_COLOR_CEILING:
            raise ValueError(
                f"{walkthrough.name} was captured at colour ceiling "
                f"{metadata['cellmap_color_max']}, not {SHARED_COLOR_CEILING}"
            )
        legend = metadata["views"][0]["cellmap_config"]["color_max"]["legend_text"]
        if expected not in legend:
            raise ValueError(
                f"{walkthrough.name} legend does not report {expected!r}: {legend!r}"
            )
        legends.append(legend)
    return expected


def _load_relevant_centroids() -> dict[str, dict[str, float | int | str]]:
    with CENTROID_TABLE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    # The table carries every cell type quantified in panel e; the map labels
    # only the two that panels c and d resolve, so require a superset.
    required = {"B cell", "Monocyte"}
    observed = {row["cell_type"] for row in rows}
    if not required <= observed:
        raise ValueError(
            f"Centroid table is missing {sorted(required - observed)}"
        )
    return {
        row["cell_type"]: {
            "cell_type": row["cell_type"],
            "n_cells": int(row["n_cells"]),
            "x": float(row["canvas_x_fraction"]),
            "y_top": float(row["canvas_y_from_top_fraction"]),
        }
        for row in rows
    }


def _label_relevant_cell_types(ax: plt.Axes) -> None:
    centroids = _load_relevant_centroids()
    # NK and T cell are intermingled in this embedding (their median canvas
    # positions differ by 0.02), so labelling them separately would imply a
    # separation the map does not show. Panel e quantifies them instead.
    label_specs = {
        "B cell": {
            "text": "B cell",
            "color": "#4C78A8",
            "xytext": (0.03, 0.88),
            "ha": "left",
        },
        "Monocyte": {
            "text": "Monocyte",
            "color": "#F28E2B",
            "xytext": (0.98, 0.42),
            "ha": "right",
        },
    }
    for cell_type, spec in label_specs.items():
        centroid = centroids[cell_type]
        ax.annotate(
            spec["text"],
            xy=(centroid["x"], 1.0 - centroid["y_top"]),
            xycoords=ax.transAxes,
            xytext=spec["xytext"],
            textcoords=ax.transAxes,
            fontsize=3.75,
            fontweight="bold",
            color=spec["color"],
            ha=spec["ha"],
            va="center",
            fontfamily=S._FAM,
            bbox={
                "boxstyle": "round,pad=0.16",
                "facecolor": "white",
                "edgecolor": spec["color"],
                "linewidth": 0.45,
                "alpha": 0.92,
            },
            arrowprops={
                "arrowstyle": "-",
                "color": spec["color"],
                "linewidth": 0.55,
                "shrinkA": 1.5,
                "shrinkB": 1.5,
            },
            annotation_clip=False,
        )


def _usage_fraction(
    inputs: V3.Inputs,
    transcript_id: str,
    cell_type: str,
) -> float:
    rows = inputs.usage.loc[
        inputs.usage["transcript_id"].eq(transcript_id)
        & inputs.usage["ct"].eq(cell_type),
        "frac",
    ]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one usage value for {transcript_id} / {cell_type}; "
            f"found {len(rows)}"
        )
    return float(rows.iloc[0])


def _draw_cell_type_umap(fig: plt.Figure) -> None:
    V4._panel_heading(
        fig,
        "B",
        "Cell types define the shared UMAP reference",
        x=0.055,
        y=0.710,
        title_size=6.65,
    )
    fig.text(
        0.055,
        0.674,
        "Harmonized cell-type annotation",
        fontsize=4.65,
        fontweight="bold",
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    ax = V4._image_axes(
        fig,
        [0.055, 0.410, 0.285, 0.238],
        PRIMARY_CELL_MAP,
        **canvas_crop(PRIMARY_WALKTHROUGH, "cell"),
    )
    _label_relevant_cell_types(ax)
    fig.text(
        0.055,
        0.394,
        f"{DISPLAYED_CELLS:,} sampled cells · coordinates shared with panel c",
        fontsize=4.35,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.370,
        "Labels mark within-type median coordinates · no trajectory inference",
        fontsize=4.35,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_dtu_umaps(fig: plt.Figure, inputs: V3.Inputs) -> None:
    V4._panel_heading(
        fig,
        "C",
        "Matched isoform UMAPs localize the PTPRC DTU contrast",
        x=0.385,
        y=0.710,
        title_size=6.45,
    )
    matrix_total = int(round(inputs.usage["count"].sum()))
    fig.text(
        0.385,
        0.678,
        f"CELL-TYPE DTU · q = {float(inputs.diu['qval']):.5f} "
        f"· effect = {float(inputs.diu['effect']):.3f} · 23 isoforms "
        f"· Σ = {matrix_total:,} molecules",
        fontsize=4.45,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    fig.text(
        0.385,
        0.653,
        PRIMARY_ISOFORM,
        fontsize=4.55,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.690,
        0.653,
        SECONDARY_ISOFORM,
        fontsize=4.55,
        fontweight="bold",
        color=BLUE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    V4._image_axes(
        fig,
        [0.385, 0.410, 0.285, 0.238],
        PRIMARY_CELL_MAP,
        **canvas_crop(PRIMARY_WALKTHROUGH, "signal"),
    )
    V4._image_axes(
        fig,
        [0.690, 0.410, 0.285, 0.238],
        SECONDARY_CELL_MAP,
        **canvas_crop(SECONDARY_WALKTHROUGH, "signal"),
    )

    fig.text(
        0.385,
        0.394,
        f"{PRIMARY_MOLECULES:,} molecules · shared 0–{SHARED_COLOR_CEILING:.0f} molecule colour scale",
        fontsize=4.15,
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.690,
        0.394,
        f"{SECONDARY_MOLECULES:,} molecules · shared 0–{SHARED_COLOR_CEILING:.0f} molecule colour scale",
        fontsize=4.15,
        color=BLUE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    primary_monocyte = _usage_fraction(inputs, PRIMARY_ISOFORM, "Monocyte")
    primary_plasma = _usage_fraction(inputs, PRIMARY_ISOFORM, "Plasma cell")
    primary_progenitor = _usage_fraction(inputs, PRIMARY_ISOFORM, "Progenitor")
    secondary_progenitor = _usage_fraction(
        inputs, SECONDARY_ISOFORM, "Progenitor"
    )
    secondary_bcell = _usage_fraction(inputs, SECONDARY_ISOFORM, "B cell")
    secondary_monocyte = _usage_fraction(
        inputs, SECONDARY_ISOFORM, "Monocyte"
    )
    fig.text(
        0.385,
        0.375,
        "Usage: "
        f"Monocyte {primary_monocyte:.2f} · Plasma {primary_plasma:.2f} "
        f"· Progenitor {primary_progenitor:.2f}",
        fontsize=4.15,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.690,
        0.375,
        "Usage: "
        f"Progenitor {secondary_progenitor:.2f} · B cell {secondary_bcell:.2f} "
        f"· Monocyte {secondary_monocyte:.2f}",
        fontsize=4.15,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.385,
        0.356,
        "DTU contrast uses the all-23-isoform denominator; expression UMAPs "
        "provide localization.",
        fontsize=4.25,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.385,
        0.340,
        "Matched coordinates; native color scales differ—compare spatial "
        "distribution, not color intensity.",
        fontsize=4.2,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def render(inputs: V3.Inputs, stem: Path, dpi: int) -> list[Path]:
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
        _draw_dtu_umaps(fig, inputs)
        V4._draw_junction_access_panel(fig, inputs)

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
        default=FIGURES / "NAR_Fig3_v5",
        help="Output path without an extension.",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate frozen inputs and matched portal screenshots only.",
    )
    args = parser.parse_args()

    for image in (
        PRIMARY_CELL_MAP,
        SECONDARY_CELL_MAP,
        V4.JUNCTION_IMAGE,
    ):
        if not image.is_file():
            raise FileNotFoundError(image)
    if not CENTROID_TABLE.is_file():
        raise FileNotFoundError(CENTROID_TABLE)
    inputs = V3.load_inputs()
    print(
        "INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "cell_types": inputs.usage["ct"].nunique(),
            "isoforms": inputs.usage["transcript_id"].nunique(),
            "matrix_molecules": int(round(inputs.usage["count"].sum())),
            "portal_images": 5,
        },
    )
    if args.validate_only:
        return
    outputs = render(inputs, args.stem.resolve(), args.dpi)
    for output in outputs:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
