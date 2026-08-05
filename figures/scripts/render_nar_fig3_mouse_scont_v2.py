#!/usr/bin/env python3
"""Render the webpage-led CRA044500 mouse use case for NAR Figure 3.

Figure contract
---------------
Core conclusion
    The scTHREAD portal connects real E6.5 cell context to two complementary
    Malat1 isoform-expression maps and a descriptive cell-type usage contrast,
    while a separate Tnrc6c case illustrates stage-associated isoform change.
Archetype
    Asymmetric mixed-modality figure with real portal screenshots as the hero.
Evidence boundary
    UMAPs localize isoform expression; usage fractions are descriptive.
    Cell-bootstrap pseudo-replicates are not independent embryos, so no
    biological-replicate P value or FDR is shown.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

import render_nar_fig3_mouse_scont as base


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
TABLES = PROJECT / "tables"
WEB = PROJECT.parent / "web"
ASE_ROOT = base.ASE_ROOT

WIDTH_MM = 183.0
HEIGHT_MM = 190.0
ISO1 = "ENSMUST00000245150"
ISO2 = "ENSMUST00000172812"
MALAT1_GENE = "ENSMUSG00000092341"

ISO1_PNG = (
    FIGURES
    / "website_walkthrough/mouse_scont_malat1_iso1/02_cell_map.png"
)
ISO2_PNG = (
    FIGURES
    / "website_walkthrough/mouse_scont_malat1_iso2/02_cell_map.png"
)
ISO1_META = (
    FIGURES
    / "website_walkthrough/mouse_scont_malat1_iso1/"
    "mouse_scont_malat1_iso1_metadata.json"
)
ISO2_META = (
    FIGURES
    / "website_walkthrough/mouse_scont_malat1_iso2/"
    "mouse_scont_malat1_iso2_metadata.json"
)
PORTAL_DTU = TABLES / "Fig3_mouse_scONT_portal_DTU_examples.tsv"
PORTAL_AUDIT = WEB / "data/mouse_scont_umap_audit.json"
PORTAL_POINTS = WEB / "data/mouse_scont_umap_points.parquet"
PORTAL_ISOFORMS = WEB / "data/mouse_scont_umap_isoform_expression.parquet"


def _load_web_inputs() -> dict[str, object]:
    paths = [
        ISO1_PNG,
        ISO2_PNG,
        ISO1_META,
        ISO2_META,
        PORTAL_DTU,
        PORTAL_AUDIT,
        PORTAL_POINTS,
        PORTAL_ISOFORMS,
    ]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    metadata = [json.loads(path.read_text()) for path in (ISO1_META, ISO2_META)]
    for payload, isoform in zip(metadata, (ISO1, ISO2)):
        config = payload["views"][0]["cellmap_config"]
        base._require(payload["gene"] == "Malat1", "Unexpected screenshot gene")
        base._require(payload["species"] == "mouse", "Unexpected screenshot species")
        base._require(payload["cellmap_source"] == "E6.5", "Screenshot is not E6.5")
        base._require(
            config["isoform_selection"]["selected_value"] == isoform,
            "Unexpected screenshot isoform",
        )
        base._require(
            config["load_state"].startswith("rendered:mouse:"),
            "Screenshot did not finish rendering",
        )

    isoform_dtu = pd.read_csv(PORTAL_DTU, sep="\t")
    malat1 = isoform_dtu.loc[isoform_dtu["gene_id"].eq(MALAT1_GENE)].copy()
    base._require(set(malat1["isoform_id"]) == {ISO1, ISO2}, "Missing Malat1 DTU rows")
    portal_audit = json.loads(PORTAL_AUDIT.read_text())
    base._require(portal_audit["cells"] == 25_621, "Unexpected portal UMAP cells")
    base._require(portal_audit["stages"]["E6.5"] == 8_367, "Unexpected E6.5 cells")
    return {
        "metadata": metadata,
        "malat1": malat1,
        "max_usage_difference": float(malat1["prop_diff"].abs().max()),
        "portal_audit": portal_audit,
        "paths": paths,
    }


def _draw_compact_cohort(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    base._panel_label(ax, "a", x=-0.035)
    base._panel_title(ax, "Mouse scONT query context")
    box = FancyBboxPatch(
        (0.0, 0.04),
        1.0,
        0.84,
        boxstyle="round,pad=0.015,rounding_size=0.025",
        facecolor=base.PALE_BLUE,
        edgecolor="#BDD3E1",
        linewidth=0.6,
    )
    ax.add_patch(box)
    base._draw_mouse(ax, (0.16, 0.66), 0.82)
    ax.text(0.30, 0.72, "CRA044500", fontsize=7, fontweight="bold", va="center")
    ax.text(0.30, 0.59, "published mouse gastrulation · ONT", fontsize=5.5, color=base.SLATE)
    ax.text(0.08, 0.43, "6 runs · E6.5 / E7.5 / E8.5", fontsize=5.6, fontweight="bold")
    ax.text(0.08, 0.31, "68,417 registered cells", fontsize=5.3, color=base.BLUE)
    ax.text(0.08, 0.20, "25,621-cell analysis embedding", fontsize=5.3, color=base.GREEN)
    ax.text(0.08, 0.09, "Portal query: mouse · Malat1 · E6.5", fontsize=5.0, color=base.SLATE)


def _draw_portal(ax: plt.Axes, path: Path, letter: str, title: str) -> None:
    image = mpimg.imread(path)
    cropped = image[:1245, :, :]
    ax.imshow(cropped)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#AAB2BC")
        spine.set_linewidth(0.55)
    base._panel_label(ax, letter, x=-0.018)
    base._panel_title(ax, title)


def _draw_malat1_dtu(ax: plt.Axes, data: dict[str, object]) -> None:
    rows = data["malat1"].set_index("isoform_id").loc[[ISO1, ISO2]]
    groups = ["Neuroectoderm", "Embryonic\nepiblast"]
    values = np.array(
        [
            rows["prop_Neuroectoderm"].to_numpy(),
            rows["prop_Epiblast"].to_numpy(),
        ]
    )
    x = np.arange(2)
    width = 0.34
    colors = [base.BLUE, base.ORANGE]
    for index, isoform in enumerate((ISO1, ISO2)):
        ax.bar(
            x + (index - 0.5) * width,
            values[:, index],
            width=width,
            color=colors[index],
            edgecolor=base.INK,
            linewidth=0.35,
            label=f"…{isoform[-6:]}",
        )
    ax.set_xticks(x, groups)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Isoform usage fraction")
    base._panel_label(ax, "d")
    base._panel_title(ax, "Descriptive Malat1 usage contrast at E6.5")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=4.8)
    ax.text(
        0.02,
        0.95,
        f"max |Δ usage| = {data['max_usage_difference']:.3f}",
        transform=ax.transAxes,
        fontsize=4.8,
        color=base.PURPLE,
        fontweight="bold",
        ha="left",
        va="top",
    )
    ax.text(
        0.02,
        -0.23,
        "Cell-bootstrap pseudo-replicates are not biological replicates; no P/FDR shown",
        transform=ax.transAxes,
        fontsize=4.3,
        color=base.SLATE,
        ha="left",
        va="top",
    )


def _draw_trace(ax: plt.Axes, data: dict[str, object]) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    base._panel_label(ax, "f", x=-0.035)
    base._panel_title(ax, "Portal-to-evidence traceability")
    base._rounded(ax, (0.0, 0.02), 1.0, 0.82, "#222936", edgecolor="#222936")
    entries = [
        ("WEB UMAP", "8,367 E6.5 cells · matched coordinates", base.BLUE),
        ("ISOFORM 1", "2,186.646 count units, all stages · E6.5 scale 0–1.347", base.SKY),
        ("ISOFORM 2", "514.452 count units, all stages · E6.5 scale 0–1", base.ORANGE),
        ("USAGE VIEW", "descriptive max |Δ usage| 0.245", base.PURPLE),
    ]
    for y, (label, text, color) in zip([0.69, 0.52, 0.35, 0.18], entries):
        ax.add_patch(plt.Circle((0.06, y), 0.022, color=color))
        ax.text(0.11, y + 0.025, label, fontsize=4.3, color="#AEB8C6", fontweight="bold")
        ax.text(0.11, y - 0.025, text, fontsize=4.8, color=base.WHITE)
    ax.text(
        0.98,
        0.06,
        "UMAP = localization · usage = descriptive source proportions",
        fontsize=4.3,
        color="#8FD1C6",
        fontweight="bold",
        ha="right",
    )


def _write_manifest(
    data: dict[str, object],
    figure_name: str,
    role: str,
) -> Path:
    inputs = {
        str(path): {"sha256": base._sha256(path)}
        for path in data["paths"]
    }
    inputs[str(base.TEMPORAL_CASES)] = {"sha256": base._sha256(base.TEMPORAL_CASES)}
    manifest = {
        "figure": figure_name,
        "role": role,
        "claim": (
            "Real E6.5 portal UMAPs localize two Malat1 isoforms; source usage "
            "fractions provide a descriptive cell-type contrast."
        ),
        "screenshots": {
            ISO1: {
                "png": str(ISO1_PNG),
                "metadata": str(ISO1_META),
                "source": "E6.5",
            },
            ISO2: {
                "png": str(ISO2_PNG),
                "metadata": str(ISO2_META),
                "source": "E6.5",
            },
        },
        "guardrails": [
            "UMAP localizes expression and does not itself test differential usage",
            "the two portal screenshots use the same E6.5 coordinates",
            "cell-bootstrap pseudo-replicates are not independent embryos",
            "no biological-replicate P value or FDR is shown",
            "Tnrc6c stage trajectories are descriptive proportions",
        ],
        "inputs": inputs,
    }
    manifest_name = (
        "Fig3_mouse_scONT_v2_manifest.json"
        if figure_name == "NAR_Fig3_mouse_scONT_v2"
        else f"{figure_name}_manifest.json"
    )
    path = TABLES / manifest_name
    path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def render(
    stem: Path,
    dpi: int,
    suptitle: str = "",
    role: str = (
        "webpage-led use of published NGDC CRA044500 mouse data; "
        "preferred PTPRC v6 unchanged"
    ),
) -> tuple[list[Path], Path]:
    family = base._configure_style()
    base_data = base._load_inputs()
    web_data = _load_web_inputs()
    manifest = _write_manifest(web_data, stem.name, role)

    fig = plt.figure(
        figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4),
        facecolor="none",
    )
    if suptitle:
        fig.text(
            0.055,
            0.988,
            suptitle,
            ha="left",
            va="top",
            fontsize=8.6,
            fontweight="bold",
            color=base.INK,
        )
    if suptitle:
        ax_a = fig.add_axes([0.055, 0.655, 0.255, 0.245])
        ax_b = fig.add_axes([0.350, 0.625, 0.615, 0.275])
    else:
        ax_a = fig.add_axes([0.055, 0.690, 0.255, 0.255])
        ax_b = fig.add_axes([0.350, 0.665, 0.615, 0.285])
    ax_c = fig.add_axes([0.055, 0.335, 0.615, 0.285])
    ax_d = fig.add_axes([0.735, 0.375, 0.220, 0.205])
    ax_e = fig.add_axes([0.080, 0.075, 0.365, 0.175])
    ax_f = fig.add_axes([0.535, 0.055, 0.420, 0.205])

    _draw_compact_cohort(ax_a)
    _draw_portal(
        ax_b,
        ISO1_PNG,
        "b",
        f"Live portal · {ISO1} expression",
    )
    _draw_portal(
        ax_c,
        ISO2_PNG,
        "c",
        f"Matched portal view · {ISO2} expression",
    )
    _draw_malat1_dtu(ax_d, web_data)
    base._draw_tnrc6c(ax_e, base_data)
    ax_e.texts[0].set_text("e")
    _draw_trace(ax_f, web_data)

    for ax in (ax_d, ax_e):
        ax.tick_params(length=2.0, width=0.5, pad=1.5)
        for spine in ax.spines.values():
            spine.set_color(base.INK)
            spine.set_linewidth(0.5)

    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [
        stem.with_suffix(".pdf"),
        stem.with_suffix(".svg"),
        stem.with_suffix(".png"),
    ]
    # Isolate the fixed-canvas export from a caller's global Matplotlib
    # configuration (render_nar_sf_all imports nar_style, whose default is
    # savefig.bbox="tight"). The delivery contract is exactly 183 × 190 mm.
    with plt.rc_context({"savefig.bbox": None}):
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
                "manifest": str(manifest),
            },
            indent=2,
        )
    )
    return outputs, manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_mouse_scONT_v2",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--suptitle", default="")
    parser.add_argument(
        "--role",
        default=(
            "webpage-led use of published NGDC CRA044500 mouse data; "
            "preferred PTPRC v6 unchanged"
        ),
    )
    args = parser.parse_args()
    render(
        args.stem.resolve(),
        args.dpi,
        suptitle=args.suptitle,
        role=args.role,
    )


if __name__ == "__main__":
    main()
