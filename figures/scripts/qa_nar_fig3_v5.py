#!/usr/bin/env python3
"""Independent data, portal-source and delivery QA for NAR Figure 3 v5."""

from __future__ import annotations

import argparse
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image, ImageChops

import qa_nar_fig3_v3 as Q3


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
ASSETS = FIGURES / "assets"
WALKTHROUGH = FIGURES / "website_walkthrough" / "ptprc_views_v3"
PRIMARY_WALKTHROUGH = (
    FIGURES / "website_walkthrough" / "ptprc_views_v4_isoform"
)
SECONDARY_WALKTHROUGH = (
    FIGURES / "website_walkthrough" / "ptprc_views_v5_dtu_umap"
)

PRIMARY_CELL_MAP = PRIMARY_WALKTHROUGH / "02_cell_map.png"
SECONDARY_CELL_MAP = SECONDARY_WALKTHROUGH / "02_cell_map.png"
CENTROID_TABLE = (
    PROJECT / "tables" / "PTPRC_umap_celltype_centroids_live_20260727.tsv"
)

PORTAL_SOURCES = {
    ASSETS / "ptprc_gene_card_live_v3.png": (
        "6ae052751cab44835d4d7680edd02ab96bba1680c74e6c5b1128bdc5a4467d04",
        (2280, 580),
    ),
    PRIMARY_CELL_MAP: (
        "91c71259a8ae9a3a09c1362f1a4393571458637ebfa9a3f5c91ed1ab7503e17f",
        (2280, 1522),
    ),
    SECONDARY_CELL_MAP: (
        "00a9ef15a6bd03cc5bd0a09d6e5ac8eccb90cd1e31f7a862898569dc11efb552",
        (2280, 1522),
    ),
    WALKTHROUGH / "01_junctions.png": (
        "ec75cff4eaebffeb898dda8bf1ffd24af79a097e425d354dd32f3c2d79b61fbd",
        (2280, 1498),
    ),
}


def _check_cell_map_metadata(
    directory: Path,
    *,
    isoform: str,
    molecules: str,
    scale: str,
) -> Path:
    metadata_path = directory / "ptprc_all_views_metadata.json"
    Q3.require(metadata_path.is_file(), f"Missing metadata: {metadata_path}")
    metadata = json.loads(metadata_path.read_text())
    Q3.require(
        metadata.get("cellmap_right_view") == "isoform",
        f"Cell map mode changed for {isoform}",
    )
    Q3.require(
        metadata.get("cellmap_isoform") == isoform,
        f"Unexpected selected isoform: {metadata.get('cellmap_isoform')}",
    )
    views = metadata.get("views", [])
    Q3.require(
        len(views) == 1 and views[0].get("view") == "cellmap",
        f"Unexpected capture views for {isoform}",
    )
    config = views[0].get("cellmap_config", {})
    Q3.require(
        config.get("right_view") == "isoform",
        f"Right view is not isoform expression for {isoform}",
    )
    Q3.require(
        config.get("isoform_selection", {})
        .get("selected_text", "")
        .startswith(f"{isoform} · {molecules} molecules"),
        f"Selected molecule count changed for {isoform}",
    )
    Q3.require(
        scale in config.get("legend_text", ""),
        f"Native legend scale changed for {isoform}",
    )
    Q3.require(
        "71,913 sampled cells" in views[0].get("text_excerpt", ""),
        f"Sampled-cell scope changed for {isoform}",
    )
    return metadata_path


def check_portal_sources() -> list[Path]:
    for path, (expected_hash, expected_size) in PORTAL_SOURCES.items():
        Q3.require(path.is_file(), f"Missing portal source: {path}")
        Q3.require(
            Q3.sha256(path) == expected_hash,
            f"Portal source changed: {path}",
        )
        with Image.open(path) as image:
            Q3.require(
                image.size == expected_size,
                f"Unexpected portal source dimensions for {path.name}: "
                f"{image.size}",
            )

    primary_metadata = _check_cell_map_metadata(
        PRIMARY_WALKTHROUGH,
        isoform="ENST00000367364",
        molecules="57,185",
        scale="0 → 5 long-read molecules",
    )
    secondary_metadata = _check_cell_map_metadata(
        SECONDARY_WALKTHROUGH,
        isoform="ENST00000697630",
        molecules="25,697",
        scale="0 → 3 long-read molecules",
    )

    with Image.open(PRIMARY_CELL_MAP) as primary:
        primary_reference = primary.crop((0, 220, 1140, 1150)).convert("RGB")
    with Image.open(SECONDARY_CELL_MAP) as secondary:
        secondary_reference = secondary.crop((0, 220, 1140, 1150)).convert(
            "RGB"
        )
    difference = ImageChops.difference(
        primary_reference,
        secondary_reference,
    )
    Q3.require(
        difference.getbbox() is None,
        "Cell-type reference UMAP differs between the two isoform captures",
    )
    print(
        "PORTAL SOURCE ASSERTIONS PASS",
        {
            "sources": len(PORTAL_SOURCES),
            "matched_cell_coordinates": True,
            "native_scales": ["0–5", "0–3"],
        },
    )
    return [primary_metadata, secondary_metadata]


def check_centroid_table() -> None:
    Q3.require(CENTROID_TABLE.is_file(), f"Missing centroid table: {CENTROID_TABLE}")
    with CENTROID_TABLE.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {
        "B cell": (5370, 0.4459564392529185, 0.09396847355421717),
        "Monocyte": (3057, 0.40947345303695687, 0.38095863791857015),
        "Plasma cell": (1678, 0.26724975338376333, 0.30693738908384044),
        "Progenitor": (2972, 0.3117128246802273, 0.45451336423591415),
    }
    Q3.require(
        {row["cell_type"] for row in rows} == set(expected),
        "Unexpected cell types in centroid table",
    )
    for row in rows:
        n_cells, x_fraction, y_fraction = expected[row["cell_type"]]
        Q3.require(
            int(row["n_cells"]) == n_cells,
            f"Cell count changed for {row['cell_type']}",
        )
        Q3.require(
            abs(float(row["canvas_x_fraction"]) - x_fraction) < 1e-15,
            f"Canvas x changed for {row['cell_type']}",
        )
        Q3.require(
            abs(float(row["canvas_y_from_top_fraction"]) - y_fraction) < 1e-15,
            f"Canvas y changed for {row['cell_type']}",
        )
    print("CENTROID TABLE ASSERTIONS PASS", {"cell_types": len(rows)})


def check_exports(stem: Path, dpi: int) -> list[Path]:
    pdf_path = stem.with_suffix(".pdf")
    svg_path = stem.with_suffix(".svg")
    png_path = stem.with_suffix(".png")
    outputs = [pdf_path, svg_path, png_path]
    for path in outputs:
        Q3.require(
            path.is_file() and path.stat().st_size > 0,
            f"Missing export: {path}",
        )

    svg_root = ET.parse(svg_path).getroot()
    image_count = sum(
        1 for node in svg_root.iter() if node.tag.endswith("image")
    )
    text_count = sum(
        1 for node in svg_root.iter() if node.tag.endswith("text")
    )
    Q3.require(
        image_count == 5,
        f"Expected five portal rasters, found {image_count}",
    )
    Q3.require(
        text_count >= 40,
        f"Too few editable SVG text elements: {text_count}",
    )

    pdf_text = Q3.command_output("pdftotext", str(pdf_path), "-")
    for required in (
        "PTPRC",
        "ENSG00000081237",
        "71,913",
        "ENST00000367364",
        "ENST00000697630",
        "0–5 long-read molecules",
        "0–3 long-read molecules",
        "CELL-TYPE DTU",
        "DTU contrast",
        "Matched coordinates",
        "compare spatial distribution, not color intensity",
        "q = 0.00493",
        "effect = 0.552",
        "98,445",
        "all-23-isoform denominator",
        "Labels mark within-type median coordinates",
        "no trajectory inference",
        "317",
        "8,994",
        "130,951",
        "2,341,574",
        "0.00348",
        "q = 1.00",
        "GET /api/gene/PTPRC/overview",
        "?species=human",
        "tables = snapshot = live API",
        "first 100 rows",
    ):
        Q3.require(
            required in pdf_text,
            f"Required vector text is missing: {required}",
        )
    for forbidden in (
        "lineage_isoform_DIU",
        '"export"',
        "RA / RO / RB",
        "GPT",
        "stage-specific DTU",
        "Top five cover",
    ):
        Q3.require(
            forbidden not in pdf_text,
            f"Legacy/pseudo content remains: {forbidden}",
        )

    pdf_info = Q3.command_output("pdfinfo", str(pdf_path))
    Q3.require(
        "Pages:           1" in pdf_info,
        "Figure PDF must contain one page",
    )
    match = re.search(
        r"Page size:\s+([0-9.]+) x ([0-9.]+) pts",
        pdf_info,
    )
    Q3.require(match is not None, "Could not parse PDF page size")
    width_pt, height_pt = map(float, match.groups())
    Q3.require(
        abs(width_pt - 518.740) < 0.1,
        f"PDF width is not 183 mm: {width_pt}",
    )
    Q3.require(
        abs(height_pt - 467.717) < 0.1,
        f"PDF height is not 165 mm: {height_pt}",
    )

    font_info = Q3.command_output("pdffonts", str(pdf_path))
    Q3.require("ArialMT" in font_info, "Arial Regular is not embedded")
    Q3.require(
        "Arial-BoldMT" in font_info,
        "Arial Bold is not embedded",
    )
    Q3.require("Type 3" not in font_info, "Type 3 fonts are not allowed")

    expected_png = (
        int(183.0 / 25.4 * dpi),
        int(165.0 / 25.4 * dpi),
    )
    with Image.open(png_path) as image:
        Q3.require(
            image.size == expected_png,
            f"Unexpected PNG dimensions: {image.size}",
        )
        Q3.require(image.mode == "RGBA", f"PNG is not RGBA: {image.mode}")
        image_dpi = image.info.get("dpi")
        Q3.require(image_dpi is not None, "PNG DPI metadata is missing")
        Q3.require(
            all(abs(value - dpi) < 0.1 for value in image_dpi),
            f"PNG DPI differs from {dpi}: {image_dpi}",
        )

    print("PDF INFO")
    print(pdf_info.strip())
    print("PDF FONTS")
    print(font_info.strip())
    print(
        "DELIVERY ASSERTIONS PASS",
        {
            "svg_images": image_count,
            "svg_editable_text": text_count,
            "png_px": expected_png,
            "dpi": dpi,
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v5",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--live-url",
        default=None,
        help="Optional live overview endpoint for no-proxy parity checking.",
    )
    parser.add_argument(
        "--live-json",
        type=Path,
        default=None,
        help="Fresh live response captured on a network-reachable host.",
    )
    args = parser.parse_args()

    Q3.check_data_and_live_api(args.live_url, args.live_json)
    metadata_files = check_portal_sources()
    check_centroid_table()
    outputs = check_exports(args.stem.resolve(), args.dpi)

    print("SHA256")
    provenance = [
        *outputs,
        *PORTAL_SOURCES.keys(),
        *metadata_files,
        CENTROID_TABLE,
        ASSETS / "ptprc_live_capture_v3.json",
        Path(__file__).with_name("capture_ptprc_all_views.py"),
        Path(__file__).with_name("render_nar_fig3_v3.py"),
        Path(__file__).with_name("render_nar_fig3_v4.py"),
        Path(__file__).with_name("render_nar_fig3_v5.py"),
        Path(__file__).with_name("qa_nar_fig3_v3.py"),
        Path(__file__),
        Path(__file__).with_name("run_nar_fig3_v5.slurm"),
    ]
    if args.live_json:
        provenance.append(args.live_json.resolve())
    for path in provenance:
        Q3.require(path.is_file(), f"Missing provenance file: {path}")
        print(Q3.sha256(path), path)
    print("ALL NAR FIGURE 3 V5 QA CHECKS PASS")


if __name__ == "__main__":
    main()
