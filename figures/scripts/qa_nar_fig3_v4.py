#!/usr/bin/env python3
"""Independent data, portal-source and delivery QA for NAR Figure 3 v4."""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image

import qa_nar_fig3_v3 as Q3


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
ASSETS = FIGURES / "assets"
WALKTHROUGH = FIGURES / "website_walkthrough" / "ptprc_views_v3"
ISOFORM_WALKTHROUGH = FIGURES / "website_walkthrough" / "ptprc_views_v4_isoform"

PORTAL_SOURCES = {
    ASSETS / "ptprc_gene_card_live_v3.png": (
        "6ae052751cab44835d4d7680edd02ab96bba1680c74e6c5b1128bdc5a4467d04",
        (2280, 580),
    ),
    ISOFORM_WALKTHROUGH / "02_cell_map.png": (
        "91c71259a8ae9a3a09c1362f1a4393571458637ebfa9a3f5c91ed1ab7503e17f",
        (2280, 1522),
    ),
    WALKTHROUGH / "01_junctions.png": (
        "ec75cff4eaebffeb898dda8bf1ffd24af79a097e425d354dd32f3c2d79b61fbd",
        (2280, 1498),
    ),
}


def check_portal_sources() -> None:
    for path, (expected_hash, expected_size) in PORTAL_SOURCES.items():
        Q3.require(path.is_file(), f"Missing portal source: {path}")
        Q3.require(Q3.sha256(path) == expected_hash, f"Portal source changed: {path}")
        with Image.open(path) as image:
            Q3.require(
                image.size == expected_size,
                f"Unexpected portal source dimensions for {path.name}: {image.size}",
            )
    metadata_path = ISOFORM_WALKTHROUGH / "ptprc_all_views_metadata.json"
    Q3.require(metadata_path.is_file(), f"Missing Cell map metadata: {metadata_path}")
    metadata = json.loads(metadata_path.read_text())
    Q3.require(
        metadata.get("cellmap_right_view") == "isoform",
        "Cell map right view is not isoform expression",
    )
    Q3.require(
        metadata.get("cellmap_isoform") == "ENST00000367364",
        "Unexpected Cell map isoform",
    )
    views = metadata.get("views", [])
    Q3.require(len(views) == 1 and views[0].get("view") == "cellmap", "Unexpected capture views")
    config = views[0].get("cellmap_config", {})
    Q3.require(config.get("right_view") == "isoform", "Cell map mode was not captured")
    Q3.require(
        "0 → 5 long-read molecules" in config.get("legend_text", ""),
        "Cell map isoform legend changed",
    )
    print(
        "PORTAL SOURCE ASSERTIONS PASS",
        {path.name: size for path, (_, size) in PORTAL_SOURCES.items()},
    )


def check_exports(stem: Path, dpi: int) -> list[Path]:
    pdf_path = stem.with_suffix(".pdf")
    svg_path = stem.with_suffix(".svg")
    png_path = stem.with_suffix(".png")
    outputs = [pdf_path, svg_path, png_path]
    for path in outputs:
        Q3.require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    svg_root = ET.parse(svg_path).getroot()
    image_count = sum(1 for node in svg_root.iter() if node.tag.endswith("image"))
    text_count = sum(1 for node in svg_root.iter() if node.tag.endswith("text"))
    Q3.require(image_count == 3, f"Expected three portal rasters, found {image_count}")
    Q3.require(text_count >= 75, f"Too few editable SVG text elements: {text_count}")

    pdf_text = Q3.command_output("pdftotext", str(pdf_path), "-")
    for required in (
        "PTPRC",
        "ENSG00000081237",
        "71,913",
        "ENST00000367364",
        "0–5 long-read molecules",
        "CELL-TYPE DTU",
        "DTU contrast",
        "317",
        "8,994",
        "130,951",
        "2,341,574",
        "0.00493",
        "0.00348",
        "q = 1.00",
        "98,445",
        "values not renormalized",
        "GET /api/gene/PTPRC/overview",
        "?species=human",
        "tables = snapshot = live API",
        "no trajectory inference",
        "first 100 rows",
    ):
        Q3.require(required in pdf_text, f"Required vector text is missing: {required}")
    for forbidden in (
        "lineage_isoform_DIU",
        '"export"',
        "RA / RO / RB",
        "GPT",
        "stage-specific DTU",
    ):
        Q3.require(forbidden not in pdf_text, f"Legacy/pseudo content remains: {forbidden}")

    pdf_info = Q3.command_output("pdfinfo", str(pdf_path))
    Q3.require("Pages:           1" in pdf_info, "Figure PDF must contain one page")
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    Q3.require(match is not None, "Could not parse PDF page size")
    width_pt, height_pt = map(float, match.groups())
    Q3.require(abs(width_pt - 518.740) < 0.1, f"PDF width is not 183 mm: {width_pt}")
    Q3.require(abs(height_pt - 467.717) < 0.1, f"PDF height is not 165 mm: {height_pt}")

    font_info = Q3.command_output("pdffonts", str(pdf_path))
    Q3.require("ArialMT" in font_info, "Arial Regular is not embedded")
    Q3.require("Arial-BoldMT" in font_info, "Arial Bold is not embedded")
    Q3.require("Type 3" not in font_info, "Type 3 fonts are not allowed")

    expected_png = (
        int(183.0 / 25.4 * dpi),
        int(165.0 / 25.4 * dpi),
    )
    with Image.open(png_path) as image:
        Q3.require(image.size == expected_png, f"Unexpected PNG dimensions: {image.size}")
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
        default=FIGURES / "NAR_Fig3_v4",
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
    check_portal_sources()
    outputs = check_exports(args.stem.resolve(), args.dpi)

    print("SHA256")
    provenance = [
        *outputs,
        *PORTAL_SOURCES.keys(),
        ASSETS / "ptprc_live_capture_v3.json",
        ISOFORM_WALKTHROUGH / "ptprc_all_views_metadata.json",
        Path(__file__).with_name("capture_ptprc_all_views.py"),
        Path(__file__).with_name("render_nar_fig3_v3.py"),
        Path(__file__).with_name("render_nar_fig3_v4.py"),
        Path(__file__).with_name("qa_nar_fig3_v3.py"),
        Path(__file__),
        Path(__file__).with_name("run_nar_fig3_v4.slurm"),
    ]
    if args.live_json:
        provenance.append(args.live_json.resolve())
    for path in provenance:
        Q3.require(path.is_file(), f"Missing provenance file: {path}")
        print(Q3.sha256(path), path)
    print("ALL NAR FIGURE 3 V4 QA CHECKS PASS")


if __name__ == "__main__":
    main()
