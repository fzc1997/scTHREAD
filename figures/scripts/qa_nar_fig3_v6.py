#!/usr/bin/env python3
"""Independent data, source and delivery QA for six-panel NAR Figure 3 v6."""

from __future__ import annotations

import argparse
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image, ImageChops

import qa_nar_fig3_v3 as Q3
import qa_nar_fig3_v5 as Q5


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
ASSETS = FIGURES / "assets"
JUNCTION_SOURCE = (
    FIGURES / "website_walkthrough" / "ptprc_views_v6_junction_vector"
)
JUNCTION_TABLE = JUNCTION_SOURCE / "junction_first100_live.tsv"
JUNCTION_METADATA = JUNCTION_SOURCE / "ptprc_all_views_metadata.json"
JUNCTION_REFERENCE = JUNCTION_SOURCE / "01_junctions.png"
JUNCTION_README = JUNCTION_SOURCE / "README.md"

EXPECTED_JUNCTION_HASHES = {
    JUNCTION_TABLE: (
        "5f794b732ae0f284c7241c08b79984062a306af9d7c7b579006a56ceae3495d2"
    ),
    JUNCTION_METADATA: (
        "7afec3fb62deed00be67aea87cd03a08b55018b81ed304d8a9bf41609ec72729"
    ),
    JUNCTION_REFERENCE: (
        "0663a9e8beee3837e70f8d19d46a5a683e0375b4453e19690a8270c52f6ace90"
    ),
}

WIDTH_MM = 183.0
HEIGHT_MM = 180.0


def check_used_portal_sources() -> tuple[list[Path], list[Path]]:
    used_rasters = [
        ASSETS / "ptprc_gene_card_live_v3.png",
        Q5.PRIMARY_CELL_MAP,
        Q5.SECONDARY_CELL_MAP,
    ]
    for path in used_rasters:
        expected_hash, expected_size = Q5.PORTAL_SOURCES[path]
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

    primary_metadata = Q5._check_cell_map_metadata(
        Q5.PRIMARY_WALKTHROUGH,
        isoform="ENST00000367364",
        molecules="57,185",
        scale="0 → 5 long-read molecules",
    )
    secondary_metadata = Q5._check_cell_map_metadata(
        Q5.SECONDARY_WALKTHROUGH,
        isoform="ENST00000697630",
        molecules="25,697",
        scale="0 → 3 long-read molecules",
    )
    with Image.open(Q5.PRIMARY_CELL_MAP) as primary:
        primary_reference = primary.crop((0, 220, 1140, 1150)).convert("RGB")
    with Image.open(Q5.SECONDARY_CELL_MAP) as secondary:
        secondary_reference = secondary.crop(
            (0, 220, 1140, 1150)
        ).convert("RGB")
    difference = ImageChops.difference(primary_reference, secondary_reference)
    Q3.require(
        difference.getbbox() is None,
        "Cell-type reference UMAP differs between the two isoform captures",
    )
    print(
        "USED PORTAL SOURCE ASSERTIONS PASS",
        {
            "rasters": len(used_rasters),
            "matched_cell_coordinates": True,
            "native_scales": ["0–5", "0–3"],
        },
    )
    return used_rasters, [primary_metadata, secondary_metadata]


def check_junction_table() -> tuple[list[Path], dict[str, str]]:
    for path, expected_hash in EXPECTED_JUNCTION_HASHES.items():
        Q3.require(path.is_file(), f"Missing junction source: {path}")
        Q3.require(
            Q3.sha256(path) == expected_hash,
            f"Junction source changed: {path}",
        )
    Q3.require(JUNCTION_README.is_file(), "Missing junction source README")
    with Image.open(JUNCTION_REFERENCE) as image:
        Q3.require(
            image.size == (4560, 2996),
            f"Unexpected 4× junction reference size: {image.size}",
        )

    metadata = json.loads(JUNCTION_METADATA.read_text())
    Q3.require(metadata.get("gene") == "PTPRC", "Unexpected junction gene")
    Q3.require(
        metadata.get("device_scale_factor") == 4.0,
        "Junction reference is not the 4× capture",
    )
    views = metadata.get("views", [])
    Q3.require(
        len(views) == 1 and views[0].get("view") == "junctions",
        "Unexpected junction capture views",
    )
    Q3.require(
        views[0].get("junction_table_rows") == 100,
        "Live junction export no longer has 100 rows",
    )
    Q3.require(
        "317 junctions" in metadata.get("browser_location", ""),
        "Filtered junction scope changed",
    )

    with JUNCTION_TABLE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        Q3.require(
            reader.fieldnames
            == ["junction", "span", "molecules", "reads", "runs", "studies"],
            f"Unexpected junction columns: {reader.fieldnames}",
        )
    Q3.require(len(rows) == 100, f"Expected 100 junction rows, found {len(rows)}")
    Q3.require(
        len({row["junction"] for row in rows}) == 100,
        "Junction identifiers are not unique",
    )
    pattern = re.compile(
        r"^chr(?P<chrom>[^:]+):(?P<strand>[+-]):"
        r"(?P<start>[0-9]+)-(?P<end>[0-9]+)$"
    )
    molecules: list[int] = []
    strands: set[str] = set()
    boundary_crossing = 0
    display_candidates: list[dict[str, str]] = []
    for row in rows:
        match = pattern.fullmatch(row["junction"])
        Q3.require(match is not None, f"Unparseable junction: {row['junction']}")
        start = int(match.group("start"))
        end = int(match.group("end"))
        span = int(row["span"].replace(",", "").removesuffix(" bp"))
        Q3.require(end - start == span, f"Span mismatch: {row['junction']}")
        Q3.require(
            start < end
            and end > 198_638_457
            and start < 198_757_476,
            f"Junction does not overlap PTPRC view: {row['junction']}",
        )
        boundary_crossing += int(
            start < 198_638_457 or end > 198_757_476
        )
        strands.add(match.group("strand"))
        values = [
            int(row[column].replace(",", ""))
            for column in ("molecules", "reads", "runs", "studies")
        ]
        Q3.require(all(value > 0 for value in values), "Non-positive support")
        molecules.append(values[0])
        if span >= 10_000:
            display_candidates.append(row)
    Q3.require(strands == {"+", "-"}, f"Unexpected strands: {strands}")
    Q3.require(
        boundary_crossing == 1,
        f"Unexpected boundary-crossing junction count: {boundary_crossing}",
    )
    Q3.require(
        molecules == sorted(molecules, reverse=True),
        "Junction rows are no longer support-sorted",
    )
    Q3.require(
        display_candidates,
        "No junction is long enough for the display example",
    )
    selected = sorted(
        display_candidates,
        key=lambda row: (
            -int(row["molecules"].replace(",", "")),
            -int(row["span"].replace(",", "").removesuffix(" bp")),
            row["junction"],
        ),
    )[0]
    Q3.require(
        selected
        == {
            "junction": "chr1:+:198639341-198692347",
            "span": "53,006 bp",
            "molecules": "201,317",
            "reads": "518,892",
            "runs": "47",
            "studies": "8",
        },
        f"Unexpected selected display junction: {selected}",
    )
    print(
        "VECTOR JUNCTION SOURCE ASSERTIONS PASS",
        {
            "rows": len(rows),
            "strands": sorted(strands),
            "boundary_crossing": boundary_crossing,
            "molecule_range": [min(molecules), max(molecules)],
            "selected_display_junction": selected["junction"],
        },
    )
    return (
        [
            JUNCTION_TABLE,
            JUNCTION_METADATA,
            JUNCTION_REFERENCE,
            JUNCTION_README,
        ],
        selected,
    )


def check_exports(
    stem: Path,
    dpi: int,
    selected: dict[str, str],
) -> list[Path]:
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
        image_count == 4,
        f"Expected four portal rasters, found {image_count}",
    )
    Q3.require(
        text_count >= 95,
        f"Too few editable SVG text elements: {text_count}",
    )

    pdf_text = Q3.command_output("pdftotext", str(pdf_path), "-")
    for required in (
        "PTPRC",
        "ENSG00000081237",
        "Shared cell-type UMAP",
        "ENST00000367364 localization",
        "ENST00000697630 localization",
        "Usage fractions confirm cell-type DTU",
        "Junction evidence remains inspectable and exportable",
        "71,913",
        "57,185",
        "25,697",
        "0–5 long-read molecules",
        "0–3 long-read molecules",
        "Expression UMAPs localize the contrast; panel e quantifies DTU",
        "compare spatial distribution, not native color intensity",
        "CELL-TYPE DTU",
        "q = 0.00493",
        "effect = 0.552",
        "all-23-isoform denominator",
        "98,445",
        "Usage fraction",
        "0.59",
        "0.40",
        "Fractions are descriptive; q/effect summarize the gene-level test",
        "no trajectory inference",
        "317",
        "100 rows shown / 317 filtered",
        "+ strand",
        "− strand",
        "selected example",
        "width / opacity scale with log10(molecules)",
        "chr1 genomic coordinate (Mb)",
        "Gold arc is a display example; selecting any arc exposes",
        "SELECTED ARC · DISPLAY EXAMPLE",
        "chr1:+",
        "198,639,341–198,692,347",
        "MOLECULES",
        selected["molecules"],
        "READS",
        selected["reads"],
        "RUNS",
        selected["runs"],
        "STUDIES",
        selected["studies"],
        "8,994",
        "2,341,574",
        "0.00348",
        "q = 1.00",
        "GET /api/gene/PTPRC/",
        "overview?species=human",
        "row: live table · overview: frozen snapshot",
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
        "tables = snapshot = live API",
        "source tables · frozen API snapshot",
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
    expected_width = WIDTH_MM / 25.4 * 72.0
    expected_height = HEIGHT_MM / 25.4 * 72.0
    Q3.require(
        abs(width_pt - expected_width) < 0.1,
        f"PDF width is not {WIDTH_MM:g} mm: {width_pt}",
    )
    Q3.require(
        abs(height_pt - expected_height) < 0.1,
        f"PDF height is not {HEIGHT_MM:g} mm: {height_pt}",
    )

    font_info = Q3.command_output("pdffonts", str(pdf_path))
    Q3.require("ArialMT" in font_info, "Arial Regular is not embedded")
    Q3.require(
        "Arial-BoldMT" in font_info,
        "Arial Bold is not embedded",
    )
    Q3.require("Type 3" not in font_info, "Type 3 fonts are not allowed")

    expected_png = (
        int(WIDTH_MM / 25.4 * dpi),
        int(HEIGHT_MM / 25.4 * dpi),
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
            "canvas_mm": [WIDTH_MM, HEIGHT_MM],
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v6",
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
    portal_files, metadata_files = check_used_portal_sources()
    junction_files, selected = check_junction_table()
    Q5.check_centroid_table()
    outputs = check_exports(args.stem.resolve(), args.dpi, selected)

    print("SHA256")
    provenance = [
        *outputs,
        *portal_files,
        *metadata_files,
        *junction_files,
        Q5.CENTROID_TABLE,
        ASSETS / "ptprc_live_capture_v3.json",
        Path(__file__).with_name("capture_ptprc_all_views.py"),
        Path(__file__).with_name("render_nar_fig3_v3.py"),
        Path(__file__).with_name("render_nar_fig3_v4.py"),
        Path(__file__).with_name("render_nar_fig3_v5.py"),
        Path(__file__).with_name("render_nar_fig3_v6.py"),
        Path(__file__).with_name("qa_nar_fig3_v3.py"),
        Path(__file__).with_name("qa_nar_fig3_v5.py"),
        Path(__file__),
        Path(__file__).with_name("run_nar_fig3_v6.slurm"),
    ]
    if args.live_json:
        provenance.append(args.live_json.resolve())
    for path in provenance:
        Q3.require(path.is_file(), f"Missing provenance file: {path}")
        print(Q3.sha256(path), path)
    print("ALL NAR FIGURE 3 V6 QA CHECKS PASS")


if __name__ == "__main__":
    main()
