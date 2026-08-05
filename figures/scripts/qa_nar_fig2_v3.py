#!/usr/bin/env python3
"""Structural and source-contract QA for corrected NAR Figure 2 v3."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

from PIL import Image
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables/p0_biological_unit_rerun"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def output(*args: str) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=ROOT / "figures/NAR_Fig2_v3",
    )
    parser.add_argument(
        "--observed-suffix",
        default="",
        help="Optional suffix for observed tables, for example _9999",
    )
    args = parser.parse_args()

    pdf = args.stem.with_suffix(".pdf")
    svg = args.stem.with_suffix(".svg")
    png = args.stem.with_suffix(".png")
    for path in (pdf, svg, png):
        require(path.is_file() and path.stat().st_size > 0, f"missing {path}")

    pdfinfo = output("pdfinfo", str(pdf))
    require("Pages:           1" in pdfinfo, "PDF is not one page")
    require(
        "Page size:       518.74 x 595.276 pts" in pdfinfo,
        "PDF is not the declared 183 x 210 mm canvas",
    )
    require("Type 3" not in output("pdffonts", str(pdf)), "PDF contains Type 3 font")

    image = Image.open(png)
    require(image.mode == "RGBA", f"PNG lacks alpha channel: {image.mode}")
    dpi = image.info.get("dpi", (0, 0))
    require(all(abs(value - 450) < 1 for value in dpi), f"PNG DPI is {dpi}")
    require(
        abs(image.width - 3242) <= 2 and abs(image.height - 3720) <= 2,
        f"PNG is not 183 x 210 mm at 450 dpi: {image.size}",
    )

    tree = ET.parse(svg)
    root = tree.getroot()
    text_nodes = [node for node in root.iter() if node.tag.endswith("text")]
    require(len(text_nodes) >= 45, f"too few editable SVG text nodes: {len(text_nodes)}")

    expected = {}
    frozen_v3 = {"ASE": (6930, 0), "DIU": (8092, 1971), "APA": (10531, 2527)}
    for analysis in ("ASE", "DIU", "APA"):
        table_path = TABLES / f"{analysis.lower()}_observed{args.observed_suffix}.tsv"
        table = pd.read_csv(table_path, sep="\t")
        expected[analysis] = (len(table), int(table.sig.sum()))
        if args.observed_suffix:
            manifest = json.loads(
                table_path.with_suffix(table_path.suffix + ".manifest.json").read_text()
            )
            require(manifest["genes_tested"] == len(table), f"{analysis} manifest rows")
            require(
                manifest["significant_genes"] == int(table.sig.sum()),
                f"{analysis} manifest significant count",
            )
    if not args.observed_suffix:
        require(expected == frozen_v3, f"v3 source counts changed: {expected}")

    validation = json.loads((TABLES / "validation_report.json").read_text())
    require(validation["status"] == "PASS", "upstream validation is not PASS")
    report = {
        "status": "PASS",
        "figure": str(pdf),
        "pdf_sha256": sha256(pdf),
        "svg_sha256": sha256(svg),
        "png_sha256": sha256(png),
        "png_pixels": [image.width, image.height],
        "png_dpi": list(dpi),
        "svg_editable_text_nodes": len(text_nodes),
        "source_counts": expected,
        "upstream_validation": str(TABLES / "validation_report.json"),
    }
    report_path = ROOT / "figures/_qa" / f"{args.stem.name}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise
