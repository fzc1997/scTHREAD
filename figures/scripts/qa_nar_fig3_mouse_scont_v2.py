#!/usr/bin/env python3
"""Independent source and delivery QA for webpage-led mouse scONT Figure 3."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

from PIL import Image
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
TABLES = PROJECT / "tables"
ISO1 = "ENSMUST00000245150"
ISO2 = "ENSMUST00000172812"
SOURCES = {
    ISO1: (
        FIGURES / "website_walkthrough/mouse_scont_malat1_iso1/02_cell_map.png",
        FIGURES / "website_walkthrough/mouse_scont_malat1_iso1/mouse_scont_malat1_iso1_metadata.json",
        "7b286a60ef422c89d9c96524884cdf247e2916345bf3a7f2bb20cfc1d0d4cd1f",
    ),
    ISO2: (
        FIGURES / "website_walkthrough/mouse_scont_malat1_iso2/02_cell_map.png",
        FIGURES / "website_walkthrough/mouse_scont_malat1_iso2/mouse_scont_malat1_iso2_metadata.json",
        "f0dc1a3fd5c2fb52df1e71c05aa0a9aaad5e53b3dac8825adbcc4480731bd731",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command(*args: str) -> str:
    return subprocess.run(
        args,
        check=True,
        text=True,
        capture_output=True,
    ).stdout


def qa(stem: Path, dpi: int) -> None:
    for isoform, (png, metadata_path, expected_hash) in SOURCES.items():
        require(sha256(png) == expected_hash, f"Screenshot hash changed: {png}")
        metadata = json.loads(metadata_path.read_text())
        config = metadata["views"][0]["cellmap_config"]
        require(metadata["gene"] == "Malat1", "Screenshot gene mismatch")
        require(metadata["species"] == "mouse", "Screenshot species mismatch")
        require(metadata["cellmap_source"] == "E6.5", "Screenshot stage mismatch")
        require(
            config["isoform_selection"]["selected_value"] == isoform,
            "Screenshot isoform mismatch",
        )
        require(
            config["load_state"].startswith("rendered:mouse:"),
            "Screenshot was not captured from a rendered state",
        )
        require(Image.open(png).size == (2280, 1522), "Screenshot dimensions changed")

    dtu = pd.read_csv(TABLES / "Fig3_mouse_scONT_portal_DTU_examples.tsv", sep="\t")
    malat1 = dtu.loc[dtu["gene_id"].eq("ENSMUSG00000092341")]
    require(set(malat1["isoform_id"]) == {ISO1, ISO2}, "Malat1 DTU source mismatch")
    require(
        abs(float(malat1["prop_diff"].abs().max()) - 0.244947468890936) < 1e-12,
        "Malat1 usage effect mismatch",
    )

    pdf, svg, png = (stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png"))
    for path in (pdf, svg, png):
        require(path.is_file() and path.stat().st_size > 10_000, f"Missing output: {path}")

    info = command("pdfinfo", str(pdf))
    require("Pages:           1" in info, "PDF must contain one page")
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", info)
    require(match is not None, "PDF page size missing")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(abs(width_mm - 183.0) < 0.15, f"Unexpected width: {width_mm}")
    require(abs(height_mm - 190.0) < 0.15, f"Unexpected height: {height_mm}")

    fonts = command("pdffonts", str(pdf))
    require("Type 3" not in fonts, "Type 3 font detected")
    require("Arial" in fonts and "TrueType" in fonts, "Arial TrueType not embedded")
    text = command("pdftotext", str(pdf), "-")
    for phrase in (
        "Mouse scONT query context",
        "CRA044500",
        "published mouse gastrulation",
        "68,417 registered cells",
        "Descriptive Malat1 usage contrast at E6.5",
        "UMAP = localization",
        "Tnrc6c reciprocal isoform switch",
    ):
        require(phrase in text, f"PDF text missing: {phrase}")
    for phrase in (
        "gene-level FDR",
        "DTU TEST",
        "Portal-to-inference traceability",
        "significant cell-type usage contrast",
        "internally generated",
    ):
        require(phrase not in text, f"Forbidden inferential wording in PDF: {phrase}")

    svg_text = svg.read_text(encoding="utf-8")
    require(svg_text.count("<image") == 2, "Expected exactly two portal rasters")
    require(svg_text.count("<text") >= 35, "Too few editable SVG text objects")

    with Image.open(png) as image:
        expected = (
            int(183.0 / 25.4 * dpi),
            int(190.0 / 25.4 * dpi),
        )
        require(image.size == expected, f"PNG size {image.size} != {expected}")
        require(image.mode == "RGBA", "PNG must retain transparency")
        dpi_info = image.info.get("dpi")
        require(dpi_info is not None, "PNG DPI missing")
        require(all(abs(value - dpi) < 1 for value in dpi_info), "PNG DPI mismatch")

    manifest_name = (
        "Fig3_mouse_scONT_v2_manifest.json"
        if stem.name == "NAR_Fig3_mouse_scONT_v2"
        else f"{stem.name}_manifest.json"
    )
    manifest = TABLES / manifest_name
    payload = json.loads(manifest.read_text())
    for path_text, record in payload["inputs"].items():
        require(sha256(Path(path_text)) == record["sha256"], f"Input hash mismatch: {path_text}")

    print(
        "DELIVERY ASSERTIONS PASS",
        {
            "canvas_mm": [round(width_mm, 2), round(height_mm, 2)],
            "png_px": Image.open(png).size,
            "svg_rasters": svg_text.count("<image"),
            "svg_editable_text": svg_text.count("<text"),
        },
    )
    print("SHA256")
    for path in (pdf, svg, png, manifest):
        print(sha256(path), path)
    print("ALL NAR FIGURE 3 MOUSE scONT V2 QA CHECKS PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()
    qa(args.stem.resolve(), args.dpi)


if __name__ == "__main__":
    main()
