#!/usr/bin/env python3
"""Reproducible delivery checks for the scTHREAD graphical abstract v2."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from PIL import Image


EXPECTED_SYSTEMS = {
    "Blood/immune": (231_103, 40, 3),
    "Neural/sensory": (296_089, 64, 9),
    "Cancer": (136_006, 15, 6),
    "Endocrine": (90_963, 6, 1),
    "Heart/vascular": (68_204, 224, 2),
    "Development/embryo": (20_797, 107, 5),
    "Reproductive": (2_331, 1, 1),
    "Other tissues": (288, 12, 4),
}
REQUIRED_TEXT = (
    "Human ovary",
    "Mouse gastrula",
    "Ovary",
    "Embryo",
    "469 samples",
    "31 studies",
    "845,781",
    ">200k isoforms",
)
FORBIDDEN_VISUAL_TERMS = (
    "benchmark",
    "method",
    "Smart-seq2",
    "differentiation",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(*args: str) -> str:
    return subprocess.run(
        args,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_catalog(component_dir: Path) -> None:
    summary_path = component_dir / "catalog_system_composition_20260726.tsv"
    studies_path = (
        component_dir / "catalog_study_biology_classification_20260726.tsv"
    )
    featured_path = component_dir / "featured_project_cohorts_20260726.tsv"

    summary = pd.read_csv(summary_path, sep="\t")
    observed_systems = set(summary["system"])
    require(
        observed_systems == set(EXPECTED_SYSTEMS),
        f"Unexpected biology categories: {sorted(observed_systems)}",
    )
    for system, (cells, samples, studies) in EXPECTED_SYSTEMS.items():
        row = summary.loc[summary["system"] == system].iloc[0]
        require(int(row["cells"]) == cells, f"{system}: cell total changed")
        require(int(row["samples"]) == samples, f"{system}: sample total changed")
        require(int(row["studies"]) == studies, f"{system}: study total changed")
    require(int(summary["cells"].sum()) == 845_781, "Catalog cell sum changed")
    require(int(summary["samples"].sum()) == 469, "Catalog sample sum changed")
    require(int(summary["studies"].sum()) == 31, "Catalog study sum changed")

    studies = pd.read_csv(studies_path, sep="\t")
    require(len(studies) == 31, "Study classification must contain 31 studies")
    require(
        set(studies["system"]) == set(EXPECTED_SYSTEMS),
        "Study classification and summary categories disagree",
    )
    require(int(studies["cells"].sum()) == 845_781, "Study-level cell sum changed")
    require(int(studies["samples"].sum()) == 469, "Study-level sample sum changed")

    featured = pd.read_csv(featured_path, sep="\t", keep_default_na=False)
    require(
        set(featured["display_label"]) == {"Human ovary", "Mouse gastrula"},
        "Featured biological cohorts changed",
    )
    require(set(featured["assay"]) == {"scONT"}, "Featured assays must be scONT")
    require(
        set(featured["included_in_469_release_denominator"]) == {"no"},
        "Featured cohorts must remain outside the frozen quantitative denominator",
    )
    print(
        "BIOLOGY SOURCE ASSERTIONS PASS",
        {
            "cells": int(summary["cells"].sum()),
            "samples": int(summary["samples"].sum()),
            "studies": int(summary["studies"].sum()),
            "systems": sorted(observed_systems),
            "featured": featured["display_label"].tolist(),
        },
    )


def check_exports(stem: Path, component_dir: Path) -> list[Path]:
    svg_path = stem.with_suffix(".svg")
    pdf_path = stem.with_suffix(".pdf")
    png_path = stem.with_suffix(".png")
    paths = [svg_path, pdf_path, png_path]
    for path in paths:
        require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    tree = ET.parse(svg_path)
    root = tree.getroot()
    svg_text = " ".join(
        node.text or "" for node in root.iter() if node.tag.endswith("text")
    )
    image_count = sum(1 for node in root.iter() if node.tag.endswith("image"))
    text_count = sum(1 for node in root.iter() if node.tag.endswith("text"))
    require(image_count == 1, f"Expected one GPT Image background, found {image_count}")
    require(text_count >= 70, f"Too few editable SVG text elements: {text_count}")

    pdf_text = command_output("pdftotext", str(pdf_path), "-")
    combined_text = f"{svg_text}\n{pdf_text}"
    for label in REQUIRED_TEXT:
        require(label in combined_text, f"Required label is missing: {label}")
    for term in FORBIDDEN_VISUAL_TERMS:
        require(
            re.search(rf"\b{re.escape(term)}\b", combined_text, re.I) is None,
            f"Forbidden visual term remains: {term}",
        )

    pdf_info = command_output("pdfinfo", str(pdf_path))
    require("Pages:           1" in pdf_info, "PDF must contain one page")
    page_match = re.search(
        r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info
    )
    require(page_match is not None, "Could not read PDF page size")
    width_pt, height_pt = map(float, page_match.groups())
    require(abs(width_pt - 518.74) < 0.1, "PDF width is not 183 mm")
    require(abs(height_pt - 221.102) < 0.1, "PDF height is not 78 mm")

    font_info = command_output("pdffonts", str(pdf_path))
    require("Arial-BoldMT" in font_info, "Arial Bold is not embedded")
    require("ArialMT" in font_info, "Arial Regular is not embedded")
    require("Type 3" not in font_info, "Type 3 fonts are not allowed")

    with Image.open(png_path) as png:
        require(png.size == (3242, 1381), f"Unexpected PNG size: {png.size}")
        dpi = png.info.get("dpi")
        require(dpi is not None, "PNG DPI metadata missing")
        require(
            all(abs(value - 450.0) < 0.1 for value in dpi),
            f"PNG is not 450 dpi: {dpi}",
        )
        print("PNG size/DPI:", png.size, dpi)

    print("PDF INFO")
    print(pdf_info.strip())
    print("PDF FONTS")
    print(font_info.strip())
    print(
        "EDITABILITY ASSERTIONS PASS",
        {"svg_image_elements": image_count, "svg_text_elements": text_count},
    )
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v2"),
    )
    parser.add_argument(
        "--component-dir",
        type=Path,
        default=Path("figures/ga_gpt2_components_v2"),
    )
    args = parser.parse_args()

    check_catalog(args.component_dir)
    exports = check_exports(args.stem, args.component_dir)
    provenance_files = [
        *exports,
        args.component_dir
        / "mcp_generation"
        / "gpt_image_20260726_163535_6a65b268_2.png",
        args.component_dir / "catalog_system_composition_20260726.tsv",
        args.component_dir / "catalog_study_biology_classification_20260726.tsv",
        args.component_dir / "featured_project_cohorts_20260726.tsv",
        Path("figures/scripts/render_ga_gpt2_hybrid.py"),
        Path("figures/scripts/render_ga_gpt2_real_components.py"),
        Path(__file__),
    ]
    print("SHA256")
    for path in provenance_files:
        require(path.is_file(), f"Missing provenance file: {path}")
        print(sha256(path), path)
    print("ALL QA CHECKS PASS")


if __name__ == "__main__":
    main()
