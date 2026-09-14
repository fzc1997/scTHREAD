#!/usr/bin/env python3
"""Machine QA for the simplified NAR graphical abstract."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

from PIL import Image
from PyPDF2 import PdfReader


WIDTH_MM = 150.0
HEIGHT_MM = 60.0


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stem", type=Path)
    args = parser.parse_args()
    stem = args.stem.resolve()
    paths = {suffix: stem.with_suffix(suffix) for suffix in (".pdf", ".svg", ".png", ".tiff")}
    require(all(path.is_file() and path.stat().st_size > 0 for path in paths.values()), "Missing export")

    pdf = PdfReader(str(paths[".pdf"]))
    require(len(pdf.pages) == 1, "PDF must be one page")
    box = pdf.pages[0].mediabox
    width_mm = float(box.width) / 72 * 25.4
    height_mm = float(box.height) / 72 * 25.4
    require(abs(width_mm - WIDTH_MM) < 0.1, "PDF width mismatch")
    require(abs(height_mm - HEIGHT_MM) < 0.1, "PDF height mismatch")
    require(abs(width_mm / height_mm - 2.5) < 0.005, "Canvas is not 5:2")

    svg = paths[".svg"].read_text()
    require(svg.count("<text") >= 12, "Expected editable SVG text")
    require("Arial" in svg, "Arial is not declared in SVG")
    sizes = [
        float(value)
        for value in re.findall(r"(?:font-size:|font: (?:[0-9]+ )?)([0-9.]+)px", svg)
    ]
    require(sizes and min(sizes) >= 12, f"SVG font below 12 pt/px: {min(sizes, default=-1)}")
    for phrase in (
        "scTHREAD",
        "Multi-study data",
        "Evidence layers",
        "Open retrieval",
        "https://scthread.ai4sc.ac.cn",
    ):
        require(phrase in svg, f"Missing phrase: {phrase}")

    with Image.open(paths[".png"]) as image:
        require(image.mode == "RGBA", "PNG must preserve transparency")
        require(image.size == (3543, 1417), f"Unexpected 600-dpi PNG size: {image.size}")
        require(image.getchannel("A").getextrema()[0] == 0, "PNG background is not transparent")

    print(f"NAR GRAPHICAL ABSTRACT QA PASS\tsha256={sha256(paths['.pdf'])}")


if __name__ == "__main__":
    main()
