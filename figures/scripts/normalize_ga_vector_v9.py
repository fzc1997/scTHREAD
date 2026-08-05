#!/usr/bin/env python3
"""Remove volatile SVG/PDF metadata from graphical abstract v9 exports."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from lxml import etree
from PyPDF2 import PdfReader, PdfWriter
from PyPDF2.generic import ArrayObject, ByteStringObject


FIXED_PDF_DATE = "D:20260728220000+08'00'"
FIXED_PDF_ID = bytes.fromhex("64f0d661c23c41bc9c82e69b47fcce7c")


def local_name(tag: object) -> str:
    if not isinstance(tag, str):
        return ""
    return tag.rsplit("}", 1)[-1]


def normalize_svg(path: Path) -> None:
    parser = etree.XMLParser(
        resolve_entities=False,
        no_network=True,
        remove_blank_text=False,
        huge_tree=False,
    )
    tree = etree.parse(str(path), parser)
    root = tree.getroot()
    for node in list(root.iter()):
        if local_name(node.tag) != "metadata":
            continue
        parent = node.getparent()
        if parent is not None:
            parent.remove(node)
    temporary = path.with_name(path.name + ".normalized.tmp")
    tree.write(
        str(temporary),
        encoding="UTF-8",
        xml_declaration=True,
        pretty_print=True,
    )
    os.replace(temporary, path)


def normalize_pdf(path: Path) -> None:
    reader = PdfReader(str(path))
    writer = PdfWriter()
    for page in reader.pages:
        writer.add_page(page)
    writer.add_metadata(
        {
            "/Title": "scTHREAD graphical abstract v9",
            "/Creator": "Python/Matplotlib/CairoSVG",
            "/Producer": "scTHREAD deterministic vector export",
            "/CreationDate": FIXED_PDF_DATE,
            "/ModDate": FIXED_PDF_DATE,
        }
    )
    fixed_id = ByteStringObject(FIXED_PDF_ID)
    writer._ID = ArrayObject([fixed_id, fixed_id])
    temporary = path.with_name(path.name + ".normalized.tmp")
    with temporary.open("wb") as handle:
        writer.write(handle)
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v9"),
    )
    args = parser.parse_args()
    svg_path = args.stem.with_suffix(".svg")
    pdf_path = args.stem.with_suffix(".pdf")
    if not svg_path.is_file() or not pdf_path.is_file():
        raise FileNotFoundError(f"Missing export for stem: {args.stem}")
    normalize_svg(svg_path)
    normalize_pdf(pdf_path)
    print(f"{svg_path}\t{svg_path.stat().st_size} bytes")
    print(f"{pdf_path}\t{pdf_path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
