#!/usr/bin/env python3
"""Apply every scTHREAD graphical-abstract patch to the collaborator's PDF in one pass.

Two independent overlays are merged onto the untouched Illustrator export:

1. the registry-backed tissue-icon strip (``add_tissue_icons_to_ga``);
2. the corrected gene-annotation card (``patch_ga_ase_card``).

The two overlays are concatenated into a **single** SVG, rasterised once and
merged once.  Two successive ``merge_page`` calls on the same page render
correctly but leave a dangling cross-reference entry -- poppler then reports
``Internal Error: xref num ... not found but needed, try to reconstruct`` -- which
is not something to hand a journal submission system.  Verified: either overlay
merged alone is clean, both merged separately is not, both merged as one is.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import cairosvg
from PyPDF2 import PdfReader, PdfWriter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAGE_W, PAGE_H = 518.74, 221.102


def load(module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, HERE / f"{module_name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def combine(first: str, second: str) -> str:
    """Splice two full-page overlay SVGs into one document, preserving order."""
    body = re.sub(r"^.*?<svg[^>]*>", "", second, flags=re.S)
    body = re.sub(r"</svg>\s*$", "", body, flags=re.S)
    if "</svg>" not in first:
        raise SystemExit("first overlay has no closing svg tag")
    return first.rsplit("</svg>", 1)[0] + body + "</svg>"


def to_pdf(svg: str, path: Path) -> PdfReader:
    path.write_text(svg, encoding="utf-8")
    pdf_path = path.with_suffix(".pdf")
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(pdf_path))
    reader = PdfReader(str(pdf_path))
    box = reader.pages[0].mediabox
    if abs(float(box.width) - PAGE_W) > 0.05 or abs(float(box.height) - PAGE_H) > 0.05:
        raise SystemExit(f"{path.name}: overlay page is {float(box.width)}x{float(box.height)}")
    return reader


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            d.update(block)
    return d.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path, help="collaborator's untouched export")
    ap.add_argument("--icon-dir", type=Path, default=ROOT / "figures/assets/tissue_icons")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    icons = load("add_tissue_icons_to_ga")
    card = load("patch_ga_ase_card")

    counts = card.load_counts()
    icon_svg = icons.build_overlay_svg(args.icon_dir)
    card_svg = card.build_overlay(counts)

    overlay_pdf = to_pdf(combine(icon_svg, card_svg), args.out.with_suffix(".overlay.svg"))

    base = PdfReader(str(args.base))
    if len(base.pages) != 1:
        raise SystemExit(f"expected a single-page base, got {len(base.pages)}")
    page = base.pages[0]
    box = page.mediabox
    if abs(float(box.width) - PAGE_W) > 0.05 or abs(float(box.height) - PAGE_H) > 0.05:
        raise SystemExit(f"base media box {float(box.width)}x{float(box.height)} is unexpected")

    page.merge_page(overlay_pdf.pages[0])
    # Illustrator's private artwork copy would otherwise be shown instead of the
    # merged content when the file is opened in Illustrator.
    for key in ("/PieceInfo", "/LastModified", "/Thumb"):
        if key in page:
            del page[key]

    writer = PdfWriter()
    writer.add_page(page)
    with args.out.open("wb") as fh:
        writer.write(fh)

    args.out.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "base_pdf": str(args.base),
                "base_sha256": sha256(args.base),
                "output_pdf": str(args.out),
                "output_sha256": sha256(args.out),
                "patches": ["tissue_icon_strip", "gene_annotation_card"],
                "icons": [
                    {"label": lab, "icon": key, "system": sys_}
                    for lab, key, sys_ in icons.STRIP
                ],
                "card_counts": {
                    a: {"significant": s, "tested": t, "fraction": round(s / t, 6)}
                    for a, (s, t) in counts.items()
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    for axis, (s, t) in counts.items():
        print(f"  {axis}: {s:,} / {t:,}  ({100 * s / t:.1f}%)")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
