#!/usr/bin/env python3
"""Add a registry-backed tissue-icon strip to the scTHREAD NAR graphical abstract.

The base PDF is the collaborator's Illustrator export
(``NAR_GraphicalAbstract_20260804.pdf``).  It is never edited in place: the new
icon strip is drawn as a self-contained vector overlay and merged on top, so
every existing object in the base file is preserved byte-for-byte.

The overlay replaces only the seven generic pictograms in the ``Biological
coverage`` panel with sixteen organ icons that correspond to tissues actually
present in the frozen 34-study / 453-run / 923,389-cell release, colour-coded by
the same six biological-system classes used by the KPI card legend directly
below the strip.

Icon sources (public repositories, permissive licences):

* healthicons -- https://github.com/resolvetosavelives/healthicons (MIT)
* Tabler Icons -- https://github.com/tabler/tabler-icons (MIT)

Both are MIT, so the published figure needs no icon credit line.

Geometry was measured from the base PDF itself (1200-dpi ink profiling plus
``pdftotext -bbox``), not assumed, so the white-out band sits strictly between
the subtitle ink and the KPI card border.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import cairosvg
from PIL import ImageFont
from PyPDF2 import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------------------------
# Base-page geometry, all in PostScript points, y measured from the page top.
# ---------------------------------------------------------------------------
PAGE_W = 518.74
PAGE_H = 221.102

# Measured ink extents in the base PDF (1200 dpi):
#   title  "Biological coverage"      y 16.000 - 22.440
#   subtitle, x 24-37 (no descender)  y 23.460 - 26.340
#   subtitle, x 37-115 (descenders)   y 23.460 - 27.120
#   old icon strip                    y 27.120 - 63.240
#   KPI card top border               y 69.180
# The strip is therefore cleared with two rectangles whose top edges fall in the
# measured blank gaps, and whose bottom edge stops 0.68 pt above the card.
CLEAR_LEFT = dict(x0=9.5, x1=37.0, y0=26.70, y1=68.50)
CLEAR_RIGHT = dict(x0=37.0, x1=137.5, y0=27.30, y1=68.50)

GRID_X0, GRID_X1 = 13.0, 133.5
N_COLS = 9
ICON_BOX = 12.0
ROW_ICON_CY = (34.0, 54.0)
ROW_LABEL_BASELINE = (43.6, 63.6)
LABEL_SIZE_MAX = 3.7
LABEL_SIZE_MIN = 3.15
LABEL_COLOR = "#3c4043"

# Six mutually exclusive release classes; hex values sampled from the KPI card
# legend dots in the base PDF so the icons key to the composition bar below.
SYSTEM_COLOR = {
    "Blood/immune": "#0b7c86",
    "Neural/sensory": "#5b4b9a",
    "Cancer": "#d55e00",
    "Endocrine": "#2f6fb2",
    "Heart/vascular": "#7a9e3a",
    "Development/embryo": "#d88b1b",
    "Species": "#443677",
}

# (label, icon key, system).  Row-major, 9 per row.
# Every organ entry is backed by at least one study in the frozen release; see
# tables/ga_tissue_icon_map.tsv for the study-level evidence.
"""
The two species markers are kept adjacent rather than one per row: the rows are
mixed (Marrow/Blood/Lymph/Heart are human, Retina/Cochlea/Spleen are mouse), so a
species icon at the head of each row would read as a row label and mislead.
"""
STRIP = [
    ("Human", "body", "Species"),
    ("Mouse", "animal-rat", "Species"),
    ("Brain", "neurology", "Neural/sensory"),
    ("Retina", "eye", "Neural/sensory"),
    ("Cochlea", "ear", "Neural/sensory"),
    ("Marrow", "bone", "Blood/immune"),
    ("Blood", "blood-drop", "Blood/immune"),
    ("Spleen", "spleen", "Blood/immune"),
    ("Lymph", "lymph-nodes", "Blood/immune"),
    ("Heart", "heart-organ", "Heart/vascular"),
    ("Lung", "lungs", "Cancer"),
    ("Kidney", "kidneys", "Cancer"),
    ("Prostate", "prostate", "Cancer"),
    ("Cell line", "cell-nuclei", "Cancer"),
    ("Islet", "pancreas", "Endocrine"),
    ("Testis", "testicles", "Endocrine"),
    ("Ovary", "female-reproductive_system", "Endocrine"),
    ("Embryo", "mouse_embryo", "Development/embryo"),
]

# Some icon sets draw a bone on the diagonal.  Rotating it flat matches the
# upright neighbours; the extra shrink keeps the rotated shape inside its cell.
ICON_ROTATE = {"bone": 45}
ICON_SCALE = {"bone": 0.88}

ARIAL = Path("/gpfs/home/fuzc/.local/share/fonts/Arial.ttf")


# ---------------------------------------------------------------------------
def load_icon(icon_dir: Path, key: str) -> tuple[str, float, float]:
    """Return (inner SVG markup, viewBox width, viewBox height) for one icon."""
    for sub in ("healthicons", "tabler", "anchor"):
        path = icon_dir / sub / f"{key}.svg"
        if path.is_file():
            break
    else:
        raise FileNotFoundError(f"icon not found: {key}")

    text = path.read_text(encoding="utf-8")
    match = re.search(r'viewBox="([\d.\-\s]+)"', text)
    if not match:
        raise ValueError(f"{path} has no viewBox")
    _, _, vb_w, vb_h = (float(v) for v in match.group(1).split())

    inner = re.sub(r"^.*?<svg[^>]*>", "", text, flags=re.S)
    inner = re.sub(r"</svg>\s*$", "", inner, flags=re.S)
    inner = re.sub(r"<!--.*?-->", "", inner, flags=re.S)
    # healthicons and Tabler both use currentColor; strip it (and any fill="none"
    # on the root) so the wrapping group's fill applies uniformly.
    inner = inner.replace('fill="currentColor"', "").replace('fill="none"', "")
    return inner.strip(), vb_w, vb_h


def fit_label_size(text: str, max_width: float) -> float:
    """Largest Arial size (<= LABEL_SIZE_MAX) whose advance width fits."""
    size = LABEL_SIZE_MAX
    while size > LABEL_SIZE_MIN:
        # Measure at 100x then scale down: PIL truetype sizes are integers.
        font = ImageFont.truetype(str(ARIAL), int(round(size * 100)))
        width = font.getlength(text) / 100.0
        if width <= max_width:
            return size
        size -= 0.05
    return LABEL_SIZE_MIN


def build_overlay_svg(icon_dir: Path) -> str:
    pitch = (GRID_X1 - GRID_X0) / N_COLS
    parts = [
        # xlink must be declared here: the icon loader strips each source <svg>
        # root, and art converted from Illustrator references its masks and clip
        # paths through xlink:href.
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{PAGE_W}pt" height="{PAGE_H}pt" '
        f'viewBox="0 0 {PAGE_W} {PAGE_H}">'
    ]

    for rect in (CLEAR_LEFT, CLEAR_RIGHT):
        parts.append(
            '<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#ffffff"/>'.format(
                x0=rect["x0"],
                y0=rect["y0"],
                w=rect["x1"] - rect["x0"],
                h=rect["y1"] - rect["y0"],
            )
        )

    parts.append('<g id="scthread-tissue-icons">')
    for index, (label, key, system) in enumerate(STRIP):
        row, col = divmod(index, N_COLS)
        cx = GRID_X0 + pitch * (col + 0.5)
        cy = ROW_ICON_CY[row]
        color = SYSTEM_COLOR[system]

        inner, vb_w, vb_h = load_icon(icon_dir, key)
        scale = ICON_BOX / max(vb_w, vb_h) * ICON_SCALE.get(key, 1.0)
        tx = cx - vb_w * scale / 2.0
        ty = cy - vb_h * scale / 2.0
        spin = ICON_ROTATE.get(key)
        rot = f" rotate({spin} {vb_w / 2:.4f} {vb_h / 2:.4f})" if spin else ""
        parts.append(
            f'<g transform="translate({tx:.4f} {ty:.4f}) scale({scale:.6f}){rot}" '
            f'fill="{color}" color="{color}">{inner}</g>'
        )

        size = fit_label_size(label, pitch - 1.0)
        parts.append(
            f'<text x="{cx:.4f}" y="{ROW_LABEL_BASELINE[row]:.4f}" '
            f'font-family="Arial" font-size="{size:.2f}" fill="{LABEL_COLOR}" '
            f'text-anchor="middle">{label}</text>'
        )
    parts.append("</g></svg>")
    return "\n".join(parts)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--icon-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--overlay-svg", type=Path)
    args = parser.parse_args()

    svg = build_overlay_svg(args.icon_dir)
    overlay_svg = args.overlay_svg or args.out.with_suffix(".overlay.svg")
    overlay_svg.write_text(svg, encoding="utf-8")

    # cairosvg treats output_width/output_height as CSS pixels (96 dpi), which
    # silently rescales the page by 0.75.  The SVG already declares its size in
    # points, so let cairosvg honour that and assert the result afterwards.
    overlay_pdf = overlay_svg.with_suffix(".pdf")
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(overlay_pdf))

    base = PdfReader(str(args.base))
    over = PdfReader(str(overlay_pdf))
    if len(base.pages) != 1:
        raise SystemExit(f"expected a single-page base, got {len(base.pages)}")

    over_box = over.pages[0].mediabox
    if (
        abs(float(over_box.width) - PAGE_W) > 0.05
        or abs(float(over_box.height) - PAGE_H) > 0.05
    ):
        raise SystemExit(
            f"overlay page {float(over_box.width)}x{float(over_box.height)} "
            f"is not {PAGE_W}x{PAGE_H}; merging would misplace every icon"
        )

    page = base.pages[0]
    box = page.mediabox
    if abs(float(box.width) - PAGE_W) > 0.05 or abs(float(box.height) - PAGE_H) > 0.05:
        raise SystemExit(
            f"base media box {float(box.width)}x{float(box.height)} "
            f"does not match the measured geometry {PAGE_W}x{PAGE_H}"
        )
    page.merge_page(over.pages[0])

    # The base was exported with "preserve Illustrator editing capabilities":
    # /PieceInfo holds a private copy of the original artwork that Illustrator
    # reads *instead of* the PDF content stream.  Left in place it would show the
    # old pictograms in Illustrator while every other viewer showed the new
    # strip, so drop it (and the stale page thumbnail) to keep one truth.
    for key in ("/PieceInfo", "/LastModified", "/Thumb"):
        if key in page:
            del page[key]

    writer = PdfWriter()
    writer.add_page(page)
    with args.out.open("wb") as handle:
        writer.write(handle)

    manifest = {
        "base_pdf": str(args.base),
        "base_sha256": sha256(args.base),
        "output_pdf": str(args.out),
        "output_sha256": sha256(args.out),
        "overlay_svg_sha256": sha256(overlay_svg),
        "n_icons": len(STRIP),
        "icons": [
            {"label": lab, "icon": key, "system": sys_} for lab, key, sys_ in STRIP
        ],
        "clear_rects": [CLEAR_LEFT, CLEAR_RIGHT],
    }
    args.out.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"wrote {args.out}")
    print(f"      {overlay_svg}")
    print(f"      {args.out.with_suffix('.manifest.json')}")


if __name__ == "__main__":
    main()
