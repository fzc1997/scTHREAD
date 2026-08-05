#!/usr/bin/env python3
"""Replace the ASE row of the graphical abstract's gene-annotation card.

Why this patch exists
---------------------
The card reported ``ASE 0/6,930 genes``.  That zero is the *cell-type-differential*
ASE test -- does allelic imbalance differ between cell types -- which is
underpowered in the two eligible marrow studies and, as the user put it, carries
no biological meaning here.  The biologically meaningful quantity, allelic
imbalance itself, is not zero: ``tables/ase_allelic_imbalance_feature_20260802.tsv``
finds **9,880 of 37,415 genes** significantly imbalanced (q < 0.05,
|deviation| >= 0.10) across 93 runs and 18 studies.  Showing the zero as "the ASE
result" misrepresents the resource.

Because ASE now reports a different question from DIU/APA, the card's title and
subtitle can no longer say "Cell-type-associated genes ... criterion" for all
three rows, so the whole card interior is redrawn rather than just the number.

Encoding
--------
Bar length is proportional to the **significant fraction**, not the raw count.
The three axes test very different numbers of genes (8,092 / 10,531 / 37,415), so
raw-count bars would imply ASE is five times richer when it merely tested 4.6x
more genes.  Numerator and denominator are printed on every row.

The base PDF is never edited in place: this draws a vector overlay and merges it,
exactly like ``add_tissue_icons_to_ga.py``.  All geometry below was measured from
the base PDF (300-dpi colour sampling plus ``pdftotext -bbox``).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import cairosvg
from PIL import ImageFont
from PyPDF2 import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[2]
PAGE_W = 518.74
PAGE_H = 221.102

# Clearance verified from the base render: the nearest ink left of the card ends
# at x 254.88 and the next starts at 397.64; above, the inbound arrow ends at
# y 168.07 and the next element starts at y 208.60.
CLEAR = dict(x0=255.6, y0=168.6, x1=392.2, y1=207.8)

# Card shell, on the stroke centre-line.
CARD = dict(x0=257.68, y0=170.88, x1=390.32, y1=204.12, rx=2.0)
CARD_FILL = "#fff8ed"
CARD_STROKE = "#efb96b"
CARD_STROKE_W = 0.55

TITLE = "Precomputed gene-level annotations"
SUBTITLE = "significant / tested genes · FDR or q < 0.05 + effect size"
FOOTNOTE = (
    "DIU · APA cell-type-differential (effect ≥ 0.20) · "
    "ASE allelic imbalance (|dev| ≥ 0.10)"
)

TITLE_SIZE, TITLE_BASE, TITLE_COLOR = 6.41, 176.23, "#202124"
SUB_SIZE, SUB_BASE, SUB_COLOR = 4.48, 180.52, "#6d6e73"
LABEL_SIZE, LABEL_CX, LABEL_COLOR = 5.13, 263.0, "#000000"
VALUE_SIZE, VALUE_COLOR = 4.28, "#202124"
FOOT_SIZE, FOOT_BASE, FOOT_COLOR = 3.45, 203.0, "#6d6e73"

BAR_X0, BAR_H, BAR_MAX = 269.0, 2.40, 75.0
VALUE_GAP = 3.2

# (axis, baseline y, bar colour).  DIU/APA colours sampled from the base PDF;
# ASE takes the coral of the "Allele-aware evidence" card so the row keys to the
# evidence layer it comes from.
ROWS = [
    ("DIU", 187.68, "#5b4b9a"),
    ("APA", 193.62, "#d88b1b"),
    ("ASE", 199.61, "#c8574e"),
]

ARIAL = Path("/gpfs/home/fuzc/.local/share/fonts/Arial.ttf")
ARIAL_BOLD = Path("/gpfs/home/fuzc/.local/share/fonts/Arial-Bold.ttf")

THREE_AXIS = ROOT / "tables/Table_S_three_axis_summary.tsv"
ASE_FEATURE = ROOT / "tables/ase_allelic_imbalance_feature_20260802.tsv"
ASE_SUMMARY = ROOT / "tables/ase_allelic_imbalance_feature_20260802.summary.json"


def read_tsv(path: Path) -> list[dict]:
    raw = path.read_text(encoding="utf-8", errors="replace").replace("\r", "")
    return [r for r in csv.DictReader(raw.split("\n"), delimiter="\t") if any(r.values())]


def load_counts() -> dict[str, tuple[int, int]]:
    """Return {axis: (significant, tested)} recomputed from the source tables."""
    counts: dict[str, tuple[int, int]] = {}

    for row in read_tsv(THREE_AXIS):
        axis = (row.get("axis") or "").strip()
        if axis in {"DIU", "APA"}:
            counts[axis] = (
                int(row["n_sig_fdr05_effect_gate"]),
                int(row["n_genes_tested"]),
            )

    # ASE is recomputed from the per-gene table rather than trusting the summary,
    # and the recomputation is asserted against the stored flag.
    rows = read_tsv(ASE_FEATURE)
    tested = len(rows)
    flagged = sum(1 for r in rows if str(r["sig_imbalanced"]).lower() in {"true", "1"})
    recomputed = sum(
        1
        for r in rows
        if r["qval"] not in ("", "NA")
        and float(r["qval"]) < 0.05
        and r["pooled_deviation"] not in ("", "NA")
        and abs(float(r["pooled_deviation"])) >= 0.10
    )
    if flagged != recomputed:
        raise SystemExit(
            f"ASE flag {flagged} disagrees with the recomputed criterion {recomputed}"
        )
    summary = json.loads(ASE_SUMMARY.read_text(encoding="utf-8"))
    if summary["sig_imbalanced_genes_q05_eff10"] != flagged:
        raise SystemExit("ASE summary JSON disagrees with the per-gene table")
    counts["ASE"] = (flagged, tested)

    missing = {"DIU", "APA", "ASE"} - set(counts)
    if missing:
        raise SystemExit(f"missing axes: {sorted(missing)}")
    return counts


def text_width(text: str, size: float, bold: bool = False) -> float:
    font = ImageFont.truetype(str(ARIAL_BOLD if bold else ARIAL), int(round(size * 100)))
    return font.getlength(text) / 100.0


def fit_size(text: str, max_width: float, start: float, floor: float) -> float:
    size = start
    while size > floor and text_width(text, size) > max_width:
        size -= 0.02
    return size


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def build_overlay(counts: dict[str, tuple[int, int]]) -> str:
    fractions = {a: s / t for a, (s, t) in counts.items()}
    top = max(fractions.values())
    cx = (CARD["x0"] + CARD["x1"]) / 2.0

    p = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{PAGE_W}pt" height="{PAGE_H}pt" viewBox="0 0 {PAGE_W} {PAGE_H}">',
        '<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#ffffff"/>'.format(
            x0=CLEAR["x0"], y0=CLEAR["y0"],
            w=CLEAR["x1"] - CLEAR["x0"], h=CLEAR["y1"] - CLEAR["y0"],
        ),
        '<g id="scthread-annotation-card">',
        f'<rect x="{CARD["x0"]}" y="{CARD["y0"]}" '
        f'width="{CARD["x1"] - CARD["x0"]:.4f}" height="{CARD["y1"] - CARD["y0"]:.4f}" '
        f'rx="{CARD["rx"]}" fill="{CARD_FILL}" stroke="{CARD_STROKE}" '
        f'stroke-width="{CARD_STROKE_W}"/>',
        f'<text x="{cx:.3f}" y="{TITLE_BASE}" font-family="Arial" font-weight="bold" '
        f'font-size="{TITLE_SIZE}" fill="{TITLE_COLOR}" text-anchor="middle">{esc(TITLE)}</text>',
        f'<text x="{cx:.3f}" y="{SUB_BASE}" font-family="Arial" font-weight="bold" '
        f'font-size="{SUB_SIZE}" fill="{SUB_COLOR}" text-anchor="middle">{esc(SUBTITLE)}</text>',
    ]

    for axis, baseline, colour in ROWS:
        sig, tested = counts[axis]
        length = BAR_MAX * fractions[axis] / top
        bar_y = baseline - 3.05
        p.append(
            f'<text x="{LABEL_CX}" y="{baseline}" font-family="Arial" font-weight="bold" '
            f'font-size="{LABEL_SIZE}" fill="{LABEL_COLOR}" text-anchor="middle">{axis}</text>'
        )
        p.append(
            f'<rect x="{BAR_X0}" y="{bar_y:.3f}" width="{length:.3f}" height="{BAR_H}" '
            f'fill="{colour}"/>'
        )
        value = f"{sig:,} / {tested:,} genes"
        p.append(
            f'<text x="{BAR_X0 + length + VALUE_GAP:.3f}" y="{baseline - 0.55:.3f}" '
            f'font-family="Arial" font-size="{VALUE_SIZE}" fill="{VALUE_COLOR}" '
            f'text-anchor="start">{esc(value)}</text>'
        )

    widest = max(
        BAR_X0 + BAR_MAX * fractions[a] / top + VALUE_GAP
        + text_width(f"{counts[a][0]:,} / {counts[a][1]:,} genes", VALUE_SIZE)
        for a, _, _ in ROWS
    )
    if widest > CARD["x1"] - 2.0:
        raise SystemExit(
            f"row text reaches x={widest:.2f}, past the card inner edge "
            f"{CARD['x1'] - 2.0:.2f}; reduce BAR_MAX or VALUE_SIZE"
        )

    inner = CARD["x1"] - CARD["x0"] - 8.0
    foot_size = fit_size(FOOTNOTE, inner, FOOT_SIZE, 2.9)
    p.append(
        f'<text x="{cx:.3f}" y="{FOOT_BASE}" font-family="Arial" font-size="{foot_size:.2f}" '
        f'fill="{FOOT_COLOR}" text-anchor="middle">{esc(FOOTNOTE)}</text>'
    )
    p.append("</g></svg>")
    return "\n".join(p)


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            d.update(block)
    return d.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    counts = load_counts()
    svg = build_overlay(counts)
    overlay_svg = args.out.with_suffix(".card.svg")
    overlay_svg.write_text(svg, encoding="utf-8")

    overlay_pdf = overlay_svg.with_suffix(".pdf")
    cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=str(overlay_pdf))

    over = PdfReader(str(overlay_pdf))
    box = over.pages[0].mediabox
    if abs(float(box.width) - PAGE_W) > 0.05 or abs(float(box.height) - PAGE_H) > 0.05:
        raise SystemExit(f"overlay page {float(box.width)}x{float(box.height)} is wrong")

    base = PdfReader(str(args.base))
    page = base.pages[0]
    page.merge_page(over.pages[0])
    for key in ("/PieceInfo", "/LastModified", "/Thumb"):
        if key in page:
            del page[key]

    writer = PdfWriter()
    writer.add_page(page)
    with args.out.open("wb") as fh:
        writer.write(fh)

    args.out.with_suffix(".card.json").write_text(
        json.dumps(
            {
                "base_pdf": str(args.base),
                "base_sha256": sha256(args.base),
                "output_sha256": sha256(args.out),
                "counts": {a: {"significant": s, "tested": t} for a, (s, t) in counts.items()},
                "fractions": {a: round(s / t, 6) for a, (s, t) in counts.items()},
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
