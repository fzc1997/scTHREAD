#!/usr/bin/env python3
"""Extract a sub-drawing from an Illustrator PDF/.ai page as standalone vector art.

The ANCHOR Extended Data Fig. 3 panel draws six overlapping mouse embryos inside
one dashed container, so a rectangular crop cannot isolate a single embryo.  This
tool instead works on the content stream: Illustrator emits every path as a
self-contained ``q <cm> ... <paint> Q`` block, so each block's bounding box can be
computed in page space and blocks outside a target rectangle can be dropped while
every colour/graphics-state operator between them is kept in order.

Two modes:

* ``--survey``  print the bounding box of every block that intersects the region,
  so the caller can pick the exact sub-drawing.
* ``--extract`` write a one-page PDF containing only the blocks whose bounding box
  is *inside* the region, with the media box tightened to the retained ink.

Nothing is redrawn or traced: the retained path operators are byte-identical to
the source, so the extracted art is the original vector artwork.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from PyPDF2 import PdfReader, PdfWriter
from PyPDF2.generic import DecodedStreamObject, NameObject, RectangleObject

NUM = r"-?\d*\.?\d+"
TOKEN = re.compile(rb"(?P<num>-?\d*\.?\d+)|(?P<name>/[^\s/\[\]<>(]+)|(?P<op>[A-Za-z'\"*]+)")

# Operators whose operands are coordinate pairs we must transform.
PATH_OPS = {"m": 1, "l": 1, "c": 3, "v": 2, "y": 2}
PAINT_OPS = {"f", "F", "f*", "B", "B*", "b", "b*", "S", "s", "n"}


def split_blocks(data: bytes) -> list[tuple[str, bytes]]:
    """Split a content stream into top-level ('q' block | 'other') chunks."""
    blocks: list[tuple[str, bytes]] = []
    depth = 0
    start = 0
    pos = 0
    plain_start = 0
    for match in re.finditer(rb"(?<![A-Za-z])(q|Q)(?![A-Za-z])", data):
        token = match.group(1)
        if token == b"q":
            if depth == 0:
                if match.start() > plain_start:
                    blocks.append(("other", data[plain_start : match.start()]))
                start = match.start()
            depth += 1
        else:
            depth -= 1
            if depth == 0:
                blocks.append(("q", data[start : match.end()]))
                plain_start = match.end()
            elif depth < 0:  # unbalanced Q; treat as plain text
                depth = 0
                plain_start = min(plain_start, match.start())
    if plain_start < len(data):
        blocks.append(("other", data[plain_start:]))
    return blocks


def block_bbox(block: bytes) -> tuple[float, float, float, float] | None:
    """Bounding box of a q-block in page space, or None if it paints nothing."""
    ctm = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    stack: list[tuple] = []
    operands: list = []
    xs: list[float] = []
    ys: list[float] = []
    painted = False

    for match in TOKEN.finditer(block):
        if match.group("num") is not None:
            operands.append(float(match.group("num")))
            continue
        if match.group("name") is not None:
            operands.append(match.group("name"))
            continue
        op = match.group("op").decode("latin-1")

        if op == "q":
            stack.append(ctm)
        elif op == "Q":
            if stack:
                ctm = stack.pop()
        elif op == "cm" and len(operands) >= 6:
            a, b, c, d, e, f = (float(v) for v in operands[-6:])
            a0, b0, c0, d0, e0, f0 = ctm
            ctm = (
                a * a0 + b * c0,
                a * b0 + b * d0,
                c * a0 + d * c0,
                c * b0 + d * d0,
                e * a0 + f * c0 + e0,
                e * b0 + f * d0 + f0,
            )
        elif op in PATH_OPS:
            pairs = [v for v in operands if isinstance(v, float)]
            for i in range(0, len(pairs) - 1, 2):
                x, y = pairs[i], pairs[i + 1]
                a, b, c, d, e, f = ctm
                xs.append(a * x + c * y + e)
                ys.append(b * x + d * y + f)
        elif op == "re" and len(operands) >= 4:
            x, y, w, h = (float(v) for v in operands[-4:])
            a, b, c, d, e, f = ctm
            for px, py in ((x, y), (x + w, y), (x, y + h), (x + w, y + h)):
                xs.append(a * px + c * py + e)
                ys.append(b * px + d * py + f)
        elif op in PAINT_OPS and op != "n":
            painted = True

        operands = []

    if not xs or not painted:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", required=True, type=Path)
    parser.add_argument("--region", required=True, nargs=4, type=float,
                        metavar=("X0", "Y0", "X1", "Y1"),
                        help="page-space rectangle, PDF coordinates (origin bottom-left)")
    parser.add_argument("--survey", action="store_true")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--pad", type=float, default=0.5)
    args = parser.parse_args()

    x0, y0, x1, y1 = args.region
    reader = PdfReader(str(args.src))
    page = reader.pages[0]
    data = page.get_contents().get_data()
    blocks = split_blocks(data)

    kept: list[bytes] = []
    boxes = []
    for kind, chunk in blocks:
        if kind == "other":
            kept.append(chunk)
            continue
        box = block_bbox(chunk)
        if box is None:
            continue
        bx0, by0, bx1, by1 = box
        intersects = not (bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1)
        inside = bx0 >= x0 and by0 >= y0 and bx1 <= x1 and by1 <= y1
        if args.survey and intersects:
            boxes.append((bx0, by0, bx1, by1, inside, len(chunk)))
        if inside:
            kept.append(chunk)

    if args.survey:
        boxes.sort(key=lambda b: (-b[3], b[0]))
        print(f"{len(boxes)} blocks intersect the region")
        print("  x0       y0       x1       y1       inside  bytes")
        for bx0, by0, bx1, by1, inside, size in boxes:
            print(f"  {bx0:8.2f} {by0:8.2f} {bx1:8.2f} {by1:8.2f} {str(inside):>6}  {size}")
        return

    if args.out is None:
        raise SystemExit("--out is required unless --survey")

    ink = [block_bbox(c) for k, c in blocks if k == "q"]
    ink = [b for b in ink if b and b[0] >= x0 and b[1] >= y0 and b[2] <= x1 and b[3] <= y1]
    if not ink:
        raise SystemExit("no blocks fell inside the region")
    ix0 = min(b[0] for b in ink) - args.pad
    iy0 = min(b[1] for b in ink) - args.pad
    ix1 = max(b[2] for b in ink) + args.pad
    iy1 = max(b[3] for b in ink) + args.pad

    writer = PdfWriter()
    writer.add_page(page)
    new_page = writer.pages[0]
    stream = DecodedStreamObject()
    stream.set_data(b"\n".join(kept))
    new_page[NameObject("/Contents")] = writer._add_object(stream)
    for key in ("/PieceInfo", "/LastModified", "/Thumb", "/Annots"):
        if key in new_page:
            del new_page[key]
    for key in ("/MediaBox", "/CropBox", "/TrimBox", "/BleedBox", "/ArtBox"):
        new_page[NameObject(key)] = RectangleObject([ix0, iy0, ix1, iy1])

    with args.out.open("wb") as handle:
        writer.write(handle)
    print(f"kept {sum(1 for b in ink)} paths -> {args.out}")
    print(f"art box: {ix0:.2f} {iy0:.2f} {ix1:.2f} {iy1:.2f}  "
          f"({ix1 - ix0:.2f} x {iy1 - iy0:.2f} pt)")


if __name__ == "__main__":
    main()
