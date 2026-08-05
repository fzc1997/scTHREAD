#!/usr/bin/env python3
"""Normalise vector art extracted from an Illustrator file into a tintable icon.

``pdftocairo -svg`` faithfully reproduces the source artwork, including its
literal colours, clip paths, masks and generated element ids.  To use that output
alongside the healthicons/Tabler icons in the graphical-abstract strip it has to
be:

1. recoloured to ``currentColor``, so one group ``color`` tints the whole icon;
2. id-namespaced, so its ``clipPath``/``mask``/``filter`` ids cannot collide with
   anything else once inlined into the composite overlay.

Geometry is untouched -- the path data stays byte-identical to the extraction.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def normalize(svg: str, prefix: str) -> str:
    # 1. every literal colour becomes currentColor; fill:none / stroke:none stay.
    svg = re.sub(r"rgb\(\s*[\d.]+%\s*,\s*[\d.]+%\s*,\s*[\d.]+%\s*\)", "currentColor", svg)

    # 2. namespace generated ids and every reference to them.
    ids = sorted(set(re.findall(r'\bid="([^"]+)"', svg)), key=len, reverse=True)
    for old in ids:
        new = f"{prefix}-{old}"
        svg = svg.replace(f'id="{old}"', f'id="{new}"')
        svg = svg.replace(f"url(#{old})", f"url(#{new})")
        svg = svg.replace(f'xlink:href="#{old}"', f'xlink:href="#{new}"')
        svg = svg.replace(f'href="#{old}"', f'href="#{new}"')

    # 3. drop the XML prolog so the fragment can be inlined.
    svg = re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", svg)
    return svg


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    svg = normalize(args.src.read_text(encoding="utf-8"), args.prefix)
    if "viewBox" not in svg:
        raise SystemExit("normalised SVG has no viewBox; the loader needs one")
    if re.search(r"rgb\(", svg):
        raise SystemExit("a literal colour survived normalisation")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(svg, encoding="utf-8")
    print(f"wrote {args.out} ({len(svg)} bytes, ids prefixed '{args.prefix}-')")


if __name__ == "__main__":
    main()
