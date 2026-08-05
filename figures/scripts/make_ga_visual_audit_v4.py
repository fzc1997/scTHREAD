#!/usr/bin/env python3
"""Create Python-only thumbnail and grayscale review assets for GA v4."""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image


def save_review_assets(source: Path, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    with Image.open(source) as image:
        rgba = image.convert("RGBA")
        for label, scale in (("thumbnail_50", 0.50), ("thumbnail_25", 0.25)):
            size = (
                max(1, round(rgba.width * scale)),
                max(1, round(rgba.height * scale)),
            )
            resized = rgba.resize(size, Image.Resampling.LANCZOS)
            path = output_dir / f"{label}.png"
            resized.save(path, dpi=(450, 450))
            outputs.append(path)

        grey = rgba.convert("L").convert("RGBA")
        size = (round(rgba.width * 0.50), round(rgba.height * 0.50))
        grey = grey.resize(size, Image.Resampling.LANCZOS)
        path = output_dir / "grayscale_50.png"
        grey.save(path, dpi=(450, 450))
        outputs.append(path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v4.png"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not args.source.is_file():
        raise FileNotFoundError(args.source)
    for output in save_review_assets(args.source, args.output_dir):
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
