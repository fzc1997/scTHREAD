#!/usr/bin/env python3
"""Render and compare the bit7 SVG against its source PNG without ImageMagick."""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import cairosvg
import numpy as np
from PIL import Image


def white_rgb(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    rgba = image.convert("RGBA").resize(size, Image.Resampling.LANCZOS)
    white = Image.new("RGBA", size, (255, 255, 255, 255))
    return Image.alpha_composite(white, rgba).convert("RGB")


def compare(reference: Path, candidate: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with Image.open(reference) as image:
        size = image.size
        ref = white_rgb(image, size)
    rendered_bytes = cairosvg.svg2png(
        url=str(candidate),
        output_width=size[0],
        output_height=size[1],
    )
    with Image.open(io.BytesIO(rendered_bytes)) as image:
        render = white_rgb(image, size)

    ref_array = np.asarray(ref, dtype=np.float32)
    render_array = np.asarray(render, dtype=np.float32)
    absolute = np.abs(ref_array - render_array)
    rmse = float(np.sqrt(np.mean(np.square(ref_array - render_array))))
    changed_pixels = int(np.count_nonzero(np.max(absolute, axis=2) > 0.5))
    total_pixels = size[0] * size[1]

    mean_diff = absolute.mean(axis=2)
    red = np.clip(mean_diff * 3.2, 0, 255)
    green = np.clip((mean_diff - 18) * 2.0, 0, 255)
    blue = np.clip((mean_diff - 72) * 1.2, 0, 255)
    heat = np.stack([red, green, blue], axis=2).astype(np.uint8)

    ref.save(out_dir / "ref.png")
    render.save(out_dir / "render.png")
    Image.blend(ref, render, 0.5).save(out_dir / "overlay.png")
    Image.fromarray(heat, mode="RGB").save(out_dir / "diff.png")
    metrics = (
        f"reference={reference}\n"
        f"candidate={candidate}\n"
        f"canvas={size[0]}x{size[1]}\n"
        f"RMSE_8bit={rmse:.6f}\n"
        f"changed_pixels_gt_0.5={changed_pixels}\n"
        f"changed_fraction={changed_pixels / total_pixels:.8f}\n"
    )
    (out_dir / "metrics.txt").write_text(metrics, encoding="utf-8")
    print(metrics, end="")
    for name in ("ref.png", "render.png", "overlay.png", "diff.png", "metrics.txt"):
        path = out_dir / name
        print(f"{path}\t{path.stat().st_size} bytes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    compare(
        args.reference.resolve(),
        args.candidate.resolve(),
        args.out_dir.resolve(),
    )


if __name__ == "__main__":
    main()
