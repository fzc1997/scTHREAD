#!/usr/bin/env python3
"""Derive Figure 3 cell-type label anchors from the live portal Cell map.

The figure crops the portal's cell-type UMAP canvas verbatim, so a label can
only point at the right cluster if its anchor is expressed in the same canvas
fractions the portal itself draws in.  This script therefore reproduces the
portal's canvas projection exactly (``drawUmapPanel`` in ``web/app.js``) from
two recorded facts: the ``/api/umap/points`` response and the canvas box
measured during the matching screenshot capture.

Anchors are within-cell-type medians, which stay inside their cluster even
when a cluster is crescent-shaped or has outlying cells.

The portal hides points flagged ``iso_outlier`` before drawing (``hideOutliers``
defaults to true in ``web/app.js``), and it derives the canvas extents from what
it draws. Those outliers are by definition the extreme points, so keeping them
would nearly halve ``baseScale`` and move every anchor by up to 0.27 of the
canvas. This script therefore applies the same filter.
"""

from __future__ import annotations

import argparse
import csv
import json
import urllib.request
from pathlib import Path
from statistics import median


PROJECT = Path(__file__).resolve().parents[2]
LABELLED_CELL_TYPES = ("B cell", "Monocyte", "Plasma cell", "Progenitor")
CANVAS_PADDING = 48.0


def fetch_points(base_url: str, gene: str, transcript: str, species: str) -> dict:
    url = (
        f"{base_url.rstrip('/')}/api/umap/points"
        f"?gene={gene}&signal=isoform&transcript={transcript}"
        f"&species={species}&map_scope=sampled"
    )
    with urllib.request.urlopen(url, timeout=300) as response:
        payload = json.load(response)
    if not payload.get("items"):
        raise RuntimeError(f"The portal returned no Cell map points for {gene}")
    return payload


def canvas_box(metadata_path: Path, role: str = "cell") -> dict[str, float]:
    metadata = json.loads(metadata_path.read_text())
    views = [view for view in metadata["views"] if view["view"] == "cellmap"]
    if len(views) != 1:
        raise ValueError(f"Expected one cellmap view in {metadata_path}")
    geometry = views[0]["cellmap_config"]["canvas_geometry"][role]
    return {"width": float(geometry["width"]), "height": float(geometry["height"])}


def canvas_fractions(
    points: list[dict],
    width: float,
    height: float,
    labelled_cell_types: tuple[str, ...] = LABELLED_CELL_TYPES,
    exclude_outliers: bool = True,
) -> dict[str, dict[str, float | int]]:
    """Project points with the portal's own unzoomed, unpanned canvas formula.

    ``exclude_outliers`` reproduces the portal's "Hide outlier points" control.
    That control changes the drawn extent, not just which dots are visible, so a
    capture taken with it enabled must be anchored on the same subset or every
    label lands in the wrong place.
    """
    if exclude_outliers:
        missing = [point for point in points if "iso_outlier" not in point]
        if missing:
            raise ValueError("The points endpoint does not report iso_outlier")
        points = [point for point in points if not point["iso_outlier"]]
        if not points:
            raise ValueError("Every point is flagged as an outlier")
    xs = [float(point["umap1"]) for point in points]
    ys = [float(point["umap2"]) for point in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    base_scale = min(
        (width - CANVAS_PADDING) / max(0.01, max_x - min_x),
        (height - CANVAS_PADDING) / max(0.01, max_y - min_y),
    )
    centre_x, centre_y = (min_x + max_x) / 2, (min_y + max_y) / 2

    grouped: dict[str, list[tuple[float, float]]] = {}
    for point in points:
        grouped.setdefault(str(point["cell_type"]), []).append(
            (float(point["umap1"]), float(point["umap2"]))
        )
    missing = [name for name in labelled_cell_types if name not in grouped]
    if missing:
        raise ValueError(f"The Cell map has no cells for {missing}")

    anchors = {}
    for cell_type in labelled_cell_types:
        cells = grouped[cell_type]
        median_x = median(value for value, _ in cells)
        median_y = median(value for _, value in cells)
        anchors[cell_type] = {
            "n_cells": len(cells),
            "median_umap1": median_x,
            "median_umap2": median_y,
            "canvas_x_fraction": (
                (median_x - centre_x) * base_scale + width / 2
            ) / width,
            "canvas_y_from_top_fraction": (
                (median_y - centre_y) * base_scale + height / 2
            ) / height,
        }
    return anchors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default="http://10.168.3.4:4173")
    parser.add_argument("--gene", default="PTPRC")
    parser.add_argument("--species", default="human")
    parser.add_argument("--transcript", default="ENST00000367364")
    parser.add_argument(
        "--capture-metadata",
        type=Path,
        default=(
            PROJECT
            / "figures"
            / "website_walkthrough"
            / "ptprc_views_v7_scale5_primary"
            / "ptprc_all_views_metadata.json"
        ),
        help="Capture whose canvas box the anchors must match.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--keep-outliers",
        dest="exclude_outliers",
        action="store_false",
        help=(
            "Anchor on every returned point. The portal hides outliers by "
            "default (hideOutliers=true in web/app.js) and derives the canvas "
            "extent from what it draws, so this only matches a capture taken "
            "with that control switched off."
        ),
    )
    parser.set_defaults(exclude_outliers=True)
    parser.add_argument(
        "--cell-types",
        nargs="+",
        default=list(LABELLED_CELL_TYPES),
        help="Cell types to anchor; defaults to the PTPRC figure's four labels.",
    )
    args = parser.parse_args()

    payload = fetch_points(args.base_url, args.gene, args.transcript, args.species)
    box = canvas_box(args.capture_metadata)
    anchors = canvas_fractions(
        payload["items"],
        box["width"],
        box["height"],
        tuple(args.cell_types),
        args.exclude_outliers,
    )
    drawn_cells = (
        sum(1 for point in payload["items"] if not point["iso_outlier"])
        if args.exclude_outliers
        else payload["count"]
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cell_type",
                "n_cells",
                "median_umap1",
                "median_umap2",
                "canvas_x_fraction",
                "canvas_y_from_top_fraction",
            ]
        )
        for cell_type, anchor in anchors.items():
            writer.writerow(
                [
                    cell_type,
                    anchor["n_cells"],
                    anchor["median_umap1"],
                    anchor["median_umap2"],
                    anchor["canvas_x_fraction"],
                    anchor["canvas_y_from_top_fraction"],
                ]
            )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "returned_cells": payload["count"],
                "displayed_cells": drawn_cells,
                "outliers_hidden": args.exclude_outliers,
                "canvas_css_px": box,
                "anchors": anchors,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
