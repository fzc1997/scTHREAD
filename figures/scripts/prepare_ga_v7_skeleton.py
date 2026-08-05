#!/usr/bin/env python3
"""Create the solid-card v7 SVG skeleton from the frozen v6 skeleton."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from lxml import etree

import compose_ga_gpt2_vector_v3 as composer


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "figures/ga_gpt2_vector_v6/chatgpt_skeleton_compact.svg"
DEFAULT_OUTPUT = (
    ROOT / "figures/ga_gpt2_vector_v7/chatgpt_skeleton_solid_cards.svg"
)
FROZEN_INPUT_SHA256 = (
    "2a550f9cbec1f92409c0eb5d0a13fb297b9f17fc817b8f1b33afcb9826aa710f"
)

CARD_STYLES = {
    "card-cells-bars": ("#E6F4F3", "#197A83", "#73B8BA"),
    "card-network": ("#EEE9F7", "#5B3F99", "#9F8AC8"),
    "card-tracks": ("#EAF2FA", "#326FA5", "#82A9D2"),
    "card-arcs": ("#FFF0E7", "#C8574E", "#E7A16E"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(tag: object) -> str:
    if not isinstance(tag, str):
        return ""
    return tag.rsplit("}", 1)[-1]


def by_id(root: etree._Element, node_id: str) -> etree._Element:
    nodes = root.xpath(f".//*[@id='{node_id}']")
    if len(nodes) != 1:
        raise ValueError(
            f"Expected one SVG node with id={node_id!r}, found {len(nodes)}"
        )
    return nodes[0]


def solidify_cards(input_svg: Path, output_svg: Path) -> dict[str, object]:
    observed = sha256(input_svg)
    if observed != FROZEN_INPUT_SHA256:
        raise ValueError(
            f"Frozen v6 skeleton changed: {observed} != {FROZEN_INPUT_SHA256}"
        )

    root = composer.parse_svg(input_svg)
    composer.validate_svg(root, require_groups=True)
    changed: dict[str, dict[str, str]] = {}
    for group_id, (outer_fill, label_fill, outer_stroke) in CARD_STYLES.items():
        group = by_id(root, group_id)
        direct_rects = [
            child for child in group
            if local_name(child.tag) == "rect"
        ]
        if len(direct_rects) != 2:
            raise ValueError(
                f"Expected two direct card rectangles in {group_id}, "
                f"found {len(direct_rects)}"
            )
        outer, label_panel = direct_rects
        outer.set("fill", outer_fill)
        outer.set("stroke", outer_stroke)
        outer.set("stroke-width", "1.6")
        label_panel.set("fill", label_fill)
        label_panel.set("stroke", label_fill)
        label_panel.set("stroke-width", "0")
        changed[group_id] = {
            "outer_fill": outer_fill,
            "label_fill": label_fill,
            "outer_stroke": outer_stroke,
        }

    composer.validate_svg(root, require_groups=True)
    if any(composer.local_name(node.tag) == "image" for node in root.iter()):
        raise ValueError("V7 skeleton unexpectedly contains a raster image")
    if any(composer.local_name(node.tag) == "text" for node in root.iter()):
        raise ValueError("V7 skeleton unexpectedly contains text")

    output_svg.parent.mkdir(parents=True, exist_ok=True)
    etree.ElementTree(root).write(
        str(output_svg),
        encoding="UTF-8",
        xml_declaration=True,
        pretty_print=True,
    )
    return {
        "input": str(input_svg),
        "output": str(output_svg),
        "input_sha256": observed,
        "output_sha256": sha256(output_svg),
        "solid_cards": changed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(solidify_cards(args.input.resolve(), args.output.resolve()))


if __name__ == "__main__":
    main()
