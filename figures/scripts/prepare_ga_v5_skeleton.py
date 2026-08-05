#!/usr/bin/env python3
"""Create a quieter v5 SVG skeleton from the frozen bit7 vector artwork."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from lxml import etree

import compose_ga_gpt2_vector_v3 as composer


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "figures/ga_gpt2_vector_v3/chatgpt_skeleton_from_scratch.svg"
)
DEFAULT_OUTPUT = (
    ROOT / "figures/ga_gpt2_vector_v5/chatgpt_skeleton_simplified.svg"
)
FROZEN_INPUT_SHA256 = (
    "e30097c8d5ffb9f4c1ddc086f75e164e8ba2dd862f6fc846b39344d5dce6f83a"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def by_id(root: etree._Element, node_id: str) -> etree._Element:
    nodes = root.xpath(f".//*[@id='{node_id}']")
    if len(nodes) != 1:
        raise ValueError(f"Expected one SVG node with id={node_id!r}, found {len(nodes)}")
    return nodes[0]


def remove_by_id(root: etree._Element, node_id: str) -> None:
    node = by_id(root, node_id)
    parent = node.getparent()
    if parent is None:
        raise ValueError(f"Cannot remove root node: {node_id}")
    parent.remove(node)


def simplify(input_svg: Path, output_svg: Path) -> dict[str, object]:
    observed = sha256(input_svg)
    if observed != FROZEN_INPUT_SHA256:
        raise ValueError(
            f"Frozen bit7 skeleton changed: {observed} != {FROZEN_INPUT_SHA256}"
        )

    root = composer.parse_svg(input_svg)
    composer.validate_svg(root, require_groups=True)

    transforms = {
        "human-icon": "translate(45 112) scale(0.62)",
        "brain-icon": "translate(125 112) scale(0.72)",
        "blood-drop": "translate(225 115) scale(0.72)",
        "heart-icon": "translate(125 178) scale(0.62)",
        "lungs-icon": "translate(200 178) scale(0.62)",
        "reproductive-icon": "translate(275 175) scale(0.62)",
        "mouse-icon": "translate(55 258) scale(0.66)",
        "embryo-icon": "translate(165 246) scale(0.72)",
    }
    for node_id, transform in transforms.items():
        by_id(root, node_id).set("transform", transform)

    for node_id in (
        "documents-icon",
        "sequencers",
        "browser-tabs",
        "browser-download",
        "bottom-routing",
    ):
        remove_by_id(root, node_id)

    source_icons = by_id(root, "source-icons")
    for child in list(source_icons):
        if composer.local_name(child.tag) != "rect":
            continue
        if child.get("x") == "16" and child.get("y") == "469":
            source_icons.remove(child)

    source_streams = by_id(root, "source-streams")
    source_streams.set("opacity", "0.82")
    for node in source_streams.iter():
        width = node.get("stroke-width")
        if width == "14":
            node.set("stroke-width", "8")
            node.set("opacity", "0.58")
        elif width == "2":
            node.set("stroke-width", "1.6")
        elif width == "1":
            node.set("stroke-width", "0.8")

    bottom_access = by_id(root, "bottom-access")
    for child in list(bottom_access):
        bottom_access.remove(child)

    query_browser = by_id(root, "query-browser")
    for child in list(query_browser):
        if composer.local_name(child.tag) != "circle":
            continue
        if child.get("cx") in {"1496", "1519", "1542"}:
            query_browser.remove(child)

    composer.validate_svg(root, require_groups=True)
    if any(composer.local_name(node.tag) == "image" for node in root.iter()):
        raise ValueError("Simplified skeleton unexpectedly contains a raster image")
    if any(composer.local_name(node.tag) == "text" for node in root.iter()):
        raise ValueError("Simplified skeleton unexpectedly contains text")

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
        "removed": [
            "documents-icon",
            "sequencers",
            "browser-tabs",
            "browser-download",
            "bottom-routing",
            "large bottom access modules",
        ],
        "retained_groups": sorted(composer.REQUIRED_GROUPS),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(simplify(args.input.resolve(), args.output.resolve()))


if __name__ == "__main__":
    main()
