#!/usr/bin/env python3
"""Create the compact, tissue-aware v6 SVG skeleton from frozen v5 artwork."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from lxml import etree

import compose_ga_gpt2_vector_v3 as composer


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "figures/ga_gpt2_vector_v5/chatgpt_skeleton_simplified.svg"
)
DEFAULT_OUTPUT = (
    ROOT / "figures/ga_gpt2_vector_v6/chatgpt_skeleton_compact.svg"
)
FROZEN_INPUT_SHA256 = (
    "9c64ce7223ff3c577c850ab722e9ce0909bf310bc2ead63470a6490a8d1fda98"
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
        raise ValueError(
            f"Expected one SVG node with id={node_id!r}, found {len(nodes)}"
        )
    return nodes[0]


def remove_by_id(root: etree._Element, node_id: str) -> None:
    node = by_id(root, node_id)
    parent = node.getparent()
    if parent is None:
        raise ValueError(f"Cannot remove root node: {node_id}")
    parent.remove(node)


def compact(input_svg: Path, output_svg: Path) -> dict[str, object]:
    observed = sha256(input_svg)
    if observed != FROZEN_INPUT_SHA256:
        raise ValueError(
            f"Frozen v5 skeleton changed: {observed} != {FROZEN_INPUT_SHA256}"
        )

    root = composer.parse_svg(input_svg)
    composer.validate_svg(root, require_groups=True)

    # Keep ovary and embryo biology visible, but subordinate every pictogram to
    # the evidence panels. Lungs are removed because the frozen registry text
    # used for v6 does not directly support a lung label.
    transforms = {
        "human-icon": "translate(42 118) scale(0.44)",
        "brain-icon": "translate(108 118) scale(0.52)",
        "blood-drop": "translate(185 122) scale(0.52)",
        "heart-icon": "translate(244 116) scale(0.44)",
        "reproductive-icon": "translate(313 116) scale(0.42)",
        "mouse-icon": "translate(47 220) scale(0.48)",
        "embryo-icon": "translate(159 211) scale(0.52)",
    }
    for node_id, transform in transforms.items():
        by_id(root, node_id).set("transform", transform)
    remove_by_id(root, "lungs-icon")

    source_icons = by_id(root, "source-icons")
    source_icons.set("opacity", "0.92")

    composer.validate_svg(root, require_groups=True)
    if any(composer.local_name(node.tag) == "image" for node in root.iter()):
        raise ValueError("Compact skeleton unexpectedly contains a raster image")
    if any(composer.local_name(node.tag) == "text" for node in root.iter()):
        raise ValueError("Compact skeleton unexpectedly contains text")

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
        "removed": ["lungs-icon"],
        "rescaled": sorted(transforms),
        "retained_groups": sorted(composer.REQUIRED_GROUPS),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(compact(args.input.resolve(), args.output.resolve()))


if __name__ == "__main__":
    main()
