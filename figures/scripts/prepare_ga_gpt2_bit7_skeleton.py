#!/usr/bin/env python3
"""Normalize the raw bit7 SVG while preserving its visible vector artwork."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from lxml import etree

import compose_ga_gpt2_vector_v3 as composer


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "figures/ga_gpt2_vector_v3/chatgpt_skeleton_bit7_raw.svg"
)
DEFAULT_OUTPUT = (
    ROOT / "figures/ga_gpt2_vector_v3/chatgpt_skeleton_from_scratch.svg"
)
EXPECTED_RAW_SHA256 = (
    "94f4ae81133c7717d91c4b550829338c2fe1be7ab1f31252ccd2e9c6f306c470"
)
ID_MAP = {
    "biological-sources": "source-icons",
    "read-ribbons": "source-streams",
    "database": "database-cylinder",
    "analytic-cards": "evidence-cards",
    "browser-window": "query-browser",
    "bottom-workflow": "bottom-access",
    "right-routing": "connectors",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(input_path: Path, output_path: Path) -> None:
    observed_sha = sha256(input_path)
    if observed_sha != EXPECTED_RAW_SHA256:
        raise RuntimeError(
            "Raw bit7 SVG checksum changed: "
            f"{observed_sha} != {EXPECTED_RAW_SHA256}"
        )
    root = composer.parse_svg(input_path)
    composer.validate_svg(root, require_groups=False)

    ids = {node.get("id"): node for node in root.iter() if node.get("id")}
    placeholder = ids.get("editable-text-placeholders")
    if placeholder is None:
        raise RuntimeError("Expected invisible text placeholder group is missing")
    placeholder.getparent().remove(placeholder)

    for old_id, new_id in ID_MAP.items():
        node = ids.get(old_id)
        if node is None:
            raise RuntimeError(f"Expected bit7 group is missing: {old_id}")
        node.set("id", new_id)

    visible_text = [
        node for node in root.iter()
        if composer.local_name(node.tag) == "text"
    ]
    if visible_text:
        raise RuntimeError(
            f"Unexpected text nodes remain in the skeleton: {len(visible_text)}"
        )
    composer.validate_svg(root, require_groups=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    etree.ElementTree(root).write(
        str(output_path),
        encoding="UTF-8",
        xml_declaration=True,
        pretty_print=True,
    )
    print(
        {
            "input": str(input_path),
            "output": str(output_path),
            "raw_sha256": observed_sha,
            "output_sha256": sha256(output_path),
            "text_nodes": len(visible_text),
            "required_groups": sorted(composer.REQUIRED_GROUPS),
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.input.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
