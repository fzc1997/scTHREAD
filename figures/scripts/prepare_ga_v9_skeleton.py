#!/usr/bin/env python3
"""Prepare the metadata-corrected, provenance-ready GA v9 vector skeleton."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

from lxml import etree


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT / "figures/ga_gpt2_vector_v7/chatgpt_skeleton_solid_cards.svg"
)
DEFAULT_OUTPUT = ROOT / "figures/ga_gpt2_vector_v9/chatgpt_skeleton_v9.svg"
INPUT_SHA256 = "a2e896aa33346ee241ba50db9adcd33374e3a08b8662afaef9e30a863d3475a8"


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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def single(root: etree._Element, node_id: str) -> etree._Element:
    nodes = root.xpath(f".//*[@id='{node_id}']")
    require(len(nodes) == 1, f"Expected one #{node_id}, found {len(nodes)}")
    return nodes[0]


def prepare(source: Path, output: Path) -> None:
    require(source.is_file(), f"Missing skeleton: {source}")
    require(
        sha256(source) == INPUT_SHA256,
        f"Frozen v7 skeleton changed: {sha256(source)}",
    )
    parser = etree.XMLParser(
        resolve_entities=False,
        no_network=True,
        remove_blank_text=False,
        recover=False,
        huge_tree=True,
    )
    tree = etree.parse(str(source), parser)
    root = tree.getroot()

    title = single(root, "svgTitle")
    desc = single(root, "svgDesc")
    title.text = "scTHREAD single-cell long-read transcriptome atlas graphical abstract"
    desc.text = (
        "Human and mouse biological coverage feeds a unified cell-resolved "
        "long-read transcriptome database, linked transcriptome evidence, "
        "PTPRC isoform views and open data access."
    )

    hierarchy = single(root, "database-hierarchy")
    hierarchy.getparent().remove(hierarchy)

    # Spread the two lower pictograms across the available biological-coverage
    # width without adding cohort-specific text labels.
    embryo = single(root, "embryo-icon")
    embryo.set("transform", "translate(260 211) scale(0.48)")

    streams = single(root, "source-streams")
    streams.set("opacity", "0.72")

    require(not root.xpath(".//*[local-name()='image']"), "Raster image remains")
    require(not root.xpath(".//*[local-name()='text']"), "Skeleton text remains")
    require(
        not root.xpath(".//*[@id='database-hierarchy']"),
        "Old unlabeled database hierarchy remains",
    )
    require(
        single(root, "embryo-icon").get("transform")
        == "translate(260 211) scale(0.48)",
        "Embryo placement was not updated",
    )
    for required in (
        "source-icons",
        "source-streams",
        "database-cylinder",
        "evidence-cards",
        "query-browser",
        "bottom-access",
        "connectors",
    ):
        single(root, required)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    tree.write(
        str(temporary),
        encoding="UTF-8",
        xml_declaration=True,
        pretty_print=True,
    )
    os.replace(temporary, output)
    print(f"GA V9 SKELETON PASS\t{output}\t{output.stat().st_size}\t{sha256(output)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.input.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
