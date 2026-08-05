#!/usr/bin/env python3
"""Compose a pure-vector ChatGPT skeleton with the real-data overlay."""

from __future__ import annotations

import argparse
import copy
import re
from pathlib import Path

import cairosvg
from lxml import etree
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SKELETON = (
    ROOT / "figures/ga_gpt2_vector_v3/chatgpt_skeleton_from_scratch.svg"
)
DEFAULT_OVERLAY = (
    ROOT
    / "figures/ga_gpt2_vector_v3/scTHREAD_ga_gpt2_vector_v3_overlay.svg"
)
DEFAULT_STEM = ROOT / "figures/scTHREAD_graphical_abstract_gpt2_hybrid_v3"
REQUIRED_GROUPS = {
    "source-icons",
    "source-streams",
    "database-cylinder",
    "evidence-cards",
    "query-browser",
    "bottom-access",
    "connectors",
}
FORBIDDEN_TAGS = {
    "image",
    "foreignObject",
    "script",
    "iframe",
    "object",
    "embed",
}
SVG_NS = "http://www.w3.org/2000/svg"


def local_name(value: object) -> str:
    if not isinstance(value, str):
        return ""
    return value.rsplit("}", 1)[-1]


def parse_svg(path: Path) -> etree._Element:
    parser = etree.XMLParser(
        resolve_entities=False,
        no_network=True,
        remove_blank_text=False,
        recover=False,
        huge_tree=True,
    )
    return etree.parse(str(path), parser).getroot()


def validate_svg(root: etree._Element, *, require_groups: bool) -> None:
    if local_name(root.tag) != "svg":
        raise ValueError("Root element is not svg")
    viewbox = [float(value) for value in root.get("viewBox", "").split()]
    if len(viewbox) != 4 or any(
        abs(observed - expected) > 0.01
        for observed, expected in zip(viewbox, [0.0, 0.0, 1916.0, 821.0])
    ):
        raise ValueError(f"Unexpected SVG viewBox: {root.get('viewBox')}")

    observed_ids: set[str] = set()
    for node in root.iter():
        tag = local_name(node.tag)
        if not tag:
            continue
        if tag in FORBIDDEN_TAGS:
            raise ValueError(f"Forbidden SVG element: {tag}")
        node_id = node.get("id")
        if node_id:
            observed_ids.add(node_id)
        for raw_name, value in node.attrib.items():
            name = local_name(raw_name).lower()
            lowered = value.lower()
            if name.startswith("on"):
                raise ValueError(f"Event handler attribute is forbidden: {name}")
            if "javascript:" in lowered or "data:" in lowered:
                raise ValueError(f"Unsafe SVG attribute value: {name}")
            if name == "href" and value and not value.startswith("#"):
                raise ValueError(f"External href is forbidden: {value}")
            if "url(" in lowered:
                for match in re.findall(r"url\(([^)]+)\)", value):
                    target = match.strip().strip("'\"")
                    if not target.startswith("#"):
                        raise ValueError(f"External CSS URL is forbidden: {target}")
        if tag == "style" and node.text:
            css = node.text
            lowered_css = css.lower()
            if "@import" in lowered_css or "javascript:" in lowered_css:
                raise ValueError("Unsafe CSS content is forbidden")
            for match in re.findall(r"url\(([^)]+)\)", css):
                target = match.strip().strip("'\"")
                if not target.startswith("#"):
                    raise ValueError(f"External CSS URL is forbidden: {target}")
    if require_groups:
        missing = REQUIRED_GROUPS - observed_ids
        if missing:
            raise ValueError(f"Required group ids are missing: {sorted(missing)}")


def prefix_overlay_ids(root: etree._Element, prefix: str = "overlay-") -> None:
    id_map: dict[str, str] = {}
    for node in root.iter():
        old_id = node.get("id")
        if old_id:
            new_id = f"{prefix}{old_id}"
            id_map[old_id] = new_id
            node.set("id", new_id)
    for node in root.iter():
        for name, value in list(node.attrib.items()):
            updated = value
            for old_id, new_id in id_map.items():
                updated = updated.replace(f"url(#{old_id})", f"url(#{new_id})")
                if updated == f"#{old_id}":
                    updated = f"#{new_id}"
            if updated != value:
                node.set(name, updated)
        if local_name(node.tag) == "style" and node.text:
            updated_text = node.text
            for old_id, new_id in id_map.items():
                updated_text = updated_text.replace(
                    f"url(#{old_id})", f"url(#{new_id})"
                )
            node.text = updated_text


def compose_svg(skeleton: Path, overlay: Path, output: Path) -> None:
    skeleton_root = parse_svg(skeleton)
    validate_svg(skeleton_root, require_groups=True)
    overlay_root = parse_svg(overlay)
    if local_name(overlay_root.tag) != "svg":
        raise ValueError("Overlay root is not svg")
    prefix_overlay_ids(overlay_root)

    final_root = copy.deepcopy(skeleton_root)
    final_root.set("width", "183mm")
    final_root.set("height", "78mm")
    final_root.set("viewBox", "0 0 1916 821")
    final_root.set("preserveAspectRatio", "none")

    nested = etree.Element(f"{{{SVG_NS}}}svg")
    nested.set("id", "real-data-overlay")
    nested.set("x", "0")
    nested.set("y", "0")
    nested.set("width", "1916")
    nested.set("height", "821")
    nested.set("viewBox", overlay_root.get("viewBox", "0 0 518.740157 221.102362"))
    nested.set("preserveAspectRatio", "none")
    for child in overlay_root:
        nested.append(copy.deepcopy(child))
    final_root.append(nested)
    validate_svg(final_root, require_groups=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    etree.ElementTree(final_root).write(
        str(output),
        encoding="UTF-8",
        xml_declaration=True,
        pretty_print=True,
    )


def export(svg_path: Path, stem: Path, dpi: int) -> list[Path]:
    pdf_path = stem.with_suffix(".pdf")
    png_path = stem.with_suffix(".png")
    cairosvg.svg2pdf(url=str(svg_path), write_to=str(pdf_path))
    width_px = round(183.0 / 25.4 * dpi)
    height_px = round(78.0 / 25.4 * dpi)
    temp_png = stem.with_name(stem.name + ".cairosvg.tmp.png")
    cairosvg.svg2png(
        url=str(svg_path),
        write_to=str(temp_png),
        output_width=width_px,
        output_height=height_px,
    )
    with Image.open(temp_png) as image:
        image.save(png_path, dpi=(dpi, dpi))
    temp_png.unlink()
    return [svg_path, pdf_path, png_path]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skeleton", type=Path, default=DEFAULT_SKELETON)
    parser.add_argument("--overlay", type=Path, default=DEFAULT_OVERLAY)
    parser.add_argument("--stem", type=Path, default=DEFAULT_STEM)
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()

    output_svg = args.stem.resolve().with_suffix(".svg")
    compose_svg(
        args.skeleton.resolve(),
        args.overlay.resolve(),
        output_svg,
    )
    for output in export(output_svg, args.stem.resolve(), args.dpi):
        print(f"{output}\t{output.stat().st_size} bytes")


if __name__ == "__main__":
    main()
