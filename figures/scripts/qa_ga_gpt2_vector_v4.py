#!/usr/bin/env python3
"""Publication, semantic and freeze QA for the pure-vector graphical abstract v4."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
from pathlib import Path

from PIL import Image

import compose_ga_gpt2_vector_v3 as composer
import qa_ga_gpt2_hybrid_v2 as source_qa


REQUIRED_TEXT = (
    "scTHREAD: a single-cell long-read transcriptome atlas",
    "Human ovary",
    "Mouse gastrula",
    "Biological breadth",
    "Endocrine",
    "Reproductive",
    "Unified database",
    "469 samples",
    "31 studies",
    "845,781",
    ">200k isoforms",
    "Four linked RNA evidence layers",
    "Gene expression",
    "Isoform usage (DIU)",
    "poly(A) site usage (APA)",
    "Allelic expression (ASE)",
    "splice junctions",
    "Explore genes & cell types",
    "Cell-type landscape",
    "Lymphoid",
    "Myeloid",
    "Erythroid",
    "Progenitor",
    "Neural",
    "Cardiovasc.",
    "Stromal",
    "74,906 portal cells",
    "PTPRC isoform usage",
    "Genome-wide events",
)
FORBIDDEN_VISUAL_TERMS = (
    "benchmark",
    "method",
    "Smart-seq2",
    "differentiation",
)
RETIRED_PHRASES = (
    "uniform IsoQuant reprocessing",
    "Uniform reprocessing",
    "Precomputed maps",
    "Query & online analysis",
    "4 aligned RNA evidence layers",
    "Multi-layer gene card",
)
V3_FROZEN_HASHES = {
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v3.svg"):
        "b8ae12b7261b0ee89dab8f2e0373c0a46eb2ce0bf8691132addfb05686e2d77e",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v3.pdf"):
        "bab9375e5817472ed126d0982d88d38ba448f5e0f3ddf604288687c775e0b99d",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v3.png"):
        "75fe30b76e47f8980e005f213e10906139d4c77ddf0c4860c184d9d1bae83655",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def command_output(*args: str) -> str:
    return subprocess.run(
        args,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_v3_freeze() -> None:
    for path, expected in V3_FROZEN_HASHES.items():
        require(path.is_file(), f"Frozen v3 file is missing: {path}")
        observed = sha256(path)
        require(
            observed == expected,
            f"Frozen v3 changed: {path} {observed} != {expected}",
        )
    print("V3 FREEZE PASS", {str(path): value for path, value in V3_FROZEN_HASHES.items()})


def check_svg(svg_path: Path, skeleton_path: Path) -> str:
    skeleton_root = composer.parse_svg(skeleton_path)
    composer.validate_svg(skeleton_root, require_groups=True)
    skeleton_text = [
        node for node in skeleton_root.iter()
        if composer.local_name(node.tag) == "text"
    ]
    require(not skeleton_text, "Frozen bit7 skeleton is no longer text-free")

    root = composer.parse_svg(svg_path)
    composer.validate_svg(root, require_groups=True)
    ids = {node.get("id") for node in root.iter() if node.get("id")}
    require("real-data-overlay" in ids, "Real-data overlay group is missing")
    image_count = sum(
        composer.local_name(node.tag) == "image" for node in root.iter()
    )
    text_nodes = [
        node for node in root.iter()
        if composer.local_name(node.tag) == "text"
    ]
    require(image_count == 0, f"Final SVG contains {image_count} raster image nodes")
    require(len(text_nodes) >= 80, f"Too few editable text nodes: {len(text_nodes)}")
    svg_text = " ".join("".join(node.itertext()) for node in text_nodes)
    print(
        "SVG EDITABILITY PASS",
        {
            "image_elements": image_count,
            "text_elements": len(text_nodes),
            "required_groups": sorted(composer.REQUIRED_GROUPS),
        },
    )
    return svg_text


def check_exports(stem: Path, skeleton_path: Path, dpi: int) -> list[Path]:
    svg_path = stem.with_suffix(".svg")
    pdf_path = stem.with_suffix(".pdf")
    png_path = stem.with_suffix(".png")
    paths = [svg_path, pdf_path, png_path]
    for path in paths:
        require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    svg_text = check_svg(svg_path, skeleton_path)
    pdf_text = command_output("pdftotext", str(pdf_path), "-")
    combined_text = f"{svg_text}\n{pdf_text}"
    for label in REQUIRED_TEXT:
        require(label in combined_text, f"Required label is missing: {label}")
    for term in FORBIDDEN_VISUAL_TERMS:
        require(
            re.search(rf"\b{re.escape(term)}\b", combined_text, re.I) is None,
            f"Forbidden visual term remains: {term}",
        )
    for phrase in RETIRED_PHRASES:
        require(phrase not in combined_text, f"Retired phrase remains: {phrase}")

    pdf_info = command_output("pdfinfo", str(pdf_path))
    require("Pages:           1" in pdf_info, "PDF must contain one page")
    page_match = re.search(
        r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info
    )
    require(page_match is not None, "Could not read PDF page size")
    width_pt, height_pt = map(float, page_match.groups())
    require(abs(width_pt - 518.74) < 0.2, f"Unexpected PDF width: {width_pt}")
    require(abs(height_pt - 221.102) < 0.2, f"Unexpected PDF height: {height_pt}")

    font_info = command_output("pdffonts", str(pdf_path))
    require("Type 3" not in font_info, "Type 3 fonts are not allowed")
    require(
        re.search(r"Arial|LiberationSans|Helvetica", font_info, re.I) is not None,
        "A publication-safe sans-serif font was not embedded",
    )
    pdf_images = command_output("pdfimages", "-list", str(pdf_path))
    require(
        re.search(r"(?m)^\s*\d+\s+\d+\s+", pdf_images) is None,
        "PDF contains an embedded raster image",
    )

    with Image.open(png_path) as png:
        expected = (
            round(183.0 / 25.4 * dpi),
            round(78.0 / 25.4 * dpi),
        )
        require(png.size == expected, f"Unexpected PNG size: {png.size}")
        metadata_dpi = png.info.get("dpi")
        require(metadata_dpi is not None, "PNG DPI metadata is missing")
        require(
            all(abs(value - dpi) < 0.2 for value in metadata_dpi),
            f"PNG DPI is not {dpi}: {metadata_dpi}",
        )
        print("PNG SIZE/DPI PASS", png.size, metadata_dpi)

    print("PDF INFO")
    print(pdf_info.strip())
    print("PDF FONTS")
    print(font_info.strip())
    print("PDF RASTER CHECK")
    print(pdf_images.strip())
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v4"),
    )
    parser.add_argument(
        "--skeleton",
        type=Path,
        default=Path(
            "figures/ga_gpt2_vector_v3/chatgpt_skeleton_from_scratch.svg"
        ),
    )
    parser.add_argument(
        "--component-dir",
        type=Path,
        default=Path("figures/ga_gpt2_components_v2"),
    )
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()

    source_qa.check_catalog(args.component_dir)
    check_v3_freeze()
    exports = check_exports(args.stem, args.skeleton, args.dpi)
    provenance_files = [
        *exports,
        args.skeleton,
        Path("figures/scripts/render_ga_gpt2_vector_v4.py"),
        Path("figures/scripts/compose_ga_gpt2_vector_v3.py"),
        Path(__file__),
        Path("docs/GRAPHICAL_ABSTRACT_VECTOR_V4_CONTRACT.md"),
        args.component_dir / "catalog_system_composition_20260726.tsv",
        args.component_dir / "ptprc_two_isoform_switch.tsv",
        args.component_dir / "atlas_umap_stratified_sample.tsv",
        args.component_dir / "three_axis_inventory_source.tsv",
    ]
    print("SHA256")
    for path in provenance_files:
        require(path.is_file(), f"Missing provenance file: {path}")
        print(sha256(path), path)
    print("ALL V4 QA CHECKS PASS")


if __name__ == "__main__":
    main()
