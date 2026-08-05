#!/usr/bin/env python3
"""Publication, semantic and reading-order QA for graphical abstract v8."""

from __future__ import annotations

import argparse
from pathlib import Path

import compose_ga_gpt2_vector_v3 as composer
import qa_ga_gpt2_hybrid_v2 as source_qa
import qa_ga_gpt2_vector_v4 as q4
import qa_ga_gpt2_vector_v6 as q6
import qa_ga_gpt2_vector_v7 as q7


REQUIRED_TEXT = (
    *q7.REQUIRED_TEXT,
    "HUMAN TISSUES",
    "MOUSE / DEVELOPMENT",
    "ASE: 673 raw P<0.05 · min q=0.639",
)
FORBIDDEN_VISUAL_TERMS = q7.FORBIDDEN_VISUAL_TERMS
RETIRED_PHRASES = q7.RETIRED_PHRASES
V7_FROZEN_HASHES = {
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v7.svg"):
        "964186243baab899a291fb06559dbbe9b8c8721e2b7d313e61054e7ee8a9ef1a",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v7.pdf"):
        "ca75387278b5dd9d07b3af2423d5e9b392f767e223512727fb547ed23bdc8b6f",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v7.png"):
        "a34662473d7efda2ad3964a47012f357a72726982bdf647ab6c6f2f8939bd607",
    Path("figures/ga_gpt2_vector_v7/chatgpt_skeleton_solid_cards.svg"):
        "a2e896aa33346ee241ba50db9adcd33374e3a08b8662afaef9e30a863d3475a8",
}

TISSUE_ROWS = (
    ("Frontal cortex", "Marrow"),
    ("Heart", "Myeloma"),
    ("Brain", "Fibroblast→HSC"),
    ("Retina", "Prostate"),
    ("Glioma", "Ovary"),
)


def text_coordinates(path: Path, labels: set[str]) -> dict[str, tuple[float, float]]:
    root = composer.parse_svg(path)
    found: dict[str, list[tuple[float, float]]] = {label: [] for label in labels}
    for node in root.iter():
        if composer.local_name(node.tag) != "text":
            continue
        value = "".join(node.itertext())
        if value not in labels:
            continue
        x = float(node.get("x", "nan"))
        y = float(node.get("y", "nan"))
        found[value].append((x, y))
    result: dict[str, tuple[float, float]] = {}
    tissue_labels = {label for pair in TISSUE_ROWS for label in pair}
    for label, positions in found.items():
        if label in tissue_labels:
            positions = [
                position for position in positions
                if position[0] < 130 and position[1] > 105
            ]
        q4.require(
            len(positions) == 1,
            f"Expected one positioned text node for {label!r}, found {positions}",
        )
        result[label] = positions[0]
    return result


def check_reading_order(svg_path: Path) -> None:
    labels = {
        *(label for pair in TISSUE_ROWS for label in pair),
        "0/11,506 passed",
        "ASE: 673 raw P<0.05 · min q=0.639",
    }
    coordinates = text_coordinates(svg_path, labels)
    row_y: list[float] = []
    for left, right in TISSUE_ROWS:
        lx, ly = coordinates[left]
        rx, ry = coordinates[right]
        q4.require(lx < rx, f"Tissue pair is not left-to-right: {left}, {right}")
        q4.require(
            abs(ly - ry) < 0.05,
            f"Tissue pair is not aligned: {left} y={ly}, {right} y={ry}",
        )
        row_y.append((ly + ry) / 2)
    q4.require(
        all(first < second for first, second in zip(row_y, row_y[1:])),
        f"Tissue rows are not top-to-bottom: {row_y}",
    )

    _, ase_y = coordinates["0/11,506 passed"]
    _, note_y = coordinates["ASE: 673 raw P<0.05 · min q=0.639"]
    q4.require(
        note_y - ase_y > 2.5,
        f"ASE note remains too close to adjusted result: {ase_y}, {note_y}",
    )
    print(
        "READING ORDER PASS",
        {
            "tissue_rows": TISSUE_ROWS,
            "row_y": row_y,
            "ASE_result_y": ase_y,
            "ASE_note_y": note_y,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v8"),
    )
    parser.add_argument(
        "--skeleton",
        type=Path,
        default=Path(
            "figures/ga_gpt2_vector_v7/chatgpt_skeleton_solid_cards.svg"
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
    q6.check_hashes("V7 FREEZE", V7_FROZEN_HASHES)
    q6.check_hashes("V6 FREEZE", q7.V6_FROZEN_HASHES)
    q6.check_hashes("RAW API", q6.RAW_API_HASHES)
    q6.check_hashes("DERIVED SOURCE", q6.DERIVED_HASHES)
    q6.check_compact_skeleton(args.skeleton)
    q7.check_solid_card_skeleton(args.skeleton)
    q6.check_frozen_sources(args.component_dir)
    q4.REQUIRED_TEXT = REQUIRED_TEXT
    q4.FORBIDDEN_VISUAL_TERMS = FORBIDDEN_VISUAL_TERMS
    q4.RETIRED_PHRASES = RETIRED_PHRASES
    exports = q4.check_exports(args.stem, args.skeleton, args.dpi)
    check_reading_order(args.stem.with_suffix(".svg"))

    provenance_files = [
        *exports,
        args.skeleton,
        Path("figures/scripts/render_ga_gpt2_vector_v8.py"),
        Path("figures/scripts/compose_ga_gpt2_vector_v3.py"),
        Path("figures/scripts/normalize_ga_vector_v8.py"),
        Path(__file__),
        Path("docs/GRAPHICAL_ABSTRACT_VECTOR_V8_CONTRACT.md"),
        args.component_dir / "catalog_system_composition_20260726.tsv",
        args.component_dir / "ptprc_two_isoform_switch.tsv",
        args.component_dir / "three_axis_inventory_source.tsv",
        *q6.RAW_API_HASHES,
        *q6.DERIVED_HASHES,
    ]
    print("SHA256")
    for path in provenance_files:
        q4.require(path.is_file(), f"Missing provenance file: {path}")
        print(q6.sha256(path), path)
    print("ALL V8 QA CHECKS PASS")


if __name__ == "__main__":
    main()
