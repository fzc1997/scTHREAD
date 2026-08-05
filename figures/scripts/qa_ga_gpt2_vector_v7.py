#!/usr/bin/env python3
"""Publication, semantic and layout QA for graphical abstract v7."""

from __future__ import annotations

import argparse
from pathlib import Path

import compose_ga_gpt2_vector_v3 as composer
import qa_ga_gpt2_hybrid_v2 as source_qa
import qa_ga_gpt2_vector_v4 as q4
import qa_ga_gpt2_vector_v6 as q6
import prepare_ga_v7_skeleton as skeleton_v7


REQUIRED_TEXT = (
    "scTHREAD: a single-cell long-read transcriptome atlas",
    "Biological coverage",
    "Tissues · disease · development · two species",
    "469 samples",
    "31 studies",
    "845,781 cells",
    "Endocrine",
    "Reproductive",
    "Representative registry datasets · cells",
    "Marrow",
    "Myeloma",
    "Frontal cortex",
    "Retina",
    "Heart",
    "Prostate",
    "Glioma",
    "Ovary",
    "Brain",
    "Fibroblast",
    "Unified database",
    "scTHREAD",
    ">200k isoforms",
    "cell resolved",
    "Four linked RNA evidence layers",
    "Gene expression",
    "across cell types",
    "Isoform usage (DIU)",
    "cell-type transcript shifts",
    "poly(A) site usage",
    "alternative 3′ ends",
    "Allelic expression",
    "ASE + splice junctions",
    "Explore genes & cell types",
    "PTPRC / CD45",
    "PTPRC isoform expression",
    "same 71,913-cell embedding",
    "ENST…67364",
    "ENST…97630",
    "0–5 molecules/cell",
    "0–3 molecules/cell",
    "all positive cells shown",
    "stratified double-zero background",
    "Open data access",
    "Browse · Search · Download",
    "PTPRC / CD45 isoform usage",
    "Genes with cell-type effects",
    "FDR < 0.05 + prespecified effect gate",
    "2,396/10,494 passed",
    "4,168/13,214 passed",
    "0/11,506 passed",
    "673 raw P<0.05",
    "min q=0.639",
)
FORBIDDEN_VISUAL_TERMS = q6.FORBIDDEN_VISUAL_TERMS
RETIRED_PHRASES = (
    *q6.RETIRED_PHRASES,
    "Representative tissues & contexts",
    "cell-resolved transcriptomes",
    "poly(A) site usage (APA)",
    "Allelic expression (ASE)",
)
V6_FROZEN_HASHES = {
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v6.svg"):
        "1b9bd736d4da4f49ca114813290797e44401e1977063f96ec8f2997394443c50",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v6.pdf"):
        "52d4f47b1cdf22aa54ddccced3db6ec01349c96a9535fa7cf47d287c6fa767ab",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v6.png"):
        "48bc2cee24d2fe8e0de0890ab06d39161c7ce93af8753c3c4493337c1ff254ce",
    Path("figures/ga_gpt2_vector_v6/chatgpt_skeleton_compact.svg"):
        "2a550f9cbec1f92409c0eb5d0a13fb297b9f17fc817b8f1b33afcb9826aa710f",
}


def direct_rects(group: object) -> list[object]:
    return [
        child for child in group
        if skeleton_v7.local_name(child.tag) == "rect"
    ]


def check_solid_card_skeleton(path: Path) -> None:
    root = composer.parse_svg(path)
    composer.validate_svg(root, require_groups=True)
    for group_id, (outer_fill, label_fill, outer_stroke) in (
        skeleton_v7.CARD_STYLES.items()
    ):
        group = root.xpath(f".//*[@id='{group_id}']")
        q4.require(len(group) == 1, f"Missing card group: {group_id}")
        rects = direct_rects(group[0])
        q4.require(
            len(rects) == 2,
            f"Expected two direct rectangles in {group_id}, found {len(rects)}",
        )
        outer, label_panel = rects
        q4.require(
            outer.get("fill") == outer_fill
            and outer.get("stroke") == outer_stroke,
            f"Outer solid fill changed in {group_id}",
        )
        q4.require(
            label_panel.get("fill") == label_fill
            and label_panel.get("stroke") == label_fill,
            f"Full-height label fill changed in {group_id}",
        )
        q4.require(
            label_panel.get("fill") != "#ffffff"
            and not label_panel.get("fill", "").startswith("url("),
            f"White/gradient empty label box remains in {group_id}",
        )
    print("SOLID EVIDENCE CARD PASS", skeleton_v7.CARD_STYLES)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v7"),
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
    q6.check_hashes("V6 FREEZE", V6_FROZEN_HASHES)
    q6.check_hashes("V5 FREEZE", q6.V5_FROZEN_HASHES)
    q6.check_hashes("RAW API", q6.RAW_API_HASHES)
    q6.check_hashes("DERIVED SOURCE", q6.DERIVED_HASHES)
    q6.check_compact_skeleton(args.skeleton)
    check_solid_card_skeleton(args.skeleton)
    q6.check_frozen_sources(args.component_dir)
    q4.REQUIRED_TEXT = REQUIRED_TEXT
    q4.FORBIDDEN_VISUAL_TERMS = FORBIDDEN_VISUAL_TERMS
    q4.RETIRED_PHRASES = RETIRED_PHRASES
    exports = q4.check_exports(args.stem, args.skeleton, args.dpi)

    provenance_files = [
        *exports,
        args.skeleton,
        Path("figures/scripts/prepare_ga_v7_skeleton.py"),
        Path("figures/scripts/render_ga_gpt2_vector_v7.py"),
        Path("figures/scripts/compose_ga_gpt2_vector_v3.py"),
        Path("figures/scripts/normalize_ga_vector_v7.py"),
        Path(__file__),
        Path("docs/GRAPHICAL_ABSTRACT_VECTOR_V7_CONTRACT.md"),
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
    print("ALL V7 QA CHECKS PASS")


if __name__ == "__main__":
    main()
