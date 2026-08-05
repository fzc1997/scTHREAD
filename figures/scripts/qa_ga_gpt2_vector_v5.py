#!/usr/bin/env python3
"""Publication, hierarchy and freeze QA for graphical abstract v5."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import compose_ga_gpt2_vector_v3 as composer
import qa_ga_gpt2_hybrid_v2 as source_qa
import qa_ga_gpt2_vector_v4 as q4


REQUIRED_TEXT = (
    "scTHREAD: a single-cell long-read transcriptome atlas",
    "Biological coverage",
    "Human ovary",
    "Mouse gastrula",
    "469 samples",
    "31 studies",
    "845,781 cells",
    "Endocrine",
    "Reproductive",
    "Unified database",
    ">200k isoforms",
    "Four linked RNA evidence layers",
    "Gene expression",
    "Isoform usage (DIU)",
    "poly(A) site usage (APA)",
    "Allelic expression (ASE)",
    "splice junctions",
    "Explore genes & cell types",
    "PTPRC-linked transcriptome evidence",
    "Cell-type landscape",
    "Lymphoid",
    "Myeloid",
    "Erythroid",
    "Progenitor",
    "Neural",
    "Cardiovasc.",
    "Stromal",
    "74,906 portal cells",
    "Open data access",
    "Browse · Search · Download",
    "PTPRC / CD45 isoform usage",
    "Genome-wide RNA events",
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
    "Linked gene evidence",
)
V4_FROZEN_HASHES = {
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v4.svg"):
        "5271d2641cef74fa8a25e00655858466ebbd13068f3532183037efbed2a65c68",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v4.pdf"):
        "933defdecb6b5dae3d60dd316d69eed15d07c67f22559134d796fbf068183a1b",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v4.png"):
        "54c1cf986a3efcb709b71d7bb9201f606d83a7afa8e0157cd156384be4f5e2d2",
}
REMOVED_SKELETON_IDS = {
    "documents-icon",
    "sequencers",
    "browser-tabs",
    "browser-download",
    "bottom-routing",
    "folder-module",
    "search-module",
    "download-module",
    "purple-blank-module",
    "orange-blank-module",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_v4_freeze() -> None:
    for path, expected in V4_FROZEN_HASHES.items():
        q4.require(path.is_file(), f"Frozen v4 file is missing: {path}")
        observed = sha256(path)
        q4.require(
            observed == expected,
            f"Frozen v4 changed: {path} {observed} != {expected}",
        )
    print(
        "V4 FREEZE PASS",
        {str(path): value for path, value in V4_FROZEN_HASHES.items()},
    )


def check_simplified_skeleton(path: Path) -> None:
    root = composer.parse_svg(path)
    composer.validate_svg(root, require_groups=True)
    ids = {node.get("id") for node in root.iter() if node.get("id")}
    present = sorted(REMOVED_SKELETON_IDS & ids)
    q4.require(not present, f"Low-information SVG modules remain: {present}")
    q4.require(
        not any(composer.local_name(node.tag) == "image" for node in root.iter()),
        "Simplified skeleton contains a raster image",
    )
    q4.require(
        not any(composer.local_name(node.tag) == "text" for node in root.iter()),
        "Simplified skeleton is not text-free",
    )
    bottom_access = [
        node for node in root.iter() if node.get("id") == "bottom-access"
    ][0]
    q4.require(
        len(bottom_access) == 0,
        "The old large bottom access modules were not fully removed",
    )
    q4.require(
        all(node.get("stroke-width") != "14" for node in root.iter()),
        "A dominant 14-unit source ribbon remains",
    )
    print(
        "REDUCED-ICON SKELETON PASS",
        {
            "sha256": sha256(path),
            "removed_ids": sorted(REMOVED_SKELETON_IDS),
            "retained_groups": sorted(composer.REQUIRED_GROUPS),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v5"),
    )
    parser.add_argument(
        "--skeleton",
        type=Path,
        default=Path(
            "figures/ga_gpt2_vector_v5/chatgpt_skeleton_simplified.svg"
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
    check_v4_freeze()
    check_simplified_skeleton(args.skeleton)
    q4.REQUIRED_TEXT = REQUIRED_TEXT
    q4.FORBIDDEN_VISUAL_TERMS = FORBIDDEN_VISUAL_TERMS
    q4.RETIRED_PHRASES = RETIRED_PHRASES
    exports = q4.check_exports(args.stem, args.skeleton, args.dpi)

    provenance_files = [
        *exports,
        args.skeleton,
        Path("figures/scripts/prepare_ga_v5_skeleton.py"),
        Path("figures/scripts/render_ga_gpt2_vector_v5.py"),
        Path("figures/scripts/compose_ga_gpt2_vector_v3.py"),
        Path("figures/scripts/normalize_ga_vector_v5.py"),
        Path(__file__),
        Path("docs/GRAPHICAL_ABSTRACT_VECTOR_V5_CONTRACT.md"),
        args.component_dir / "catalog_system_composition_20260726.tsv",
        args.component_dir / "ptprc_two_isoform_switch.tsv",
        args.component_dir / "atlas_umap_stratified_sample.tsv",
        args.component_dir / "three_axis_inventory_source.tsv",
    ]
    print("SHA256")
    for path in provenance_files:
        q4.require(path.is_file(), f"Missing provenance file: {path}")
        print(sha256(path), path)
    print("ALL V5 QA CHECKS PASS")


if __name__ == "__main__":
    main()
