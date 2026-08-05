#!/usr/bin/env python3
"""Publication, source-integrity and freeze QA for graphical abstract v6."""

from __future__ import annotations

import os

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

import compose_ga_gpt2_vector_v3 as composer
import qa_ga_gpt2_hybrid_v2 as source_qa
import qa_ga_gpt2_vector_v4 as q4


REQUIRED_TEXT = (
    "scTHREAD: a single-cell long-read transcriptome atlas",
    "Biological coverage",
    "469 samples",
    "31 studies",
    "845,781 cells",
    "Endocrine",
    "Reproductive",
    "Representative tissues & contexts",
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
    ">200k isoforms",
    "Four linked RNA evidence layers",
    "Gene expression",
    "Isoform usage (DIU)",
    "poly(A) site usage (APA)",
    "Allelic expression (ASE)",
    "splice junctions",
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
FORBIDDEN_VISUAL_TERMS = (
    "benchmark",
    "method",
    "Smart-seq2",
    "differentiation",
)
RETIRED_PHRASES = (
    "Human ovary (scONT)",
    "Mouse gastrula (scONT)",
    "Genome-wide RNA events",
    "Genome-wide events",
    "Cell-type landscape",
    "PTPRC-linked transcriptome evidence",
    "74,906 portal cells",
    "uniform IsoQuant reprocessing",
    "Uniform reprocessing",
    "Precomputed maps",
    "Query & online analysis",
    "4 aligned RNA evidence layers",
    "Multi-layer gene card",
    "Linked gene evidence",
)
V5_FROZEN_HASHES = {
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v5.svg"):
        "3d33b881dd54c4dd10614529e79dc407a72567a8e14f8a5d4bb7b1ba667b01f4",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v5.pdf"):
        "77c027b1f849805aa25efecbc71de29bdabe3f74929b494fdebd7c90b23a604f",
    Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v5.png"):
        "143ac3103dea74d188af5ce823bf761a1ef10414c7c178c7aca9be3f181d2f29",
    Path("figures/ga_gpt2_vector_v5/chatgpt_skeleton_simplified.svg"):
        "9c64ce7223ff3c577c850ab722e9ce0909bf310bc2ead63470a6490a8d1fda98",
}
RAW_API_HASHES = {
    Path("figures/ga_gpt2_vector_v6/raw_api/ENST00000367364.json"):
        "fbf4086d97602b183b96119e6657598f2391a98f26ae96501ac3b8160a08a086",
    Path("figures/ga_gpt2_vector_v6/raw_api/ENST00000697630.json"):
        "e6c60b580817af64057f8b9db811b017c4c4e4ba6f5526c54d4356a4c3c8bdfb",
}
DERIVED_HASHES = {
    Path("figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_points.tsv"):
        "686360acf4550507fc8bc5bdc4c57568ef4b27d5697680d6572d9496b29b6364",
    Path("figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_stats.json"):
        "547f7d5f6d2fa5fb95312e4066bb6d15f1ab2ce67b91b39d892484c5fb48730c",
    Path("figures/ga_gpt2_vector_v6/registry_tissue_examples.tsv"):
        "e3a8cfdb8945f1e65b2801340742c6d18f0508cdd51dd1fb1d2184060aaa4ae1",
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
    "lungs-icon",
}
REGISTRY_FROZEN = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
ASE_FULL = Path(
    os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD") + "/results/paper1/f2_grammar/"
    "figdata/ase_interaction.tsv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_hashes(label: str, expected: dict[Path, str]) -> None:
    observed: dict[str, str] = {}
    for path, target in expected.items():
        q4.require(path.is_file(), f"{label} file is missing: {path}")
        value = sha256(path)
        q4.require(value == target, f"{label} changed: {path} {value} != {target}")
        observed[str(path)] = value
    print(f"{label} HASH PASS", observed)


def check_compact_skeleton(path: Path) -> None:
    root = composer.parse_svg(path)
    composer.validate_svg(root, require_groups=True)
    ids = {node.get("id") for node in root.iter() if node.get("id")}
    present = sorted(REMOVED_SKELETON_IDS & ids)
    q4.require(not present, f"Removed SVG modules remain: {present}")
    q4.require(
        not any(composer.local_name(node.tag) == "image" for node in root.iter()),
        "Compact skeleton contains a raster image",
    )
    q4.require(
        not any(composer.local_name(node.tag) == "text" for node in root.iter()),
        "Compact skeleton is not text-free",
    )
    source_icons = [node for node in root.iter() if node.get("id") == "source-icons"][0]
    scales: dict[str, float] = {}
    for node in source_icons.iter():
        node_id = node.get("id")
        transform = node.get("transform", "")
        match = (
            re.search(r"scale\(([0-9.]+)\)", transform)
            if transform
            else None
        )
        if node_id and match:
            scales[node_id] = float(match.group(1))
    q4.require(scales, "No scaled source pictograms were found")
    q4.require(
        max(scales.values()) <= 0.52,
        f"A source pictogram remains too large: {scales}",
    )
    print(
        "COMPACT SKELETON PASS",
        {
            "sha256": sha256(path),
            "pictogram_scales": scales,
            "retained_groups": sorted(composer.REQUIRED_GROUPS),
        },
    )


def check_frozen_sources(component_dir: Path) -> None:
    points_path = Path(
        "figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_points.tsv"
    )
    stats_path = Path(
        "figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_stats.json"
    )
    tissues_path = Path(
        "figures/ga_gpt2_vector_v6/registry_tissue_examples.tsv"
    )
    points = pd.read_csv(points_path, sep="\t")
    stats = json.loads(stats_path.read_text())
    q4.require(len(points) == stats["retained_cells"] == 11_472, "UMAP sample changed")
    union = (points["expr_67364"] > 0) | (points["expr_97630"] > 0)
    both = (points["expr_67364"] > 0) & (points["expr_97630"] > 0)
    q4.require(int(union.sum()) == 8_432, "Positive-cell union changed")
    q4.require(int(both.sum()) == 976, "Double-positive count changed")
    q4.require(
        int((points["sample_reason"] == "double_zero_background").sum()) == 3_040,
        "Double-zero background sample changed",
    )
    q4.require(
        stats["ENST00000367364"]["display_cap"] == 5
        and stats["ENST00000697630"]["display_cap"] == 3,
        "Isoform display caps changed",
    )

    tissues = pd.read_csv(tissues_path, sep="\t")
    expected_tissues = [
        "Marrow",
        "Myeloma",
        "Frontal cortex",
        "Retina",
        "Heart",
        "Prostate",
        "Glioma",
        "Ovary",
        "Brain",
        "Fibroblast→HSC",
    ]
    q4.require(
        tissues["display"].tolist() == expected_tissues,
        f"Registry tissue examples changed: {tissues['display'].tolist()}",
    )
    registry = pd.read_csv(REGISTRY_FROZEN, sep="\t", dtype=str)
    cells = pd.to_numeric(registry["isoquant_cells"], errors="coerce").fillna(0)
    done = registry["isoquant_status"].fillna("").str.lower().eq("done")
    q4.require(int(done.sum()) == 469, "Frozen completed-run count changed")
    q4.require(int(cells[done].sum()) == 845_781, "Frozen cell count changed")
    q4.require(
        int(registry.loc[done, "gse"].nunique()) == 31,
        "Frozen completed-study count changed",
    )

    inventory = pd.read_csv(
        component_dir / "three_axis_inventory_source.tsv", sep="\t"
    ).set_index("axis")
    expected_inventory = {
        "DIU": (2_396, 10_494),
        "APA": (4_168, 13_214),
        "ASE": (0, 11_506),
    }
    for axis, (passed, tested) in expected_inventory.items():
        q4.require(
            int(inventory.loc[axis, "n_sig_fdr05_effect_gate"]) == passed
            and int(inventory.loc[axis, "n_genes_tested"]) == tested,
            f"{axis} inventory changed",
        )
    ase = pd.read_csv(ASE_FULL, sep="\t")
    q4.require(len(ase) == 11_506, "ASE tested-gene count changed")
    q4.require(int(ase["sig"].sum()) == 0, "ASE adjusted significant count changed")
    q4.require(int((ase["pval"] < 0.05).sum()) == 673, "ASE raw-P count changed")
    q4.require(int((ase["qval"] < 0.05).sum()) == 0, "ASE q<0.05 count changed")
    q4.require(
        abs(float(ase["qval"].min()) - 0.6385836385836386) < 1e-12,
        "ASE minimum q changed",
    )
    print(
        "BIOLOGICAL SOURCE PASS",
        {
            "embedding_cells": stats["embedding_cells"],
            "displayed_cells": len(points),
            "positive_either_isoform": int(union.sum()),
            "registry_examples": expected_tissues,
            "inventory": expected_inventory,
            "ASE_raw_p_lt_0.05": 673,
            "ASE_min_q": float(ase["qval"].min()),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v6"),
    )
    parser.add_argument(
        "--skeleton",
        type=Path,
        default=Path("figures/ga_gpt2_vector_v6/chatgpt_skeleton_compact.svg"),
    )
    parser.add_argument(
        "--component-dir",
        type=Path,
        default=Path("figures/ga_gpt2_components_v2"),
    )
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()

    source_qa.check_catalog(args.component_dir)
    check_hashes("V5 FREEZE", V5_FROZEN_HASHES)
    check_hashes("RAW API", RAW_API_HASHES)
    check_hashes("DERIVED SOURCE", DERIVED_HASHES)
    check_compact_skeleton(args.skeleton)
    check_frozen_sources(args.component_dir)
    q4.REQUIRED_TEXT = REQUIRED_TEXT
    q4.FORBIDDEN_VISUAL_TERMS = FORBIDDEN_VISUAL_TERMS
    q4.RETIRED_PHRASES = RETIRED_PHRASES
    exports = q4.check_exports(args.stem, args.skeleton, args.dpi)

    provenance_files = [
        *exports,
        args.skeleton,
        Path("figures/scripts/audit_ga_v6_sources.py"),
        Path("figures/scripts/prepare_ga_v6_sources.py"),
        Path("figures/scripts/prepare_ga_v6_skeleton.py"),
        Path("figures/scripts/render_ga_gpt2_vector_v6.py"),
        Path("figures/scripts/compose_ga_gpt2_vector_v3.py"),
        Path("figures/scripts/normalize_ga_vector_v6.py"),
        Path(__file__),
        Path("docs/GRAPHICAL_ABSTRACT_VECTOR_V6_CONTRACT.md"),
        args.component_dir / "catalog_system_composition_20260726.tsv",
        args.component_dir / "ptprc_two_isoform_switch.tsv",
        args.component_dir / "three_axis_inventory_source.tsv",
        REGISTRY_FROZEN,
        ASE_FULL,
        *RAW_API_HASHES,
        *DERIVED_HASHES,
    ]
    print("SHA256")
    for path in provenance_files:
        q4.require(path.is_file(), f"Missing provenance file: {path}")
        print(sha256(path), path)
    print("ALL V6 QA CHECKS PASS")


if __name__ == "__main__":
    main()
