#!/usr/bin/env python3
"""Scientific, vector, geometry and delivery QA for graphical abstract v9."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

from lxml import etree
from PIL import Image

import compose_ga_gpt2_vector_v3 as composer
import prepare_ga_v9_sources as source_prep


ROOT = Path(__file__).resolve().parents[2]
V9_DATA = ROOT / "figures/ga_gpt2_vector_v9"
PTPRC_HASHES = {
    ROOT / "figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_points.tsv":
        "686360acf4550507fc8bc5bdc4c57568ef4b27d5697680d6572d9496b29b6364",
    ROOT / "figures/ga_gpt2_vector_v6/ptprc_two_isoform_umap_stats.json":
        "547f7d5f6d2fa5fb95312e4066bb6d15f1ab2ce67b91b39d892484c5fb48730c",
    ROOT / "figures/ga_gpt2_components_v2/ptprc_two_isoform_switch.tsv":
        "0e5c6352c6c6f7a95c581c3d7ae58f7edeb21f4c7219444e8d29f0688bf727d3",
}

REQUIRED_TEXT = (
    "scTHREAD: a single-cell long-read transcriptome atlas",
    "https://scthread.ai4sc.ac.cn",
    "Biological coverage",
    "Tissues · cell lines · disease · development · two species",
    "HUMAN TISSUES / CELLS",
    "MOUSE",
    "DEVELOPMENT",
    "453 runs · 34 datasets · 923,389 cells",
    "Human + mouse · ONT + PacBio",
    "Blood/immune",
    "Neural",
    "Cancer/models",
    "Endocrine",
    "Heart",
    "Development",
    "Registry-backed tissues & contexts · cells",
    "Frontal cortex",
    "CCUS marrow",
    "iNeuron trajectory",
    "Cross-platform BM mix",
    "Mouse gastrulation",
    "Myeloma marrow",
    "Splicing dynamics",
    "Mouse testis",
    "iPSC to cortical neuron",
    "Human ovary",
    "Unified database",
    "scTHREAD",
    ">200k isoforms",
    "cell-resolved",
    "study → run → cell → transcript",
    "Transcriptome evidence",
    "Gene expression",
    "across cell types",
    "Isoform usage",
    "cell-type shifts · DIU",
    "poly(A) sites",
    "alternative 3′ ends",
    "Allele-aware evidence",
    "ASE counts · junctions",
    "REF",
    "ALT",
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
    "Browse",
    "Search",
    "Download",
    "PTPRC isoform usage (% of all 23)",
    "Cell-type-associated genes",
    "FDR < 0.05 + effect-size criterion",
    "2,008 / 8,092 genes",
    "2,558 / 10,531 genes",
    "0 / 6,930 genes",
    "ASE: 538 nominal P < 0.05; 0 after FDR + effect filter (min q = 0.231)",
    "allele-aware expression",
)

FORBIDDEN_TERMS = (
    "benchmark",
    "method",
    "Genome-wide RNA events",
    "Cell-type landscape",
    "cell-resolved transcriptomes",
    "Four linked RNA evidence layers",
    "Allelic expression",
    "1,971/8,092",
    "2,527/10,531",
    "673 raw P<0.05",
    "min q=0.639",
    "0/11,506",
    "Genes with cell-type effects",
    "prespecified effect gate",
    "joint-gate hits",
)


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


def check_source_bundle() -> list[Path]:
    manifest_path = V9_DATA / "source_manifest_v9.json"
    require(manifest_path.is_file(), "V9 source manifest is missing")
    manifest = json.loads(manifest_path.read_text())
    require(
        manifest["catalog"]
        == {"runs": 453, "studies": 34, "cells": 923_389},
        f"Catalog manifest changed: {manifest['catalog']}",
    )
    require(
        manifest["systems"] == source_prep.EXPECTED_SYSTEMS,
        "System manifest changed",
    )
    require(
        manifest["inventory"] == source_prep.EXPECTED_INVENTORY,
        "Inventory manifest changed",
    )
    outputs: list[Path] = []
    for raw_path, record in manifest["outputs"].items():
        path = Path(raw_path)
        require(path.is_file(), f"Missing frozen v9 source: {path}")
        require(path.stat().st_size == record["bytes"], f"Size changed: {path}")
        require(
            source_prep.sha256(path) == record["sha256"],
            f"Hash changed: {path}",
        )
        outputs.append(path)
    for path, expected in PTPRC_HASHES.items():
        require(path.is_file(), f"Missing PTPRC source: {path}")
        require(source_prep.sha256(path) == expected, f"PTPRC source changed: {path}")
        outputs.append(path)

    composition = source_prep.pd.read_csv(
        V9_DATA / "catalog_system_composition_v9.tsv", sep="\t"
    )
    observed = {
        row.system: {
            "cells": int(row.cells),
            "samples": int(row.samples),
            "studies": int(row.studies),
        }
        for row in composition.itertuples(index=False)
    }
    require(observed == source_prep.EXPECTED_SYSTEMS, "Composition TSV changed")
    inventory = source_prep.pd.read_csv(
        V9_DATA / "three_axis_inventory_9999_v9.tsv", sep="\t"
    ).set_index("axis")
    for axis, expected in source_prep.EXPECTED_INVENTORY.items():
        require(
            int(inventory.loc[axis, "n_genes_tested"]) == expected["tested"],
            f"{axis} tested count changed",
        )
        require(
            int(inventory.loc[axis, "n_sig_fdr05_effect_gate"])
            == expected["passed"],
            f"{axis} passed count changed",
        )
        require(
            int(inventory.loc[axis, "n_raw_p_lt_0.05"]) == expected["raw"],
            f"{axis} raw-P count changed",
        )
        require(
            int(inventory.loc[axis, "permutations"]) == 9999,
            f"{axis} is not the 9,999-permutation result",
        )
    require(
        abs(float(inventory.loc["ASE", "min_q"]) - 0.231) < 1e-12,
        "ASE minimum q changed",
    )
    print(
        "V9 SOURCE BUNDLE PASS",
        {
            "catalog": manifest["catalog"],
            "systems": observed,
            "inventory": source_prep.EXPECTED_INVENTORY,
            "ASE_min_q": float(inventory.loc["ASE", "min_q"]),
        },
    )
    return [manifest_path, *outputs]


def check_skeleton(path: Path) -> None:
    root = composer.parse_svg(path)
    composer.validate_svg(root, require_groups=True)
    require(
        not [node for node in root.iter() if composer.local_name(node.tag) == "image"],
        "V9 skeleton contains a raster image",
    )
    require(
        not [node for node in root.iter() if composer.local_name(node.tag) == "text"],
        "V9 skeleton contains visible text",
    )
    ids = [node.get("id") for node in root.iter() if node.get("id")]
    require(len(ids) == len(set(ids)), "V9 skeleton has duplicate IDs")
    require("database-hierarchy" not in ids, "Old database hierarchy remains")
    title = root.xpath(".//*[@id='svgTitle']")
    desc = root.xpath(".//*[@id='svgDesc']")
    require(len(title) == len(desc) == 1, "Accessibility metadata is incomplete")
    require(
        "reconstructed from a PNG" not in "".join(title[0].itertext()),
        "Stale PNG reconstruction title remains",
    )
    require(
        "sequencing instruments" not in "".join(desc[0].itertext()),
        "Stale sequencing-instrument description remains",
    )
    embryo = root.xpath(".//*[@id='embryo-icon']")
    require(
        len(embryo) == 1
        and embryo[0].get("transform") == "translate(260 211) scale(0.48)",
        "Embryo pictogram placement changed",
    )
    print(
        "V9 SKELETON PASS",
        {
            "sha256": source_prep.sha256(path),
            "ids": len(ids),
            "database_hierarchy": 0,
            "embryo_transform": embryo[0].get("transform"),
        },
    )


def check_svg(path: Path) -> str:
    root = composer.parse_svg(path)
    composer.validate_svg(root, require_groups=True)
    ids = [node.get("id") for node in root.iter() if node.get("id")]
    require(len(ids) == len(set(ids)), "Final SVG has duplicate IDs")
    require("real-data-overlay" in ids, "Real-data overlay is missing")
    require("database-hierarchy" not in ids, "Old database hierarchy leaked")
    image_nodes = [
        node for node in root.iter() if composer.local_name(node.tag) == "image"
    ]
    text_nodes = [
        node for node in root.iter() if composer.local_name(node.tag) == "text"
    ]
    require(not image_nodes, f"Final SVG contains {len(image_nodes)} raster images")
    require(len(text_nodes) >= 100, f"Too few editable SVG text nodes: {len(text_nodes)}")

    observed_ids = set(ids)
    missing_references: list[str] = []
    for node in root.iter():
        for raw_name, value in node.attrib.items():
            name = composer.local_name(raw_name)
            for target in re.findall(r"url\\(#([^)]+)\\)", value):
                if target not in observed_ids:
                    missing_references.append(f"{name}:{target}")
            if name == "href" and value.startswith("#") and value[1:] not in observed_ids:
                missing_references.append(f"href:{value[1:]}")
    require(not missing_references, f"Broken SVG references: {missing_references[:10]}")
    svg_text = " ".join("".join(node.itertext()) for node in text_nodes)
    print(
        "V9 SVG STRUCTURE PASS",
        {
            "image_elements": 0,
            "editable_text_elements": len(text_nodes),
            "unique_ids": len(ids),
            "missing_references": 0,
        },
    )
    return svg_text


def parse_pdf_lines(pdf_path: Path) -> tuple[float, float, list[dict[str, object]]]:
    raw = command_output("pdftotext", "-bbox-layout", str(pdf_path), "-")
    parser = etree.XMLParser(resolve_entities=False, no_network=True)
    root = etree.fromstring(raw.encode("utf-8"), parser)
    page = root.xpath("//*[local-name()='page']")
    require(len(page) == 1, "Expected one bbox-layout page")
    width = float(page[0].get("width"))
    height = float(page[0].get("height"))
    lines: list[dict[str, object]] = []
    for line in root.xpath("//*[local-name()='line']"):
        words = ["".join(node.itertext()) for node in line if isinstance(node.tag, str)]
        text = " ".join(word for word in words if word)
        lines.append(
            {
                "text": text,
                "x_min": float(line.get("xMin")),
                "y_min": float(line.get("yMin")),
                "x_max": float(line.get("xMax")),
                "y_max": float(line.get("yMax")),
            }
        )
    return width, height, lines


def find_line(lines: list[dict[str, object]], phrase: str) -> dict[str, object]:
    matches = [line for line in lines if phrase in str(line["text"])]
    require(len(matches) == 1, f"Expected one PDF line for {phrase!r}: {matches}")
    return matches[0]


def intersects(first: dict[str, object], second: dict[str, object]) -> bool:
    return not (
        float(first["x_max"]) <= float(second["x_min"])
        or float(second["x_max"]) <= float(first["x_min"])
        or float(first["y_max"]) <= float(second["y_min"])
        or float(second["y_max"]) <= float(first["y_min"])
    )


def check_geometry(pdf_path: Path) -> None:
    width, height, lines = parse_pdf_lines(pdf_path)
    require(abs(width - 518.740157) < 0.01, f"Unexpected bbox width: {width}")
    require(abs(height - 221.102362) < 0.01, f"Unexpected bbox height: {height}")
    for line in lines:
        require(
            float(line["x_min"]) >= -0.02
            and float(line["y_min"]) >= -0.02
            and float(line["x_max"]) <= width + 0.02
            and float(line["y_max"]) <= height + 0.02,
            f"PDF text exceeds page: {line}",
        )

    title = find_line(lines, "scTHREAD: a single-cell long-read transcriptome atlas")
    require(float(title["y_min"]) >= 1.0, f"Page title lacks top safety: {title}")

    url = find_line(lines, "https://scthread.ai4sc.ac.cn")
    evidence = find_line(lines, "Transcriptome evidence")
    require(not intersects(url, evidence), "URL overlaps transcriptome-evidence heading")
    require(
        float(evidence["y_max"]) < 20.45,
        f"Evidence heading still covers the first card: {evidence}",
    )

    umap_title = find_line(lines, "PTPRC isoform expression")
    require(
        float(umap_title["y_max"]) < 54.10,
        f"PTPRC isoform title still intersects the UMAP frame: {umap_title}",
    )

    heatmap_label_names = ("B", "DC", "Ery", "Mo", "NK", "Pl", "Pro", "T")
    heatmap_labels: list[dict[str, object]] = []
    for label in heatmap_label_names:
        matches = [
            line
            for line in lines
            if str(line["text"]) == label
            and 100.0 < float(line["x_min"]) < 250.0
            and 195.0 < float(line["y_min"]) < 206.0
        ]
        require(
            len(matches) == 1,
            f"Heatmap label {label!r} not found exactly once: {matches}",
        )
        heatmap_labels.append(matches[0])
    heatmap_y_max = max(float(line["y_max"]) for line in heatmap_labels)
    require(
        max(float(line["y_min"]) for line in heatmap_labels)
        - min(float(line["y_min"]) for line in heatmap_labels)
        < 0.05,
        f"Heatmap labels are not aligned: {heatmap_labels}",
    )
    require(
        heatmap_y_max < 203.85,
        f"Heatmap labels exceed the card: {heatmap_labels}",
    )

    alt = find_line(lines, "ALT")
    require(float(alt["y_max"]) < 141.75, f"ALT label exceeds evidence card: {alt}")
    export = find_line(lines, "Export · CSV · JSON · API")
    require(
        float(export["y_max"]) < 192.80,
        f"Export label exceeds browser frame: {export}",
    )
    print(
        "V9 TARGETED GEOMETRY PASS",
        {
            "title_y_min": title["y_min"],
            "evidence_y_max": evidence["y_max"],
            "umap_title_y_max": umap_title["y_max"],
            "heatmap_labels_y_max": heatmap_y_max,
            "ALT_y_max": alt["y_max"],
            "export_y_max": export["y_max"],
        },
    )


def check_exports(stem: Path, skeleton: Path, dpi: int) -> list[Path]:
    svg_path = stem.with_suffix(".svg")
    pdf_path = stem.with_suffix(".pdf")
    png_path = stem.with_suffix(".png")
    paths = [svg_path, pdf_path, png_path]
    for path in paths:
        require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    svg_text = check_svg(svg_path)
    pdf_text = command_output("pdftotext", str(pdf_path), "-")
    combined = f"{svg_text}\n{pdf_text}"
    for required in REQUIRED_TEXT:
        require(required in combined, f"Required visible text is missing: {required}")
    for forbidden in FORBIDDEN_TERMS:
        require(
            re.search(rf"(?<!\\w){re.escape(forbidden)}(?!\\w)", combined, re.I)
            is None,
            f"Forbidden/stale visible text remains: {forbidden}",
        )

    pdf_info = command_output("pdfinfo", str(pdf_path))
    require("Pages:           1" in pdf_info, "PDF must contain one page")
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    require(match is not None, "Could not read PDF page size")
    width_pt, height_pt = map(float, match.groups())
    require(abs(width_pt - 518.74) < 0.2, f"Unexpected PDF width: {width_pt}")
    require(abs(height_pt - 221.102) < 0.2, f"Unexpected PDF height: {height_pt}")

    font_info = command_output("pdffonts", str(pdf_path))
    require("Type 3" not in font_info, "Type 3 fonts are not allowed")
    require(
        re.search(r"Arial|LiberationSans|Helvetica", font_info, re.I) is not None,
        "Publication-safe sans-serif font is missing",
    )
    image_info = command_output("pdfimages", "-list", str(pdf_path))
    require(
        re.search(r"(?m)^\s*\d+\s+\d+\s+", image_info) is None,
        "PDF contains a raster image",
    )

    with Image.open(png_path) as png:
        expected_size = (
            round(183.0 / 25.4 * dpi),
            round(78.0 / 25.4 * dpi),
        )
        require(png.size == expected_size, f"Unexpected PNG size: {png.size}")
        metadata_dpi = png.info.get("dpi")
        require(metadata_dpi is not None, "PNG DPI metadata is missing")
        require(
            all(abs(value - dpi) < 0.2 for value in metadata_dpi),
            f"PNG DPI is not {dpi}: {metadata_dpi}",
        )
    check_geometry(pdf_path)
    print("V9 EXPORT PASS", {"dpi": dpi, "png_size": expected_size})
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=Path("figures/scTHREAD_graphical_abstract_gpt2_hybrid_v9"),
    )
    parser.add_argument(
        "--skeleton",
        type=Path,
        default=V9_DATA / "chatgpt_skeleton_v9.svg",
    )
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()

    source_files = check_source_bundle()
    check_skeleton(args.skeleton)
    exports = check_exports(args.stem, args.skeleton, args.dpi)
    provenance = [
        *exports,
        args.skeleton,
        *source_files,
        ROOT / "figures/scripts/prepare_ga_v9_sources.py",
        ROOT / "figures/scripts/prepare_ga_v9_skeleton.py",
        ROOT / "figures/scripts/render_ga_gpt2_vector_v9.py",
        ROOT / "figures/scripts/compose_ga_gpt2_vector_v3.py",
        ROOT / "figures/scripts/normalize_ga_vector_v9.py",
        Path(__file__),
        ROOT / "docs/GRAPHICAL_ABSTRACT_VECTOR_V9_CONTRACT.md",
    ]
    print("SHA256")
    for path in provenance:
        require(path.is_file(), f"Missing provenance file: {path}")
        print(source_prep.sha256(path), path)
    print("ALL V9 QA CHECKS PASS")


if __name__ == "__main__":
    main()
