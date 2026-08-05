#!/usr/bin/env python3
"""Independent source and export QA for NAR_Fig3_mouse_scONT_v1."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image
import pandas as pd

import render_nar_fig3_mouse_scont as R


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def output(*args: str) -> str:
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def check_sources() -> list[Path]:
    study = pd.read_csv(R.STUDY_ATLAS, sep="\t")
    own = study.loc[study["gse"].eq("OWN_ASE_scONT")]
    require(len(own) == 1, "OWN_ASE_scONT row is not unique")
    require(int(own.iloc[0]["n_runs"]) == 6, "Run count drift")
    require(int(own.iloc[0]["n_cells"]) == 68417, "Registered cell count drift")

    comp = pd.read_csv(R.CELL_COMPOSITION, sep="\t")
    require(int(comp["Total"].sum()) == 25621, "Analysis cell denominator drift")
    require(len(comp) == 24, "Cell-type denominator drift")

    lineage = pd.read_csv(R.LINEAGE_DTU, sep="\t")
    require(
        dict(zip(lineage["celltype"], lineage["sig_genes"]))
        == {
            "Neuroectoderm": 732,
            "Presomitic_Mesoderm": 335,
            "Visceral_Endoderm_Yolk_Sac": 853,
        },
        "Lineage DTU counts drift",
    )

    temporal = pd.read_csv(R.TEMPORAL_CASES, sep="\t")
    tnrc6c = temporal.loc[temporal["gene_name"].eq("Tnrc6c")]
    require(len(tnrc6c) == 2, "Tnrc6c source rows drift")
    require(
        abs(float(tnrc6c["prop_diff_E65_E85"].abs().min()) - 0.996113433433989)
        < 1e-12,
        "Tnrc6c effect drift",
    )

    as_dtu = pd.read_csv(R.AS_DTU, sep="\t")
    require(len(as_dtu) == 5, "AS-DTU gene count drift")
    require((as_dtu["padj"] < 0.05).all(), "AS-DTU FDR guard failed")

    source_outputs = [
        R.TABLES / "Fig3_mouse_scONT_cell_composition.tsv",
        R.TABLES / "Fig3_mouse_scONT_lineage_DTU.tsv",
        R.TABLES / "Fig3_mouse_scONT_Tnrc6c_switch.tsv",
        R.TABLES / "Fig3_mouse_scONT_AS_DTU.tsv",
        R.TABLES / "Fig3_mouse_scONT_manifest.json",
    ]
    for path in source_outputs:
        require(path.is_file() and path.stat().st_size > 0, f"Missing source export: {path}")
    manifest = json.loads(source_outputs[-1].read_text(encoding="utf-8"))
    require(manifest["project_cohort"]["registered_cells"] == 68417, "Manifest cohort drift")
    require(manifest["analysis_subset"]["cells"] == 25621, "Manifest subset drift")
    require(len(manifest["inputs"]) == 7, "Manifest input count drift")
    for raw_path, metadata in manifest["inputs"].items():
        path = Path(raw_path)
        require(path.is_file(), f"Manifest input missing: {path}")
        require(R._sha256(path) == metadata["sha256"], f"Input hash mismatch: {path}")
    return source_outputs


def check_exports(stem: Path, dpi: int) -> list[Path]:
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    png = stem.with_suffix(".png")
    outputs = [pdf, svg, png]
    for path in outputs:
        require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    pdf_info = output("pdfinfo", str(pdf))
    require("Pages:           1" in pdf_info, "PDF must be one page")
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    require(match is not None, "Could not parse PDF dimensions")
    width_pt, height_pt = map(float, match.groups())
    require(abs(width_pt - R.WIDTH_MM / 25.4 * 72) < 0.1, "PDF width drift")
    require(abs(height_pt - R.HEIGHT_MM / 25.4 * 72) < 0.1, "PDF height drift")

    fonts = output("pdffonts", str(pdf))
    require("ArialMT" in fonts, "Arial regular is not embedded")
    require("Arial-BoldMT" in fonts, "Arial bold is not embedded")
    require("Type 3" not in fonts, "Type 3 fonts are forbidden")

    pdf_text = output("pdftotext", str(pdf), "-")
    normalized_text = re.sub(r"\s+", " ", pdf_text)
    required_text = (
        "OWN_ASE_scONT mouse gastrulation cohort",
        "6 runs · 68,417 registered project cells",
        "25,621-cell DTU analysis subset",
        "25,621 cells · 24 cell types",
        "Stage-associated isoform shifts across lineages",
        "Cell-bootstrap pseudo-replicates",
        "feature-level BH q < 0.05",
        "Genes with ≥1 significant isoform",
        "732 genes",
        "335 genes",
        "853 genes",
        "isoform–comparison hits",
        "Tnrc6c reciprocal isoform switch",
        "Δ usage = ±0.996",
        "Allelic-ratio heterogeneity at E6.5",
        "Gene-level max–min screen",
        "not isoform-level AS-DTU",
        "182,689 transcripts",
        "Alternative mouse use case",
    )
    for text in required_text:
        require(text in normalized_text, f"Required vector text missing: {text}")
    for forbidden in ("ENSG00000081237", "71,913 sampled cells"):
        require(forbidden not in normalized_text, f"Human Figure 3 content leaked: {forbidden}")

    root = ET.parse(svg).getroot()
    image_count = sum(1 for node in root.iter() if node.tag.endswith("image"))
    text_count = sum(1 for node in root.iter() if node.tag.endswith("text"))
    require(image_count == 0, f"Expected fully vector SVG, found {image_count} rasters")
    require(text_count >= 75, f"Too few editable SVG text objects: {text_count}")

    expected_px = (
        int(R.WIDTH_MM / 25.4 * dpi),
        int(R.HEIGHT_MM / 25.4 * dpi),
    )
    with Image.open(png) as image:
        require(image.size == expected_px, f"PNG dimensions drift: {image.size}")
        require(image.mode == "RGBA", f"PNG must be RGBA: {image.mode}")
        png_dpi = image.info.get("dpi")
        require(png_dpi is not None, "PNG DPI metadata missing")
        require(all(abs(value - dpi) < 0.1 for value in png_dpi), f"PNG DPI drift: {png_dpi}")

    print(pdf_info.strip())
    print(fonts.strip())
    print(
        "DELIVERY ASSERTIONS PASS",
        {
            "canvas_mm": [R.WIDTH_MM, R.HEIGHT_MM],
            "png_px": expected_px,
            "svg_rasters": image_count,
            "svg_editable_text": text_count,
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=R.FIGURES / "NAR_Fig3_mouse_scONT_v1",
    )
    parser.add_argument("--dpi", type=int, default=450)
    args = parser.parse_args()
    sources = check_sources()
    exports = check_exports(args.stem.resolve(), args.dpi)
    print("SHA256")
    for path in [
        *exports,
        *sources,
        Path(R.__file__).resolve(),
        Path(__file__).resolve(),
        Path(__file__).with_name("run_nar_fig3_mouse_scont.slurm"),
    ]:
        require(path.is_file(), f"Missing provenance file: {path}")
        print(R._sha256(path), path)
    print("ALL NAR FIGURE 3 MOUSE scONT QA CHECKS PASS")


if __name__ == "__main__":
    main()
