#!/usr/bin/env python3
"""Independent data, live-API and delivery QA for NAR Figure 3 v3."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


PROJECT = Path(__file__).resolve().parents[2]
SC_ROOT = PROJECT.parent
FIGURES = PROJECT / "figures"
ASSETS = FIGURES / "assets"
TABLES = PROJECT / "tables"
FIGDATA = SC_ROOT / "results" / "paper1" / "figdata"
F2DATA = SC_ROOT / "results" / "paper1" / "f2_grammar" / "figdata"

EXPECTED_CAPTURE_SHA256 = (
    "6ae052751cab44835d4d7680edd02ab96bba1680c74e6c5b1128bdc5a4467d04"
)
EXPECTED_CELL_TOTALS = {
    "B cell": 11_944,
    "Dendritic cell": 881,
    "Erythroid": 2_705,
    "Monocyte": 38_213,
    "NK": 22_412,
    "Plasma cell": 381,
    "Progenitor": 1_430,
    "T cell": 20_479,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(*args: str) -> str:
    return subprocess.run(
        args,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout


def fetch_json_without_proxy(url: str) -> dict:
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with opener.open(request, timeout=30) as response:
        require(response.status == 200, f"Live API returned HTTP {response.status}")
        return json.load(response)


def assert_close(observed: object, expected: object, label: str) -> None:
    require(
        np.isclose(float(observed), float(expected), rtol=0.0, atol=1e-12),
        f"{label} differs: {observed!r} versus {expected!r}",
    )


def selected_api_values(payload: dict) -> dict[str, object]:
    values: dict[str, object] = {
        "gene.gid": payload["gene"]["gid"],
        "gene.gname": payload["gene"]["gname"],
    }
    for analysis in ("diu", "apa", "ase"):
        for field in ("pval", "effect", "qval", "sig"):
            values[f"analyses.{analysis}.{field}"] = payload["analyses"][analysis][
                field
            ]
    for layer in ("isoforms", "pas", "junctions"):
        for field in ("features", "molecules", "studies"):
            values[f"coverage.{layer}.{field}"] = payload["coverage"][layer][field]
    return values


def compare_payloads(left: dict, right: dict, label: str) -> None:
    left_values = selected_api_values(left)
    right_values = selected_api_values(right)
    require(left_values.keys() == right_values.keys(), f"{label}: selected keys differ")
    for key in left_values:
        lvalue, rvalue = left_values[key], right_values[key]
        if isinstance(lvalue, bool) or isinstance(rvalue, bool):
            require(bool(lvalue) == bool(rvalue), f"{label}: {key} differs")
        elif isinstance(lvalue, (int, float)) and isinstance(rvalue, (int, float)):
            assert_close(lvalue, rvalue, f"{label}: {key}")
        else:
            require(lvalue == rvalue, f"{label}: {key} differs")


def check_data_and_live_api(
    live_url: str | None,
    live_json: Path | None,
) -> None:
    snapshot = json.loads((TABLES / "PTPRC_api_overview.json").read_text())
    require(snapshot["gene"]["gid"] == "ENSG00000081237", "Unexpected PTPRC gene ID")
    require(snapshot["gene"]["gname"] == "PTPRC", "Unexpected PTPRC gene symbol")

    source_specs = {
        "diu": (
            F2DATA / "diu_celltype.tsv",
            "gene",
            "ENSG00000081237",
            ("pval", "effect", "qval", "sig", "n_iso"),
        ),
        "apa": (
            F2DATA / "apa_celltype.tsv",
            "gene",
            "ENSG00000081237",
            ("pval", "effect", "qval", "sig", "n_pas"),
        ),
        "ase": (
            F2DATA / "ase_interaction.tsv",
            "gene",
            "PTPRC",
            ("pval", "effect", "qval", "sig"),
        ),
    }
    for analysis, (path, key, value, fields) in source_specs.items():
        table = pd.read_csv(path, sep="\t")
        rows = table.loc[table[key].astype(str).str.upper().eq(value.upper())]
        require(len(rows) == 1, f"Expected one {analysis} source row")
        row = rows.iloc[0]
        api_row = snapshot["analyses"][analysis]
        for field in fields:
            if field == "sig":
                require(
                    str(row[field]).lower() == str(api_row[field]).lower(),
                    f"{analysis}.{field} differs",
                )
            else:
                assert_close(row[field], api_row[field], f"{analysis}.{field}")

    usage = pd.read_csv(FIGDATA / "ptprc_isoform_usage.tsv", sep="\t")
    require(usage["transcript_id"].nunique() == 23, "Expected 23 PTPRC isoforms")
    require(usage["ct"].nunique() == 8, "Expected eight PTPRC cell types")
    totals = usage.groupby("ct")["count"].sum().round().astype(int).to_dict()
    require(totals == EXPECTED_CELL_TOTALS, f"Cell-type totals changed: {totals}")
    require(sum(totals.values()) == 98_445, "PTPRC matrix total changed")
    row_sums = usage.groupby("ct")["frac"].sum().to_numpy()
    require(
        np.allclose(row_sums, 1.0, rtol=0.0, atol=1e-12),
        "Usage rows are not normalized over all 23 isoforms",
    )
    expected_frac = usage["count"] / usage.groupby("ct")["count"].transform("sum")
    require(
        np.allclose(usage["frac"], expected_frac, rtol=0.0, atol=1e-12),
        "Usage fractions and molecule counts disagree",
    )
    top = usage.groupby("transcript_id")["count"].sum().nlargest(5).index
    top_coverage = (
        usage.loc[usage["transcript_id"].isin(top)].groupby("ct")["count"].sum()
        / usage.groupby("ct")["count"].sum()
    )
    require(top_coverage.min() > 0.82, "Top-five row coverage fell below 82%")
    require(top_coverage.max() < 0.89, "Top-five row coverage rose above 89%")

    capture_path = ASSETS / "ptprc_gene_card_live_v3.png"
    require(sha256(capture_path) == EXPECTED_CAPTURE_SHA256, "Portal crop changed")
    with Image.open(capture_path) as capture:
        require(capture.size == (2280, 580), f"Unexpected portal crop: {capture.size}")
    metadata = json.loads((ASSETS / "ptprc_live_capture_v3.json").read_text())
    metadata_text = f"{metadata['browser_location']}\n{metadata['evidence_text']}"
    for text in (
        "317 junctions",
        "8,994 junctions",
        "2,341,574 molecules",
        "Isoform usage (DIU) · effect 0.552 · q 0.00493",
        "Poly(A) usage (APA) · effect 0.201 · q 0.00348",
        "ASE cell-type interaction · effect 0.170 · q 1.00",
    ):
        require(text in metadata_text, f"Capture metadata lacks {text!r}")

    live: dict | None = None
    live_label: str | None = None
    if live_json:
        require(live_json.is_file(), f"Live API snapshot is missing: {live_json}")
        live = json.loads(live_json.read_text())
        live_label = str(live_json)
    elif live_url:
        live = fetch_json_without_proxy(live_url)
        live_label = live_url
    if live is not None:
        compare_payloads(snapshot, live, "frozen snapshot versus live API")
        require(live.get("species") == "human", "Live API species is not human")
        require(live["gene"].get("assembly") == "GRCh38", "Live API assembly changed")
        print("LIVE API PARITY PASS", live_label)
    print(
        "DATA ASSERTIONS PASS",
        {
            "matrix_molecules": sum(totals.values()),
            "cell_types": len(totals),
            "isoforms": usage["transcript_id"].nunique(),
            "top5_row_coverage": [
                round(float(top_coverage.min()), 4),
                round(float(top_coverage.max()), 4),
            ],
        },
    )


def check_exports(stem: Path, dpi: int) -> list[Path]:
    pdf_path = stem.with_suffix(".pdf")
    svg_path = stem.with_suffix(".svg")
    png_path = stem.with_suffix(".png")
    outputs = [pdf_path, svg_path, png_path]
    for path in outputs:
        require(path.is_file() and path.stat().st_size > 0, f"Missing export: {path}")

    svg_root = ET.parse(svg_path).getroot()
    image_count = sum(1 for node in svg_root.iter() if node.tag.endswith("image"))
    text_count = sum(1 for node in svg_root.iter() if node.tag.endswith("text"))
    require(image_count == 1, f"Expected one portal raster, found {image_count}")
    require(text_count >= 95, f"Too few editable SVG text elements: {text_count}")

    pdf_text = command_output("pdftotext", str(pdf_path), "-")
    for required in (
        "PTPRC",
        "ENSG00000081237",
        "8,994",
        "130,951",
        "2,341,574",
        "0.004933",
        "0.003478",
        "1.000000",
        "98,445",
        "values not renormalized",
        "tables = snapshot = live API",
    ):
        require(required in pdf_text, f"Required vector text is missing: {required}")
    for forbidden in ("lineage_isoform_DIU", '"export"', "RA / RO / RB", "GPT"):
        require(forbidden not in pdf_text, f"Legacy/pseudo content remains: {forbidden}")

    pdf_info = command_output("pdfinfo", str(pdf_path))
    require("Pages:           1" in pdf_info, "Figure PDF must contain one page")
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    require(match is not None, "Could not parse PDF page size")
    width_pt, height_pt = map(float, match.groups())
    require(abs(width_pt - 518.740) < 0.1, f"PDF width is not 183 mm: {width_pt}")
    require(abs(height_pt - 467.717) < 0.1, f"PDF height is not 165 mm: {height_pt}")

    font_info = command_output("pdffonts", str(pdf_path))
    require("ArialMT" in font_info, "Arial Regular is not embedded")
    require("Arial-BoldMT" in font_info, "Arial Bold is not embedded")
    require("Type 3" not in font_info, "Type 3 fonts are not allowed")

    expected_png = (
        round(183.0 / 25.4 * dpi),
        round(165.0 / 25.4 * dpi),
    )
    with Image.open(png_path) as image:
        require(image.size == expected_png, f"Unexpected PNG dimensions: {image.size}")
        require(image.mode == "RGBA", f"PNG is not RGBA: {image.mode}")
        image_dpi = image.info.get("dpi")
        require(image_dpi is not None, "PNG DPI metadata is missing")
        require(
            all(abs(value - dpi) < 0.1 for value in image_dpi),
            f"PNG DPI differs from {dpi}: {image_dpi}",
        )

    print("PDF INFO")
    print(pdf_info.strip())
    print("PDF FONTS")
    print(font_info.strip())
    print(
        "DELIVERY ASSERTIONS PASS",
        {
            "svg_images": image_count,
            "svg_editable_text": text_count,
            "png_px": expected_png,
            "dpi": dpi,
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v3",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument(
        "--live-url",
        default=None,
        help="Optional live overview endpoint for no-proxy parity checking.",
    )
    parser.add_argument(
        "--live-json",
        type=Path,
        default=None,
        help="Fresh live response captured on a network-reachable host.",
    )
    args = parser.parse_args()

    check_data_and_live_api(args.live_url, args.live_json)
    outputs = check_exports(args.stem.resolve(), args.dpi)
    print("SHA256")
    provenance = [
        *outputs,
        ASSETS / "ptprc_gene_card_live_v3.png",
        ASSETS / "ptprc_live_capture_v3.json",
        Path(__file__).with_name("render_nar_fig3_v3.py"),
        Path(__file__),
    ]
    if args.live_json:
        provenance.append(args.live_json.resolve())
    for path in provenance:
        require(path.is_file(), f"Missing provenance file: {path}")
        print(sha256(path), path)
    print("ALL NAR FIGURE 3 V3 QA CHECKS PASS")


if __name__ == "__main__":
    main()
