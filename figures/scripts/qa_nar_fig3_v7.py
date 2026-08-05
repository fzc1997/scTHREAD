#!/usr/bin/env python3
"""Independent source, statistics and export QA for corrected Figure 3 v7."""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

import qa_nar_fig3_v3 as Q3
import qa_nar_fig3_v5 as Q5
import qa_nar_fig3_v6 as Q6
import render_nar_fig3_v7 as V7


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
P0 = PROJECT / "tables" / "p0_biological_unit_rerun"
WIDTH_MM = 183.0
HEIGHT_MM = 180.0


def check_corrected_statistics(
    live_json: Path, observed_suffix: str = ""
) -> list[Path]:
    inputs = V7.load_corrected_inputs(
        observed_suffix,
        live_json if observed_suffix else None,
    )
    expected = (
        {
            "diu": (0.0008786102062975, 0.3382009946583691, True),
            "apa": (0.0006010844748858, 0.20355870762310796, True),
            "ase": (0.9487670103092783, 0.19381173597187495, False),
        }
        if observed_suffix
        else {
            "diu": (0.005629468025011201, 0.3382009946583691, True),
            "apa": (0.003928483764182047, 0.20355870762310796, True),
            "ase": (0.9553322261980175, 0.19381173597187495, False),
        }
    )
    for analysis, (qval, effect, significant) in expected.items():
        row = getattr(inputs, analysis)
        Q3.require(np.isclose(float(row.qval), qval), f"{analysis} q mismatch")
        Q3.require(np.isclose(float(row.effect), effect), f"{analysis} effect mismatch")
        Q3.require(bool(row.sig) is significant, f"{analysis} significance mismatch")

    live = json.loads(live_json.read_text())
    Q3.require(live["gene"]["gid"] == V7.V3.GID, "Live endpoint is not PTPRC")
    for analysis in expected:
        for field in ("pval", "effect", "qval"):
            Q3.require(
                np.isclose(
                    float(live["analyses"][analysis][field]),
                    float(inputs.api["analyses"][analysis][field]),
                ),
                f"Live {analysis}.{field} differs from corrected source",
            )
        Q3.require(
            bool(live["analyses"][analysis]["sig"])
            is bool(inputs.api["analyses"][analysis]["sig"]),
            f"Live {analysis}.sig differs from corrected source",
        )
    if observed_suffix:
        Q3.require(
            inputs.api["coverage"]["junctions"]
            == {"features": 8994, "molecules": 2341560, "studies": 12},
            "Live junction coverage was not propagated into the figure inputs",
        )

    tests = pd.read_csv(
        P0 / "study_stratified_replication" / "study_stratified_targeted_tests.tsv",
        sep="\t",
    )
    diu = tests[
        tests.analysis.eq("diu")
        & tests.gene_name.eq("PTPRC")
        & tests.eligible.astype(str).str.lower().eq("true")
    ].set_index("study")
    Q3.require(set(diu.index) == {"GSE276974", "GSE307660"}, "Missing DIU study")
    Q3.require(np.isclose(diu.loc["GSE276974", "pval"], 0.0099), "CCUS P mismatch")
    Q3.require(np.isclose(diu.loc["GSE307660", "pval"], 0.0001), "Myeloma P mismatch")
    Q3.require(
        diu.loc["GSE276974", "effect_equal_donor"] >= 0.20
        and diu.loc["GSE307660", "effect_equal_donor"] >= 0.20,
        "PTPRC DIU effect does not clear 0.20 in both studies",
    )
    concordance = pd.read_csv(
        P0 / "study_stratified_replication" / "cross_study_concordance.tsv",
        sep="\t",
    )
    row = concordance[
        concordance.analysis.eq("diu") & concordance.gene.eq("PTPRC")
    ]
    Q3.require(len(row) == 1, "Expected one PTPRC DIU concordance row")
    row = row.iloc[0]
    Q3.require(
        0.25 < row.spearman_signed_contrasts < 0.35,
        "Unexpected PTPRC signed-contrast concordance",
    )
    Q3.require(
        0.60 < row.sign_agreement_fraction < 0.63,
        "Unexpected PTPRC sign agreement",
    )
    print(
        "CORRECTED STATISTICS AND LIVE PARITY PASS",
        {
            "diu": expected["diu"],
            "apa": expected["apa"],
            "ase": expected["ase"],
            "study_p": diu.pval.to_dict(),
            "study_effect": diu.effect_equal_donor.to_dict(),
            "signed_contrast_spearman": row.spearman_signed_contrasts,
            "sign_agreement": row.sign_agreement_fraction,
        },
    )
    return [
        P0 / f"diu_observed{observed_suffix}.tsv",
        P0 / f"apa_observed{observed_suffix}.tsv",
        P0 / f"ase_observed{observed_suffix}.tsv",
        P0 / "validation_report.json",
        P0 / "study_stratified_replication" / "study_stratified_targeted_tests.tsv",
        P0 / "study_stratified_replication" / "study_stratified_profiles.tsv",
        P0 / "study_stratified_replication" / "cross_study_concordance.tsv",
        P0 / "study_stratified_replication" / "study_stratified_replication_manifest.json",
        live_json,
    ]


def check_exports(
    stem: Path,
    dpi: int,
    selected: dict[str, str],
    observed_suffix: str = "",
) -> list[Path]:
    outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
    for path in outputs:
        Q3.require(path.is_file() and path.stat().st_size > 0, f"Missing {path}")

    svg_root = ET.parse(outputs[1]).getroot()
    image_count = sum(node.tag.endswith("image") for node in svg_root.iter())
    text_count = sum(node.tag.endswith("text") for node in svg_root.iter())
    Q3.require(image_count == 4, f"Expected four portal rasters, found {image_count}")
    Q3.require(text_count >= 95, f"Too few editable SVG text nodes: {text_count}")

    pdf_text = Q3.command_output("pdftotext", str(outputs[0]), "-")
    q_tokens = (
        ("q = 0.00088", "q = 0.00060", "q = 0.95")
        if observed_suffix
        else ("q = 0.00563", "q = 0.00393", "q = 0.96")
    )
    required = (
        "PTPRC",
        "Usage fractions and study-stratified DTU",
        q_tokens[0],
        "effect = 0.338",
        q_tokens[1],
        "effect = 0.204",
        q_tokens[2],
        "effect = 0.194",
        "Targeted DIU: CCUS P=0.0099",
        "myeloma P=0.0001",
        "Signed cell-type contrasts:",
        "ρ=0.29",
        "sign agreement=62%",
        "Fractions use all 23 isoforms",
        "Junction evidence remains inspectable and exportable",
        "Gold arc is a display example; selecting any arc exposes",
        "SELECTED ARC · DISPLAY EXAMPLE",
        "198,639,341–198,692,347",
        selected["molecules"],
        selected["reads"],
        selected["runs"],
        selected["studies"],
        "row: live table · overview: frozen snapshot",
        "GET /api/gene/PTPRC/",
    )
    if observed_suffix:
        required += ("2,341,560 molecules",)
    for text in required:
        Q3.require(text in pdf_text, f"Required text missing: {text}")
    forbidden = (
        "q = 0.00493",
        "effect = 0.552",
        "q = 0.00348",
        "q = 1.00",
        "highly concordant",
    )
    if observed_suffix:
        forbidden += (
            "q = 0.00563",
            "q = 0.00393",
            "q = 0.96",
            "2,341,574 molecules",
        )
    for text in forbidden:
        Q3.require(text not in pdf_text, f"Legacy/overstated text remains: {text}")

    pdf_info = Q3.command_output("pdfinfo", str(outputs[0]))
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    Q3.require(match is not None, "Could not parse PDF page size")
    width_pt, height_pt = map(float, match.groups())
    Q3.require(abs(width_pt - WIDTH_MM / 25.4 * 72) < 0.1, "PDF width mismatch")
    Q3.require(abs(height_pt - HEIGHT_MM / 25.4 * 72) < 0.1, "PDF height mismatch")
    fonts = Q3.command_output("pdffonts", str(outputs[0]))
    Q3.require("ArialMT" in fonts and "Arial-BoldMT" in fonts, "Arial not embedded")
    Q3.require("Type 3" not in fonts, "Type 3 font found")

    expected_px = (
        int(WIDTH_MM / 25.4 * dpi),
        int(HEIGHT_MM / 25.4 * dpi),
    )
    with Image.open(outputs[2]) as image:
        Q3.require(image.size == expected_px, f"PNG size mismatch: {image.size}")
        Q3.require(image.mode == "RGBA", f"PNG mode mismatch: {image.mode}")
        image_dpi = image.info.get("dpi")
        Q3.require(
            image_dpi is not None and all(abs(value - dpi) < 0.1 for value in image_dpi),
            f"PNG DPI mismatch: {image_dpi}",
        )
    print(
        "V7 DELIVERY ASSERTIONS PASS",
        {
            "canvas_mm": [WIDTH_MM, HEIGHT_MM],
            "png_px": expected_px,
            "svg_images": image_count,
            "svg_editable_text": text_count,
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, default=FIGURES / "NAR_Fig3_v7")
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--live-json", type=Path, required=True)
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=V7.CORRECTED_SNAPSHOT,
        help="Renderer-generated PTPRC snapshot",
    )
    parser.add_argument(
        "--observed-suffix",
        default="",
        help="Optional suffix for observed tables, for example _9999",
    )
    args = parser.parse_args()

    corrected_sources = check_corrected_statistics(
        args.live_json, args.observed_suffix
    )
    portal_files, metadata_files = Q6.check_used_portal_sources()
    junction_files, selected = Q6.check_junction_table()
    Q5.check_centroid_table()
    outputs = check_exports(
        args.stem.resolve(),
        args.dpi,
        selected,
        args.observed_suffix,
    )

    run_script = Path(__file__).with_name(
        "run_nar_fig3_v8.slurm"
        if args.stem.name == "NAR_Fig3_v8"
        else "run_nar_fig3_v7.slurm"
    )
    provenance = [
        *outputs,
        *corrected_sources,
        *portal_files,
        *metadata_files,
        *junction_files,
        Q5.CENTROID_TABLE,
        args.snapshot,
        Path(__file__).with_name("render_nar_fig3_v6.py"),
        Path(__file__).with_name("render_nar_fig3_v7.py"),
        Path(__file__),
        run_script,
    ]
    print("SHA256")
    for path in provenance:
        Q3.require(path.is_file(), f"Missing provenance file: {path}")
        print(Q3.sha256(path), path)
    report = {
        "status": "PASS",
        "figure": str(outputs[0]),
        "pdf_sha256": Q3.sha256(outputs[0]),
        "svg_sha256": Q3.sha256(outputs[1]),
        "png_sha256": Q3.sha256(outputs[2]),
        "live_json_sha256": Q3.sha256(args.live_json),
    }
    report_path = FIGURES / "_qa" / f"{args.stem.name}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print("ALL NAR FIGURE 3 V7 QA CHECKS PASS")


if __name__ == "__main__":
    main()
