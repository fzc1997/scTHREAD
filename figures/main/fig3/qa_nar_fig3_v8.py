#!/usr/bin/env python3
"""Independent QA for NAR Figure 3 v9 (CD45 variable-exon evidence).

Checks three things the renderer cannot check for itself:

1. the exon-inclusion measurement is internally consistent and matches what the
   figure prints, recomputed from the per-run tables rather than trusted;
2. the delivered PDF/SVG/PNG carry the new claims and none of the superseded
   ones — in particular no trace of the two truncation-sink transcripts that
   v8 panels c/d/e were built on;
3. the shared portal sources and gene-level statistics still validate.
"""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

import qa_nar_fig3_v3 as Q3
import qa_nar_fig3_v5 as Q5
import qa_nar_fig3_v6 as Q6
import qa_nar_fig3_v7 as Q7
import render_nar_fig3_v8 as V9


FIGURES = Path(__file__).resolve().parents[1]
WIDTH_MM = V9.WIDTH_MM
HEIGHT_MM = V9.HEIGHT_MM

# Superseded: the two most-abundant models that v8 c/d/e were built on. A
# 3-exon fragment and a retained_intron model, both truncation sinks.
RETIRED_TRANSCRIPTS = ("ENST00000367364", "ENST00000697630")


def check_exon_tables() -> list[Path]:
    """Recompute the figure's numbers from the per-run counts."""
    sources = [
        V9.PROJECT / "tables" / "cd45_exon_inclusion_GSE307660.tsv",
        V9.PROJECT / "tables" / "cd45_exon_inclusion_GSE276974.tsv",
        V9.ACCEPTOR_TABLE,
        V9.PER_RUN_TABLE,
    ]
    for path in sources:
        Q3.require(path.is_file(), f"Missing CD45 source: {path}")

    # Recompute through the builder's own functions rather than a private copy
    # of the logic, so a change to the thresholds or the alias map cannot pass
    # QA by being duplicated identically in two places.
    import sys
    sys.path.insert(0, str(V9.PROJECT / "scripts"))
    import build_cd45_figdata as B

    raw = pd.concat([pd.read_csv(p, sep="\t") for p in sources[:2]])
    frame = V9.load_acceptors()
    recomputed_frame = B.load_counts(list(sources[:2]))
    strata = B.surviving_strata(recomputed_frame)
    kept = recomputed_frame.merge(
        strata[["study", "run", "cell_type"]], on=["study", "run", "cell_type"]
    )

    # Every displayed fraction must be reproducible from the raw read counts.
    for row in frame.itertuples():
        group = kept[kept.study.eq(row.study) & kept.cell_type.eq(row.cell_type)]
        total = group.reads.sum()
        got = group[group.acceptor.eq(row.acceptor)].reads.sum() / total
        Q3.require(
            np.isclose(got, row.frac, atol=1e-9),
            f"{row.study}/{row.cell_type}/{row.acceptor}: table says "
            f"{row.frac:.6f}, recomputed {got:.6f}",
        )
        Q3.require(
            int(total) == int(row.total_reads),
            f"{row.study}/{row.cell_type}: total {row.total_reads} != {total}",
        )

    per_run = pd.read_csv(V9.PER_RUN_TABLE, sep="\t")
    Q3.require(
        set(per_run.cell_type) <= set(V9.CELL_TYPE_ORDER),
        "Per-run table carries a cell type the figure does not show",
    )
    Q3.require(
        (per_run.n >= B.MIN_STRATUM_READS).all(),
        "A per-run stratum is below the stratum threshold",
    )

    # The legend says "the four cell types carrying at least N reads in at least
    # K runs of both studies". That is an exhaustive claim, so it is not enough
    # to confirm the shown ones qualify -- no excluded cell type may qualify
    # either. Recompute the selection from the raw counts and demand equality.
    thresholds = V9.THRESHOLDS
    Q3.require(
        (thresholds["min_stratum_reads"], thresholds["min_cell_type_reads"],
         thresholds["min_cell_type_runs"])
        == (B.MIN_STRATUM_READS, B.MIN_CELL_TYPE_READS, B.MIN_CELL_TYPE_RUNS),
        "The manifest thresholds disagree with build_cd45_figdata.py",
    )
    selected, audit = B.select_cell_types(strata, sorted(set(recomputed_frame.study)))
    Q3.require(
        set(selected) == set(V9.CELL_TYPE_ORDER),
        f"Selection drifted: figure shows {sorted(V9.CELL_TYPE_ORDER)}, "
        f"thresholds select {sorted(selected)}",
    )

    # Every printed read/run count must come from the same filtered set as the
    # per-run strata, so the two accountings on the figure reconcile.
    per_run_check = pd.read_csv(V9.PER_RUN_TABLE, sep="\t")
    for cell_type in V9.CELL_TYPE_ORDER:
        for study in V9.STUDIES:
            rows = frame[frame.study.eq(study) & frame.cell_type.eq(cell_type)]
            Q3.require(len(rows) > 0, f"{cell_type} missing from {study}")
            group = per_run_check[
                per_run_check.study.eq(study) & per_run_check.cell_type.eq(cell_type)
            ]
            Q3.require(
                int(rows.total_reads.iloc[0]) == int(group.n.sum()),
                f"{cell_type}/{study}: printed reads {int(rows.total_reads.iloc[0])} "
                f"!= per-run total {int(group.n.sum())}",
            )
            Q3.require(
                int(rows.runs.iloc[0]) == int(group.run.nunique()),
                f"{cell_type}/{study}: printed runs {int(rows.runs.iloc[0])} "
                f"!= per-run runs {int(group.run.nunique())}",
            )
            Q3.require(
                int(rows.total_reads.iloc[0]) >= B.MIN_CELL_TYPE_READS,
                f"{cell_type}/{study} is below the read threshold",
            )
            Q3.require(
                int(rows.runs.iloc[0]) >= B.MIN_CELL_TYPE_RUNS,
                f"{cell_type}/{study} is below the run threshold",
            )

    # The headline contrast must actually hold in both studies.
    for study in V9.STUDIES:
        b_skip = frame[
            frame.study.eq(study)
            & frame.cell_type.eq("B cell")
            & frame.acceptor.eq("exon 7 (skip all)")
        ].frac
        mono_skip = frame[
            frame.study.eq(study)
            & frame.cell_type.eq("Monocyte")
            & frame.acceptor.eq("exon 7 (skip all)")
        ].frac
        b_value = float(b_skip.iloc[0]) if len(b_skip) else 0.0
        Q3.require(
            float(mono_skip.iloc[0]) > b_value * 3,
            f"{study}: monocyte skip-all is not clearly above B cell",
        )

    totals = frame.groupby(["study", "cell_type"]).total_reads.first().sum()
    print(
        "CD45 EXON TABLE ASSERTIONS PASS",
        {
            "informative_reads": int(totals),
            "runs": int(raw.run.nunique()),
            "studies": sorted(set(raw.study)),
            "per_run_strata": len(per_run),
            "selected_by_threshold": sorted(selected),
            "excluded": sorted(r["cell_type"] for r in audit if not r["selected"]),
        },
    )
    return sources + [V9.INCLUSION_MANIFEST]


def check_exports(stem: Path, dpi: int, selected: dict[str, str]) -> list[Path]:
    outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
    for path in outputs:
        Q3.require(path.is_file() and path.stat().st_size > 0, f"Missing {path}")

    svg_root = ET.parse(outputs[1]).getroot()
    image_count = sum(node.tag.endswith("image") for node in svg_root.iter())
    text_count = sum(node.tag.endswith("text") for node in svg_root.iter())
    # b keeps its portal raster and a keeps the gene card; c/d are now vector.
    Q3.require(image_count == 2, f"Expected two portal rasters, found {image_count}")
    Q3.require(text_count >= 95, f"Too few editable SVG text nodes: {text_count}")

    pdf_text = Q3.command_output("pdftotext", str(outputs[0]), "-")
    required = (
        "PTPRC",
        "B cell: exon-3 splice choice",
        "Monocyte: exon-3 splice choice",
        "Variable-exon usage across cell types and studies",
        "CD45RA/RB/RC versus CD45RO",
        "198, 141 and 144 bp",
        "Reads, not molecules",
        "informative reads",
        "skip all (RO)",
        "GSE307660",
        "GSE276974",
        "q = 0.00088",
        "effect = 0.338",
        "Junction evidence remains inspectable and exportable",
        selected["molecules"],
        "GET /api/gene/PTPRC/",
    )
    for text in required:
        Q3.require(text in pdf_text, f"Required text missing: {text}")

    forbidden = (
        *RETIRED_TRANSCRIPTS,
        "localization",
        "shared 0–5 molecule colour scale",
        "0–3 long-read molecules",
        "native color intensity",
        "Usage fractions and study-stratified DTU",
    )
    for text in forbidden:
        Q3.require(text not in pdf_text, f"Superseded text remains: {text}")

    pdf_info = Q3.command_output("pdfinfo", str(outputs[0]))
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", pdf_info)
    Q3.require(match is not None, "Could not parse PDF page size")
    width_pt, height_pt = map(float, match.groups())
    Q3.require(abs(width_pt - WIDTH_MM / 25.4 * 72) < 0.1, "PDF width mismatch")
    Q3.require(abs(height_pt - HEIGHT_MM / 25.4 * 72) < 0.1, "PDF height mismatch")
    fonts = Q3.command_output("pdffonts", str(outputs[0]))
    Q3.require("ArialMT" in fonts and "Arial-BoldMT" in fonts, "Arial not embedded")
    Q3.require("Type 3" not in fonts, "Type 3 font found")

    expected_px = (int(WIDTH_MM / 25.4 * dpi), int(HEIGHT_MM / 25.4 * dpi))
    with Image.open(outputs[2]) as image:
        Q3.require(image.size == expected_px, f"PNG size mismatch: {image.size}")
        Q3.require(image.mode == "RGBA", f"PNG mode mismatch: {image.mode}")
    print(
        "V9 DELIVERY ASSERTIONS PASS",
        {
            "canvas_mm": [WIDTH_MM, HEIGHT_MM],
            "png_px": expected_px,
            "svg_rasters": image_count,
            "svg_editable_text": text_count,
        },
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, default=FIGURES / "NAR_Fig3_v9")
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--live-json", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, default=V9.V7.CORRECTED_SNAPSHOT)
    parser.add_argument("--observed-suffix", default="")
    args = parser.parse_args()

    corrected = Q7.check_corrected_statistics(args.live_json, args.observed_suffix)
    exon_sources = check_exon_tables()
    junction_files, selected = Q6.check_junction_table()
    Q5.check_centroid_table()
    outputs = check_exports(args.stem.resolve(), args.dpi, selected)

    provenance = [
        *outputs,
        *corrected,
        *exon_sources,
        *junction_files,
        Q5.CENTROID_TABLE,
        args.snapshot,
        V9.PROJECT / "scripts" / "measure_cd45_exon_inclusion.py",
        Path(__file__).with_name("render_nar_fig3_v8.py"),
        Path(__file__),
    ]
    print("SHA256")
    for path in provenance:
        Q3.require(path.is_file(), f"Missing provenance file: {path}")
        print(Q3.sha256(path), path)
    print("ALL NAR FIGURE 3 V9 QA CHECKS PASS")


if __name__ == "__main__":
    main()
