#!/usr/bin/env python3
"""Independent delivery checks for the refreshed NAR supplementary figures."""

from __future__ import annotations

import os

import re
import subprocess
from pathlib import Path

import pandas as pd
from PIL import Image


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
REGISTRY = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def command(*args: str) -> str:
    return subprocess.run(
        args,
        check=True,
        text=True,
        capture_output=True,
    ).stdout


def pdf_text(path: Path) -> str:
    return command("pdftotext", str(path), "-")


def main() -> None:
    registry = pd.read_csv(REGISTRY, sep="\t", dtype=str)
    done = registry.loc[
        registry["isoquant_status"].fillna("").str.lower().eq("done")
    ].copy()
    done["isoquant_cells"] = pd.to_numeric(
        done["isoquant_cells"],
        errors="coerce",
    )
    require(len(done) == 469, "Frozen registry must contain 469 membership records")
    require(done["gse"].nunique() == 31, "Frozen registry must contain 31 studies")
    require(
        int(done["isoquant_cells"].fillna(0).sum()) == 845_781,
        "Frozen registry must sum to 845,781 isoquant_cells",
    )

    for number in range(1, 11):
        pdf = FIGURES / f"NAR_SF{number}.pdf"
        png = FIGURES / f"NAR_SF{number}.png"
        require(pdf.is_file() and pdf.stat().st_size > 10_000, f"Missing {pdf}")
        require(png.is_file() and png.stat().st_size > 10_000, f"Missing {png}")
        require("Pages:           1" in command("pdfinfo", str(pdf)), f"{pdf} not one page")
        require("Type 3" not in command("pdffonts", str(pdf)), f"Type 3 font in {pdf}")
        with Image.open(png) as image:
            require(image.width >= 2_000, f"PNG too narrow: {png}")
            require(image.mode == "RGBA", f"PNG lacks transparency: {png}")

    sf1_text = pdf_text(FIGURES / "NAR_SF1.pdf")
    for phrase in (
        "469-run membership snapshot",
        "845,781",
        "Frozen-release runs by species",
        "Frozen-release isoquant_cells by species",
    ):
        require(phrase in sf1_text, f"SF1 text missing: {phrase}")
    for stale in (
        "844,864",
        "818,202",
        "460 IsoQuant-complete",
        "469 IsoQuant-complete",
    ):
        require(stale not in sf1_text, f"SF1 retains stale value: {stale}")

    sf10 = FIGURES / "NAR_SF10.pdf"
    sf10_text = pdf_text(sf10)
    for phrase in (
        "Cross-species utility",
        "previously published mouse gastrulation dataset",
        "CRA044500",
        "Mouse scONT query context",
        "Descriptive Malat1 usage contrast at E6.5",
        "no P/FDR",
    ):
        require(phrase in sf10_text, f"SF10 text missing: {phrase}")
    for stale in (
        "Catalog composition by species and platform",
        "Cell-type-specific Malat1 DTU at E6.5",
        "gene-level FDR",
        "DTU TEST",
        "Portal-to-inference traceability",
        "internally generated",
    ):
        require(stale not in sf10_text, f"SF10 retains stale wording: {stale}")
    page = command("pdfinfo", str(sf10))
    match = re.search(r"Page size:\s+([0-9.]+) x ([0-9.]+) pts", page)
    require(match is not None, "SF10 page size missing")
    width_mm = float(match.group(1)) * 25.4 / 72
    height_mm = float(match.group(2)) * 25.4 / 72
    require(abs(width_mm - 183.0) < 0.15, f"Unexpected SF10 width: {width_mm}")
    require(abs(height_mm - 190.0) < 0.15, f"Unexpected SF10 height: {height_mm}")

    print(
        "ALL NAR SUPPLEMENTARY REFRESH QA CHECKS PASS",
        {
            "frozen_release": [469, 31, 845_781],
            "figures": 10,
            "sf10_mm": [round(width_mm, 2), round(height_mm, 2)],
        },
    )


if __name__ == "__main__":
    main()
