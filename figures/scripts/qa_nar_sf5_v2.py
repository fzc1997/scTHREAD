#!/usr/bin/env python3
"""Validate 9,999-permutation NAR Supplementary Figure 5 v2."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
P0 = ROOT / "tables/p0_biological_unit_rerun"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=ROOT / "figures/NAR_SF5_v2",
    )
    args = parser.parse_args()
    pdf = args.stem.with_suffix(".pdf")
    svg = args.stem.with_suffix(".svg")
    png = args.stem.with_suffix(".png")
    for path in (pdf, svg, png):
        require(path.is_file() and path.stat().st_size > 0, f"missing {path}")
    text = subprocess.run(
        ["pdftotext", str(pdf), "-"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    for axis in ("ASE", "DIU", "APA"):
        table = pd.read_csv(P0 / f"{axis.lower()}_observed_9999.tsv", sep="\t")
        token = f"{int(table.sig.astype(bool).sum()):,}/{len(table):,}"
        require(token in text, f"{axis} count missing: {token}")
    for token in ("1,971/8,092", "2,527/10,531"):
        require(token not in text, f"legacy count remains: {token}")
    require(
        "Negative-control diagnostics" in text,
        "negative-control interpretation is missing",
    )
    print("NAR SF5 V2 QA PASS")


if __name__ == "__main__":
    main()
