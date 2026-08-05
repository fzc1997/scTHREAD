#!/usr/bin/env python3
"""Check corrected NAR supplementary figures against frozen source invariants."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def pdf_text(path: Path) -> str:
    return subprocess.run(
        ["pdftotext", str(path), "-"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figdir", type=Path, default=Path("figures"))
    parser.add_argument(
        "--capture-metadata",
        type=Path,
        default=Path(
            "figures/website_walkthrough/corrected_20260727/"
            "ptprc_live_capture_v3.json"
        ),
    )
    args = parser.parse_args()

    expected = {
        1: ("469", "845,781", "human", "mouse"),
        5: ("1,971", "8,092", "2,527", "10,531", "6,930", "19 donors"),
        8: ("Corrected gene card",),
        9: ("MS4A1", "q=0.0503", "eff=0.18", "q=0.0039", "eff=0.34", "not tested"),
    }
    forbidden = {
        1: ("unknown",),
        5: ("2,396", "10,494", "4,168", "13,214", "11,506", "60 genes"),
        9: ("q=0.0049", "eff=0.55"),
    }
    for number, tokens in expected.items():
        pdf = args.figdir / f"NAR_SF{number}.pdf"
        require(pdf.is_file() and pdf.stat().st_size > 0, f"Missing {pdf}")
        text = pdf_text(pdf)
        for token in tokens:
            require(token in text, f"SF{number} missing {token!r}")
        for token in forbidden.get(number, ()):
            require(token not in text, f"SF{number} retains legacy token {token!r}")

    metadata = json.loads(args.capture_metadata.read_text(encoding="utf-8"))
    location = metadata["browser_location"]
    evidence = metadata["evidence_text"]
    require("loading" not in location.lower(), "Corrected capture was taken while loading")
    for token in (
        "effect 0.338",
        "q 0.00563",
        "effect 0.204",
        "q 0.00393",
        "effect 0.194",
        "q 0.955",
    ):
        require(token in evidence, f"Capture metadata missing {token!r}")

    walkthrough_path = (
        args.capture_metadata.parent / "portal_walkthrough_current_metadata.json"
    )
    walkthrough = json.loads(walkthrough_path.read_text(encoding="utf-8"))
    require(len(walkthrough["pages"]) == 6, "Walkthrough must contain six current pages")
    require(
        all(item["expected_text_present"] for item in walkthrough["pages"]),
        "A current walkthrough page failed its content assertion",
    )
    download = next(
        item
        for item in walkthrough["pages"]
        if item["file"] == "05_download_current.png"
    )
    require(
        "Frozen 469-run release TSV" in download["body_text_excerpt"],
        "Download capture does not expose the frozen manuscript release",
    )

    print("CORRECTED SUPPLEMENTARY QA PASS")


if __name__ == "__main__":
    main()
