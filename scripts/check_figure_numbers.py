#!/usr/bin/env python3
"""Fail if a delivered figure prints a superseded denominator.

The manuscript text has been guarded against these numbers for a while, but
nothing looked inside the figures. Supplementary Figure 3 shipped a screenshot
reading "469 IsoQuant-complete runs from 31 studies and 845,781 called cells"
while its own footnote claimed the frozen 453 / 34 / 923,389 release, and
Supplementary Figure 1 printed per-study cell counts summing to 452,213. Both
passed every check in the build.

Only what a reader can actually see counts. For SVGs that means the text
nodes: Figure 1 embeds raster assets as base64, and "469" turns up inside that
payload in every build.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "docs/Figure_Script_mapping.tsv"

# Superseded denominators and the artefacts that produced them. A figure that
# prints any of these is describing a release the paper does not report.
DEAD = {
    "469": "run count of the candidate manifest",
    "845,781": "cell count of the candidate manifest",
    "850,938": "cell count of the 434-run candidate",
    "434 run": "run count of the 434-run candidate",
    "31 studies": "study count of the candidate manifest",
    "30 studies": "study count of the 434-run candidate",
    "452,213": "per-study cell sum from the superseded cohort table",
    "10,494": "run-as-donor isoform-usage gene universe",
    "13,214": "run-as-donor poly(A)-site gene universe",
    "1,045,231": "cells summed over run rows instead of per study",
    "candidate_20260801": "superseded release identifier",
    "candidate_20260727": "superseded release identifier",
    "IsoQuant-complete": "wording of the candidate manifest download",
    # PTPRC statistics from the run-as-donor analysis, superseded by the
    # biological-unit rerun. These shipped inside the SF3 gene-card tile for a
    # week: the figure read q 0.00563 where the paper reports q = 0.000879.
    "q 0.00563": "run-as-donor PTPRC isoform-usage q-value",
    "q 0.00393": "run-as-donor PTPRC poly(A)-site q-value",
    "q 0.955": "run-as-donor PTPRC allelic-interaction q-value",
    "q 0.00493": "older still PTPRC isoform-usage q-value",
    "q 0.00348": "older still PTPRC poly(A)-site q-value",
}


def visible_text(path: Path) -> str:
    if path.suffix.lower() == ".pdf":
        out = subprocess.run(
            ["pdftotext", str(path), "-"], capture_output=True, text=True
        )
        return out.stdout
    raw = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix.lower() == ".svg":
        # text nodes only; embedded base64 raster data is not something a
        # reader sees, and it contains arbitrary digit runs
        return "\n".join(re.findall(r">([^<>]*)<", raw))
    return raw


def screenshot_text() -> list[tuple[str, str]]:
    """What the portal screenshots in Supplementary Figure 3 actually say.

    pdftotext cannot see this. The tiles are raster screenshots embedded in the
    PDF, so scanning the PDF's text objects returns nothing from them - which is
    exactly how the 469-run Download tile survived every check. The capture
    script records each page's rendered body text beside the PNG it wrote, and
    that record is the only machine-readable account of what the reader sees.
    """
    renderer = (ROOT / "figures/scripts/render_nar_sf_compact.py").read_text(
        encoding="utf-8"
    )
    used = sorted(set(re.findall(r'"([A-Za-z0-9_]+)/[^"]+\.png"', renderer)))
    out: list[tuple[str, str]] = []
    for directory in used:
        base = ROOT / "figures/website_walkthrough" / directory
        # the walkthrough metadata covers the six full pages; the two PTPRC
        # crops are captured by a different script and record themselves here
        for name, label in (
            ("portal_walkthrough_current_metadata.json", "pages"),
            ("ptprc_live_capture_v3.json", "PTPRC crops"),
        ):
            meta = base / name
            if meta.is_file():
                out.append((f"SF3 {label} ({directory})", meta.read_text(encoding="utf-8")))
    return out


def main() -> int:
    problems: list[str] = []
    print(f"{'figure':<28} {'file':<52} status")
    print("-" * 96)
    for row in csv.DictReader(MAPPING.open(newline=""), delimiter="\t"):
        if "superseded" in row["figure"]:
            continue
        path = ROOT / row["file"]
        if not path.is_file():
            problems.append(f"{row['figure']}: missing file {row['file']}")
            print(f"{row['figure']:<28} {row['file'][:51]:<52} MISSING")
            continue
        text = visible_text(path)
        hits = sorted(token for token in DEAD if token in text)
        if hits:
            for token in hits:
                problems.append(
                    f"{row['figure']} prints the dead number {token!r} "
                    f"({DEAD[token]}): {row['file']}"
                )
        print(f"{row['figure']:<28} {row['file'][:51]:<52} "
              f"{'PRINTS ' + ', '.join(hits) if hits else 'clean'}")

    shots = screenshot_text()
    if not shots:
        problems.append(
            "no capture metadata found for the screenshots the SF3 renderer uses; "
            "their contents cannot be checked"
        )
    for label, text in shots:
        hits = sorted(token for token in DEAD if token in text)
        for token in hits:
            problems.append(
                f"{label} show the dead number {token!r} ({DEAD[token]})"
            )
        print(f"{label:<28} {'capture metadata':<52} "
              f"{'SHOW ' + ', '.join(hits) if hits else 'clean'}")

    if problems:
        print("\nPROBLEMS")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print("\nNO DELIVERED FIGURE PRINTS A SUPERSEDED DENOMINATOR")
    return 0


if __name__ == "__main__":
    sys.exit(main())
