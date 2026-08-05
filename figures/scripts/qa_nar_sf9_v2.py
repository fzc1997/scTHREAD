#!/usr/bin/env python3
"""Validate 9,999-permutation NAR Supplementary Figure 9 v2."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=ROOT / "figures/NAR_SF9_v2",
    )
    args = parser.parse_args()
    outputs = [args.stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
    for path in outputs:
        require(path.is_file() and path.stat().st_size > 0, f"missing {path}")
    text = subprocess.run(
        ["pdftotext", str(outputs[0]), "-"],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    for token in ("MS4A1", "q=0.0370", "q=0.0006", "eff=0.18", "eff=0.34"):
        require(token in text, f"missing {token}")
    for token in ("q=0.0503", "q=0.0039"):
        require(token not in text, f"legacy token remains: {token}")
    print("NAR SF9 V2 QA PASS")


if __name__ == "__main__":
    main()
