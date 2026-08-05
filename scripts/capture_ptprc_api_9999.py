#!/usr/bin/env python3
"""Capture and validate the live LAN PTPRC overview after 9,999 permutations."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from urllib.request import urlopen


EXPECTED = {
    "diu": (0.0008786102062975, 0.3382009946583691, True),
    "apa": (0.0006010844748858, 0.20355870762310796, True),
    "ase": (0.9487670103092783, 0.19381173597187495, False),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--url",
        default=(
            "http://10.168.3.4:4173/api/gene/PTPRC/"
            "overview?species=human"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "tables/PTPRC_api_overview_lan_9999_20260727.json",
    )
    args = parser.parse_args()
    with urlopen(args.url, timeout=30) as response:
        payload = json.load(response)
    if payload["gene"]["gid"] != "ENSG00000081237":
        raise RuntimeError("LAN endpoint did not return PTPRC")
    for axis, (qval, effect, significant) in EXPECTED.items():
        row = payload["analyses"][axis]
        if not math.isclose(float(row["qval"]), qval, rel_tol=1e-12):
            raise RuntimeError(f"{axis} q mismatch")
        if not math.isclose(float(row["effect"]), effect, rel_tol=1e-12):
            raise RuntimeError(f"{axis} effect mismatch")
        if bool(row["sig"]) is not significant:
            raise RuntimeError(f"{axis} significance mismatch")
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
