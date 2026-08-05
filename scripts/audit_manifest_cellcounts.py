#!/usr/bin/env python3
"""Recount every manifest run's called cells from its IsoQuant barcode file.

The registry convention is isoquant_cells == number of barcodes in
<run>.transcript_grouped_tag_CB_counts.barcodes.tsv. This script checks the
frozen manifest against the files on disk so a stale registry value cannot
survive into the manuscript denominators.
"""
from __future__ import annotations

import os

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong"))
MANIFEST = ROOT / "tables/release_candidate_20260801/release_candidate_manifest.tsv"
REGISTRY = BASE / "docs/sample_registry.tsv"
OUT = ROOT / "tables/modality_audit"


def barcode_count(run: str, isoquant_path: str | float) -> int | None:
    if not isinstance(isoquant_path, str) or not isoquant_path:
        return None
    for candidate in (
        BASE / isoquant_path / run / f"{run}.transcript_grouped_tag_CB_counts.barcodes.tsv",
        BASE / isoquant_path / f"{run}.transcript_grouped_tag_CB_counts.barcodes.tsv",
    ):
        if candidate.is_file():
            with candidate.open() as handle:
                return sum(1 for line in handle if line.strip())
    return None


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    registry = pd.read_csv(REGISTRY, sep="\t", dtype=str).set_index("srr")
    rows = []
    for _, record in manifest.iterrows():
        run = record["srr"]
        path = registry.at[run, "isoquant_path"] if run in registry.index else None
        observed = barcode_count(run, path)
        manifest_cells = pd.to_numeric(record["isoquant_cells"], errors="coerce")
        rows.append(
            {
                "srr": run,
                "gse": record["gse"],
                "manifest_cells": manifest_cells,
                "ondisk_barcodes": observed,
                "isoquant_path": path,
                "status": (
                    "no_barcode_file"
                    if observed is None
                    else "match"
                    if observed == manifest_cells
                    else "mismatch"
                ),
            }
        )
    frame = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / "manifest_cellcount_audit.tsv", sep="\t", index=False)

    counts = frame["status"].value_counts().to_dict()
    mismatched = frame[frame["status"] == "mismatch"]
    summary = {
        "manifest": str(MANIFEST.relative_to(ROOT)),
        "runs": int(len(frame)),
        "status_counts": counts,
        "manifest_total": int(frame["manifest_cells"].fillna(0).sum()),
        "recounted_total_where_file_exists": int(
            frame.loc[frame["ondisk_barcodes"].notna(), "ondisk_barcodes"].sum()
        ),
        "mismatches": mismatched.drop(columns=["status"]).to_dict("records"),
    }
    (OUT / "manifest_cellcount_audit.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str)[:4000])


if __name__ == "__main__":
    main()
