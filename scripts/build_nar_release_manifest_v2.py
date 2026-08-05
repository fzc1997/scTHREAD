#!/usr/bin/env python3
"""Build the release manifest under the single-cell long-read scope rule.

Scope is a property of the data, not of the archive it came from: every
single-cell long-read dataset belongs in the release whether it arrived through
GEO, GSA or a sequencing provider. That rule, and the authoritative per-study
cell count that goes with it, are decided upstream and recorded in
`authoritative_cellcount_clean_*.tsv`; this script does not re-derive either.

Two consequences shape the output:

* The study set comes from the authoritative table. Datasets ruled out there
  (genomic DNA, short-read-only) never enter, so no modality audit is repeated.
* Cell counts are authoritative **per study**, because the winning method is
  often an author-supplied total or a STARsolo call that has no per-run
  decomposition. Run rows still carry their own pipeline counts for provenance,
  and the two are reconciled and reported rather than silently summed.

Writes <output-dir>/release_manifest.tsv (one row per run),
<output-dir>/release_study_cells.tsv (one row per study, the cell authority)
and <output-dir>/release_summary.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

RELEASE_ID = "scTHREAD_NAR_release_20260803"
MANIFEST_COLUMNS = [
    "release", "srr", "gse", "species", "platform", "ena_assay",
    "biology_group", "description", "pipeline", "isoquant_status",
    "isoquant_cells", "isoquant_cells_method", "isoquant_transcripts",
    "annotation_status", "annotation_method", "metadata_resolution",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authoritative", type=Path, required=True,
                        help="upstream per-study cell-count table")
    parser.add_argument("--registry", type=Path, required=True,
                        help="live sample registry supplying run membership")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    authority = pd.read_csv(args.authoritative, sep="\t")
    for column in ("gse", "cell_count", "method"):
        if column not in authority.columns:
            raise ValueError(f"authoritative table lacks a {column} column")
    if authority["gse"].duplicated().any():
        raise ValueError("authoritative table must hold one row per study")

    # The registry is written with stray carriage returns; splitting on \n alone
    # doubles the row count, so let pandas handle the line terminator and strip.
    registry = pd.read_csv(args.registry, sep="\t", dtype=str)
    registry["gse"] = registry["gse"].fillna("").str.strip()
    registry["isoquant_status"] = registry["isoquant_status"].fillna("").str.strip()

    in_scope = set(authority["gse"])
    runs = registry[
        registry["gse"].isin(in_scope)
        & registry["isoquant_status"].str.lower().eq("done")
    ].drop_duplicates("srr").copy()

    missing = sorted(in_scope - set(runs["gse"]))
    if missing:
        raise ValueError(f"authoritative studies with no completed run: {missing}")

    species = set(runs["species"].dropna())
    if not species <= {"human", "mouse"}:
        raise ValueError(f"release contains species outside human/mouse: {species}")

    runs["release"] = RELEASE_ID
    for column in MANIFEST_COLUMNS:
        if column not in runs.columns:
            runs[column] = ""
    manifest = runs[MANIFEST_COLUMNS].sort_values(["gse", "srr"]).reset_index(drop=True)

    # Study-level cell authority, with the per-run pipeline sum kept alongside so
    # the gap between the two is visible instead of being asserted away.
    per_run = manifest.assign(
        cells=pd.to_numeric(manifest["isoquant_cells"], errors="coerce")
    )
    study = (
        per_run.groupby("gse")
        .agg(n_runs=("srr", "nunique"),
             species=("species", lambda s: "|".join(sorted(set(s.dropna())))),
             platforms=("platform", lambda s: "|".join(sorted(set(s.dropna())))),
             pipeline_cell_sum=("cells", "sum"))
        .reset_index()
        .merge(authority.rename(columns={"cell_count": "cells",
                                         "method": "cell_count_method"}),
               on="gse", how="left")
    )
    study["pipeline_cell_sum"] = study["pipeline_cell_sum"].astype("Int64")
    study["cells_agree_with_pipeline"] = (
        study["cells"] == study["pipeline_cell_sum"]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "release_manifest.tsv"
    study_path = args.output_dir / "release_study_cells.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    study.to_csv(study_path, sep="\t", index=False)

    by_species = study.groupby("species").agg(
        studies=("gse", "nunique"), runs=("n_runs", "sum"), cells=("cells", "sum"))
    summary = {
        "release": RELEASE_ID,
        "scope": "all single-cell long-read datasets, any source archive",
        "authoritative_table": str(args.authoritative),
        "authoritative_sha256": sha256(args.authoritative),
        "n_studies": int(study["gse"].nunique()),
        "n_runs": int(manifest["srr"].nunique()),
        "n_cells": int(study["cells"].sum()),
        "by_species": {
            index: {k: int(v) for k, v in row.items()}
            for index, row in by_species.iterrows()
        },
        "by_platform": {
            str(k): int(v)
            for k, v in manifest["platform"].value_counts().items()
        },
        "cell_count_method": {
            str(k): int(v)
            for k, v in study["cell_count_method"].value_counts().items()
        },
        "studies_where_pipeline_sum_differs": int(
            (~study["cells_agree_with_pipeline"]).sum()
        ),
        "manifest_sha256": sha256(manifest_path),
        "study_table_sha256": sha256(study_path),
    }
    (args.output_dir / "release_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
