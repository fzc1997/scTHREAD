#!/usr/bin/env python3
"""Build the immutable-run candidate manifest used by the NAR draft."""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


# Membership frozen on 27 July 2026 from the .bak_2338 registry.
EXPECTED = {
    "runs": 469,
    "studies": 31,
    "cells": 845_781,
    "historical_unresolved_species": 5,
}

# A modality audit (2026-08-01) showed that GSE140890 contributes 35 genomic-DNA
# MinION libraries, not transcriptomes. NCBI SRA reports Assay Type OTHER,
# LibrarySource GENOMIC and LibrarySelection other for every one of them
# (SRP233378 / PRJNA591948; the study's 77 RNA runs are short-read Illumina and
# were never in the snapshot). They are excluded from the released manifest.
RELEASE_ID = "scTHREAD_NAR_candidate_20260801"
NCBI_EVIDENCE = Path("tables/modality_audit/GSE140890_SRP233378_SraRunTable_ncbi.csv")
EXPECTED_RNA = {
    "runs": 434,
    "studies": 30,
    "cells": 850_938,
    "excluded_runs": 35,
    "excluded_studies": ["GSE140890"],
}

# A second audit (2026-08-01) found one frozen row pointing at an IsoQuant
# directory that does not exist, so its cell and transcript counts were stale.
# Rather than rewrite the frozen registry, rows whose frozen output directory is
# missing while the current registry's directory exists are re-read from the
# current registry and verified against the files on disk.
ISOQUANT_ROOT = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong"))
EXPECTED_STALE_PATH_RUNS = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-registry", type=Path, required=True)
    parser.add_argument("--current-registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    frozen = pd.read_csv(args.frozen_registry, sep="\t", low_memory=False)
    current = pd.read_csv(args.current_registry, sep="\t", low_memory=False)
    selected = frozen[
        frozen["isoquant_status"].astype(str).str.lower().eq("done")
    ].copy()

    if len(selected) != EXPECTED["runs"]:
        raise ValueError(f"expected {EXPECTED['runs']} runs, observed {len(selected)}")
    if selected["srr"].duplicated().any():
        raise ValueError("duplicate run identifiers in frozen candidate")
    if selected["gse"].nunique() != EXPECTED["studies"]:
        raise ValueError("study count does not match frozen contract")
    cell_total = int(pd.to_numeric(selected["isoquant_cells"], errors="raise").sum())
    if cell_total != EXPECTED["cells"]:
        raise ValueError("cell total does not match frozen contract")

    unresolved = selected["species"].astype(str).eq("?")
    if int(unresolved.sum()) != EXPECTED["historical_unresolved_species"]:
        raise ValueError("unexpected number of unresolved historical species rows")

    cell_count_drift: list[dict] = []
    current_by_run = current.set_index("srr", verify_integrity=True)
    for idx, row in selected.loc[unresolved].iterrows():
        run = row["srr"]
        if run not in current_by_run.index:
            raise ValueError(f"unresolved run missing from current registry: {run}")
        now = current_by_run.loc[run]
        if str(now["isoquant_status"]).lower() != "done":
            raise ValueError(f"resolved run no longer complete: {run}")
        if float(now["isoquant_cells"]) != float(row["isoquant_cells"]):
            # The frozen value stays authoritative; only species/description are
            # taken from the current registry, so a later re-count is recorded
            # rather than adopted.
            cell_count_drift.append({
                "srr": run,
                "frozen_cells": float(row["isoquant_cells"]),
                "current_cells": float(now["isoquant_cells"]),
            })
        if now["species"] not in {"human", "mouse"}:
            raise ValueError(f"current species is still unresolved for {run}")
        selected.at[idx, "species"] = now["species"]
        selected.at[idx, "description"] = now["description"]

    if set(selected["species"]) != {"human", "mouse"}:
        raise ValueError("candidate release contains species outside human/mouse")

    # --- stale frozen rows: adopt the current registry where the files prove it ---
    stale_rows = []
    for idx, row in selected.iterrows():
        run = str(row["srr"])
        frozen_dir = row.get("isoquant_path")
        if not isinstance(frozen_dir, str) or not frozen_dir:
            continue
        if (ISOQUANT_ROOT / frozen_dir).is_dir():
            continue
        if run not in current_by_run.index:
            continue
        current_dir = current_by_run.at[run, "isoquant_path"]
        if not isinstance(current_dir, str) or not (ISOQUANT_ROOT / current_dir).is_dir():
            raise ValueError(f"frozen IsoQuant directory missing and unrecoverable: {run}")

        barcodes = (
            ISOQUANT_ROOT / current_dir / run
            / f"{run}.transcript_grouped_tag_CB_counts.barcodes.tsv"
        )
        if not barcodes.is_file():
            barcodes = (
                ISOQUANT_ROOT / current_dir
                / f"{run}.transcript_grouped_tag_CB_counts.barcodes.tsv"
            )
        if not barcodes.is_file():
            raise ValueError(f"cannot verify corrected cell count for {run}")
        with barcodes.open() as handle:
            observed_cells = sum(1 for line in handle if line.strip())
        current_cells = int(float(current_by_run.at[run, "isoquant_cells"]))
        if observed_cells != current_cells:
            raise ValueError(
                f"current registry disagrees with the barcode file for {run}: "
                f"{current_cells} vs {observed_cells}"
            )

        stale_rows.append({
            "srr": run,
            "frozen_isoquant_path": frozen_dir,
            "current_isoquant_path": current_dir,
            "frozen_cells": int(float(row["isoquant_cells"])),
            "corrected_cells": observed_cells,
            "frozen_transcripts": row["isoquant_transcripts"],
            "corrected_transcripts": current_by_run.at[run, "isoquant_transcripts"],
            "evidence": str(barcodes),
        })
        selected.at[idx, "isoquant_cells"] = current_cells
        selected.at[idx, "isoquant_transcripts"] = current_by_run.at[run, "isoquant_transcripts"]
        selected.at[idx, "isoquant_path"] = current_dir

    if len(stale_rows) != EXPECTED_STALE_PATH_RUNS:
        raise ValueError(
            f"expected {EXPECTED_STALE_PATH_RUNS} stale-path rows, found {len(stale_rows)}"
        )

    # --- modality audit: keep transcriptomic libraries only -------------------
    modality = current.set_index("srr")
    source = (
        selected["srr"].map(modality["ena_library_source"]).fillna("").str.upper()
    )
    strategy = (
        selected["srr"].map(modality["ena_library_strategy"]).fillna("").str.upper()
    )
    nontranscriptomic = source.eq("GENOMIC") | strategy.eq("OTHER")
    excluded = selected.loc[nontranscriptomic.values].copy()

    if int(len(excluded)) != EXPECTED_RNA["excluded_runs"]:
        raise ValueError(
            f"expected {EXPECTED_RNA['excluded_runs']} non-transcriptomic runs, "
            f"observed {len(excluded)}"
        )
    if sorted(excluded["gse"].unique()) != EXPECTED_RNA["excluded_studies"]:
        raise ValueError(f"unexpected excluded studies: {sorted(excluded['gse'].unique())}")

    # every exclusion must be confirmed by the archived NCBI SRA Run Table
    if not NCBI_EVIDENCE.exists():
        raise ValueError(f"NCBI modality evidence missing: {NCBI_EVIDENCE}")
    ncbi = pd.read_csv(NCBI_EVIDENCE).set_index("Run")
    unconfirmed = [
        run
        for run in excluded["srr"]
        if run not in ncbi.index or ncbi.at[run, "LibrarySource"] != "GENOMIC"
    ]
    if unconfirmed:
        raise ValueError(f"exclusions not confirmed by NCBI SRA: {unconfirmed}")

    selected = selected.loc[~nontranscriptomic.values].copy()
    if len(selected) != EXPECTED_RNA["runs"]:
        raise ValueError(f"expected {EXPECTED_RNA['runs']} RNA runs, observed {len(selected)}")
    if selected["gse"].nunique() != EXPECTED_RNA["studies"]:
        raise ValueError("RNA study count does not match the audited contract")
    rna_cells = int(pd.to_numeric(selected["isoquant_cells"], errors="raise").sum())
    if rna_cells != EXPECTED_RNA["cells"]:
        raise ValueError(f"RNA cell total {rna_cells} does not match the audited contract")

    keep = [
        "srr",
        "gse",
        "species",
        "platform",
        "biology_group",
        "description",
        "pipeline",
        "isoquant_status",
        "isoquant_cells",
        "isoquant_cells_method",
        "isoquant_transcripts",
        "annotation_status",
        "annotation_method",
        "annotation_cells",
    ]
    manifest = selected[keep].sort_values(["gse", "srr"]).reset_index(drop=True)
    manifest.insert(0, "release_candidate", RELEASE_ID)
    manifest["metadata_resolution"] = "frozen"
    manifest.loc[
        manifest["srr"].isin(selected.loc[unresolved, "srr"]),
        "metadata_resolution",
    ] = "species_and_description_from_current_registry_exact_run_match"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "release_candidate_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    summary = {
        "release_candidate": RELEASE_ID,
        "status": "candidate_not_public_release",
        "source_registry": {
            "path": str(args.frozen_registry),
            "sha256": sha256(args.frozen_registry),
        },
        "metadata_resolution_registry": {
            "path": str(args.current_registry),
            "sha256": sha256(args.current_registry),
            "resolved_run_count": int(unresolved.sum()),
            "cell_count_drift_recorded_not_adopted": cell_count_drift,
        },
        "counts": {
            "runs": int(len(manifest)),
            "studies": int(manifest["gse"].nunique()),
            "cells": int(manifest["isoquant_cells"].sum()),
            "human_runs": int((manifest["species"] == "human").sum()),
            "mouse_runs": int((manifest["species"] == "mouse").sum()),
            "ont_runs": int((manifest["platform"] == "ONT").sum()),
            "pacbio_runs": int((manifest["platform"] == "PacBio").sum()),
        },
        "manifest": {
            "path": str(manifest_path),
            "rows": int(len(manifest)),
            "sha256": sha256(manifest_path),
        },
        "stale_row_correction": {
            "date": "2026-08-01",
            "rule": (
                "frozen IsoQuant directory missing while the current registry's "
                "directory exists -> adopt the current cells/transcripts/path"
            ),
            "verification": "corrected cell count re-counted from the barcode file",
            "rows": stale_rows,
        },
        "modality_audit": {
            "date": "2026-08-01",
            "rule": "drop runs whose NCBI/ENA library source is GENOMIC or strategy is OTHER",
            "evidence": str(NCBI_EVIDENCE),
            "evidence_sha256": sha256(NCBI_EVIDENCE),
            "excluded_runs": int(len(excluded)),
            "excluded_studies": sorted(excluded["gse"].unique()),
            "membership_before_audit": {
                "runs": EXPECTED["runs"], "studies": EXPECTED["studies"], "cells": EXPECTED["cells"],
            },
        },
        "scope_note": (
            "Run membership was frozen on 27 July 2026 and restricted on 1 August 2026 "
            "to transcriptomic libraries. It is not the final public release, license "
            "or download bundle."
        ),
    }
    summary_path = args.output_dir / "release_candidate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
