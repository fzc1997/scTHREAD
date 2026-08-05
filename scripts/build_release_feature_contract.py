#!/usr/bin/env python3
"""Build a frozen run-by-layer file-availability contract for the NAR release."""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD"))
NAR_ROOT = PROJECT / "NAR_database"
SEQ_ROOT = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong"))
DEFAULT_MANIFEST = (
    NAR_ROOT
    / "tables/release_20260803/release_manifest.tsv"
)
DEFAULT_REGISTRY = SEQ_ROOT / "docs/sample_registry.tsv"
DEFAULT_OUTPUT_DIR = NAR_ROOT / "tables/release_content"
ANALYSIS_ROOT = PROJECT / "results/paper1/f2_grammar"

CORE_PATTERNS = {
    "gene_counts": "{srr}.gene_counts.tsv",
    "transcript_counts": "{srr}.transcript_counts.tsv",
    "read_assignments": "{srr}.read_assignments.tsv.gz",
    "transcript_models": "{srr}.transcript_models.gtf",
    "extended_annotation": "{srr}.extended_annotation.gtf",
    "discovered_transcript_counts": "{srr}.discovered_transcript_counts.tsv",
}
ANALYSIS_PATTERNS = {
    "diu": (ANALYSIS_ROOT / "agg_diu", "{srr}.diu.parquet"),
    "apa": (ANALYSIS_ROOT / "agg_apa", "{srr}.apa.parquet"),
    "junction": (ANALYSIS_ROOT / "agg_jct", "{srr}.jct_ct.parquet"),
    # Availability metadata only. No ASE content or detectability analysis is run.
    "ase": (ANALYSIS_ROOT / "agg_ase_ct", "{srr}.ase_ct.parquet"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def existing_file(candidates: list[Path]) -> Path | None:
    for path in candidates:
        if path.is_file():
            return path.resolve()
    return None


def existing_pattern_file(
    search_dirs: list[Path], pattern: str, srr: str
) -> tuple[Path | None, str]:
    for directory in search_dirs:
        stems = list(dict.fromkeys([srr, directory.name]))
        for stem in stems:
            path = directory / pattern.format(srr=stem)
            if path.is_file():
                return path.resolve(), stem
    return None, ""


def file_fields(prefix: str, path: Path | None) -> dict[str, object]:
    return {
        f"{prefix}_available": path is not None,
        f"{prefix}_path": str(path) if path else "",
        f"{prefix}_bytes": path.stat().st_size if path else 0,
    }


def resolve_isoquant_dir(raw_path: object, srr: str) -> tuple[Path | None, list[Path]]:
    if pd.isna(raw_path) or not str(raw_path).strip():
        return None, []
    base = Path(str(raw_path))
    if not base.is_absolute():
        base = SEQ_ROOT / base
    candidates = list(
        dict.fromkeys(
            [
                base.resolve(),
                (base / srr).resolve(),
                (base / base.name).resolve(),
            ]
        )
    )
    existing = [path for path in candidates if path.is_dir()]
    return base.resolve(), existing


def build_row(row: pd.Series, registry_row: pd.Series | None) -> dict[str, object]:
    srr = str(row["srr"])
    documented = pd.notna(row["isoquant_transcripts"])
    raw_path = registry_row["isoquant_path"] if registry_row is not None else ""
    registered_dir, search_dirs = resolve_isoquant_dir(raw_path, srr)

    result: dict[str, object] = {
        "release_candidate": row["release"],
        "srr": srr,
        "gse": row["gse"],
        "species": row["species"],
        "platform": row["platform"],
        "isoquant_cells": int(row["isoquant_cells"]),
        "isoquant_transcripts_documented": documented,
        "isoquant_transcripts_registry": (
            int(row["isoquant_transcripts"]) if documented else ""
        ),
        "isoquant_path_registry": str(raw_path) if not pd.isna(raw_path) else "",
        "isoquant_path_resolved": str(registered_dir) if registered_dir else "",
        "isoquant_dir_available": bool(search_dirs),
        "isoquant_search_dir": str(search_dirs[-1]) if search_dirs else "",
        "isoquant_file_stem": "",
    }

    for layer, pattern in CORE_PATTERNS.items():
        path, stem = existing_pattern_file(search_dirs, pattern, srr)
        if stem and not result["isoquant_file_stem"]:
            result["isoquant_file_stem"] = stem
        result.update(file_fields(layer, path))

    for layer, (directory, pattern) in ANALYSIS_PATTERNS.items():
        path = existing_file([directory / pattern.format(srr=srr)])
        result.update(file_fields(layer, path))

    core = [
        bool(result[f"{layer}_available"])
        for layer in (
            "gene_counts",
            "transcript_counts",
            "read_assignments",
            "transcript_models",
        )
    ]
    if not documented:
        status = "membership_only"
    elif all(core):
        status = "isoquant_core_complete"
    elif any(core):
        status = "isoquant_core_partial"
    else:
        status = "isoquant_core_missing"
    result["file_contract_status"] = status
    return result


def build(
    manifest_path: Path,
    registry_path: Path,
    output_dir: Path,
    selected_srr: list[str],
) -> None:
    manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
    registry = pd.read_csv(registry_path, sep="\t", low_memory=False)
    # Cells are authoritative per study, not per run, so only membership is
    # checked here; the cell total is verified against release_study_cells.tsv.
    if (
        len(manifest) != 453
        or manifest["srr"].nunique() != 453
        or manifest["gse"].nunique() != 34
    ):
        raise RuntimeError("Frozen release invariants failed")
    study_cells = pd.read_csv(
        manifest_path.parent / "release_study_cells.tsv", sep="\t"
    )
    if int(study_cells["cells"].sum()) != 923_389:
        raise RuntimeError("Release cell total changed")
    if registry["srr"].duplicated().any():
        raise RuntimeError("Current registry has duplicate run identifiers")

    if selected_srr:
        unknown = sorted(set(selected_srr) - set(manifest["srr"]))
        if unknown:
            raise ValueError(f"Requested runs are not in the manifest: {unknown}")
        manifest = manifest.loc[manifest["srr"].isin(selected_srr)].copy()

    registry_by_run = registry.set_index("srr", verify_integrity=True)
    rows = []
    for _, row in manifest.sort_values(["gse", "srr"]).iterrows():
        registry_row = (
            registry_by_run.loc[row["srr"]]
            if row["srr"] in registry_by_run.index
            else None
        )
        rows.append(build_row(row, registry_row))
    contract = pd.DataFrame(rows)

    if not selected_srr:
        if len(contract) != 453:
            raise RuntimeError("Full contract does not contain 453 rows")

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".smoke" if selected_srr else ""
    contract_path = output_dir / f"run_layer_contract{suffix}.tsv"
    summary_path = output_dir / f"run_layer_contract_summary{suffix}.tsv"
    manifest_json = output_dir / f"run_layer_contract_manifest{suffix}.json"

    contract.to_csv(contract_path, sep="\t", index=False)
    summary = (
        contract.groupby(
            ["file_contract_status", "species", "platform"],
            dropna=False,
            as_index=False,
        )
        .agg(
            n_runs=("srr", "nunique"),
            n_studies=("gse", "nunique"),
            n_cells=("isoquant_cells", "sum"),
            gene_counts_files=("gene_counts_available", "sum"),
            transcript_counts_files=("transcript_counts_available", "sum"),
            junction_files=("junction_available", "sum"),
            diu_files=("diu_available", "sum"),
            apa_files=("apa_available", "sum"),
            ase_files_metadata_only=("ase_available", "sum"),
        )
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    payload = {
        "status": "smoke" if selected_srr else "formal",
        "ase_policy": "availability metadata only; detectability analysis deferred",
        "source_manifest": {
            "path": str(manifest_path.resolve()),
            "sha256": sha256(manifest_path),
        },
        "source_registry": {
            "path": str(registry_path.resolve()),
            "sha256": sha256(registry_path),
        },
        "outputs": {
            "contract": {
                "path": str(contract_path.resolve()),
                "rows": len(contract),
                "sha256": sha256(contract_path),
            },
            "summary": {
                "path": str(summary_path.resolve()),
                "rows": len(summary),
                "sha256": sha256(summary_path),
            },
        },
        "status_counts": contract["file_contract_status"].value_counts().to_dict(),
    }
    manifest_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--srr", action="append", default=[])
    args = parser.parse_args()
    build(
        args.manifest.resolve(),
        args.registry.resolve(),
        args.output_dir.resolve(),
        args.srr,
    )


if __name__ == "__main__":
    main()
