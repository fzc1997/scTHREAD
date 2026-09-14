#!/usr/bin/env python3
"""Resolve run metadata and assembly-matched resources for evidence reducers."""

from __future__ import annotations

import csv
import os
from functools import lru_cache
from pathlib import Path


SCLONG = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "."))
REGISTRY = SCLONG / "docs/sample_registry.tsv"
ANNOTATION = (
    SCLONG
    / "metadata/annotation_recovery/release_v1/active_cell_annotation.parquet"
)
OVERRIDES = Path(__file__).resolve().parents[1] / "ann_override.tsv"
REFERENCE_ROOT = Path(os.environ.get("SCTHREAD_REFERENCE_ROOT", "."))

RESOURCES = {
    "human": {
        "assembly": "GRCh38",
        "gene_prefix": "ENSG",
        "transcript_prefix": "ENST",
        "gtf": REFERENCE_ROOT / "10x_ref/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz",
        "polyasite": REFERENCE_ROOT / "endpoint_atlas/hg38/atlas.clusters.2.0.GRCh38.96.bed.gz",
    },
    "mouse": {
        "assembly": "GRCm39",
        "gene_prefix": "ENSMUSG",
        "transcript_prefix": "ENSMUST",
        "gtf": REFERENCE_ROOT / "10x_ref/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz",
        "polyasite": REFERENCE_ROOT / "endpoint_atlas/GRCm39/PolyASite2.0_mm39.bed",
    },
}


@lru_cache(maxsize=1)
def registry_by_run() -> dict[str, dict[str, str]]:
    with REGISTRY.open(newline="") as handle:
        return {
            row["srr"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("srr")
        }


def run_config(run: str) -> tuple[dict[str, str], dict[str, object]]:
    row = registry_by_run().get(run)
    if row is None:
        raise SystemExit(f"RUN_NOT_IN_REGISTRY {run}")
    species = row.get("species", "").strip().lower()
    if species not in RESOURCES:
        raise SystemExit(f"UNSUPPORTED_SPECIES {run} {species or 'missing'}")
    resources = RESOURCES[species]
    for key in ("gtf", "polyasite"):
        if not Path(resources[key]).is_file():
            raise SystemExit(f"MISSING_{key.upper()} {resources[key]}")
    return row, resources


def isoquant_matrix_files(run: str, row: dict[str, str]) -> dict[str, Path]:
    raw_path = row.get("isoquant_path", "").strip()
    if not raw_path:
        raise SystemExit(f"MISSING_ISOQUANT_PATH {run}")
    root = Path(raw_path)
    if not root.is_absolute():
        root = SCLONG / root
    prefixes = (root / run, root / run / run)
    suffixes = {
        "matrix": "transcript_grouped_tag_CB_counts.matrix.mtx",
        "features": "transcript_grouped_tag_CB_counts.features.tsv",
        "barcodes": "transcript_grouped_tag_CB_counts.barcodes.tsv",
    }
    for prefix in prefixes:
        files = {
            key: Path(f"{prefix}.{suffix}")
            for key, suffix in suffixes.items()
        }
        if all(path.is_file() for path in files.values()):
            return files
    raise SystemExit(f"MISSING_TRANSCRIPT_MATRIX {run} {root}")


@lru_cache(maxsize=1)
def annotation_overrides() -> dict[str, str]:
    if not OVERRIDES.is_file():
        return {}
    with OVERRIDES.open(newline="") as handle:
        return {
            row["run"]: row["pattern"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def annotation_pattern(run: str, row: dict[str, str]) -> str:
    override = annotation_overrides().get(run)
    if override is not None:
        return override
    if row.get("annotation_scope", "").strip().lower() == "run_exact":
        return f"%::{run}::%"
    return "%"


def resolved_read_facts(run: str) -> Path:
    path = (
        SCLONG
        / "results/science_first/phase1_resolved_reads_v4"
        / run
        / "resolved_read_facts.parquet"
    )
    if not path.is_file():
        raise SystemExit(f"MISSING_RESOLVED_READ_FACTS {run}")
    return path
