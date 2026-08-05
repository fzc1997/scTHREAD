#!/usr/bin/env python3
"""Recount positive IsoQuant transcript and gene detections in the frozen release."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "tables/release_content/run_layer_contract.tsv"
DEFAULT_OUTPUT_DIR = ROOT / "tables/release_content"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_transcript(feature_id: str) -> str:
    if feature_id.startswith(("ENST", "ENSMUST")):
        return "reference_ensembl"
    if feature_id.startswith("transcript"):
        return "isoquant_novel"
    return "other"


def scan_counts(path: Path, feature_kind: str) -> dict[str, object]:
    feature_counts: defaultdict[str, float] = defaultdict(float)
    n_rows = 0
    duplicate_headers_skipped = 0

    # A zero-byte table is IsoQuant reporting nothing of this kind for the run,
    # not a malformed file: it contributes no features and no molecules.
    if path.stat().st_size == 0:
        empty: dict[str, object] = {
            "enumerated_features": 0, "positive_features": 0,
            "positive_molecules": 0.0, "positive_ids": set(),
            "duplicate_headers_skipped": 0,
        }
        if feature_kind == "transcript":
            for label in ("reference_ensembl", "isoquant_novel", "other"):
                empty[f"{label}_features"] = 0
                empty[f"{label}_molecules"] = 0.0
        return empty

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["feature_id", "count"]:
            raise RuntimeError(f"Unexpected count-table schema in {path}: {reader.fieldnames}")
        for row in reader:
            feature_id = row["feature_id"]
            if feature_id == "feature_id" and row["count"] == "count":
                duplicate_headers_skipped += 1
                continue
            if feature_id.startswith("__"):
                continue
            n_rows += 1
            try:
                count = float(row["count"])
            except ValueError as error:
                raise RuntimeError(
                    f"Non-numeric count in {path}: feature_id={feature_id!r}, "
                    f"count={row['count']!r}"
                ) from error
            feature_counts[feature_id] += count

    positive_ids = {
        feature_id for feature_id, count in feature_counts.items() if count > 0
    }
    positive_molecules = sum(
        count for count in feature_counts.values() if count > 0
    )

    result: dict[str, object] = {
        "enumerated_features": n_rows,
        "positive_features": len(positive_ids),
        "positive_molecules": positive_molecules,
        "positive_ids": positive_ids,
        "duplicate_headers_skipped": duplicate_headers_skipped,
    }
    if feature_kind == "transcript":
        for label in ("reference_ensembl", "isoquant_novel", "other"):
            ids = {
                feature_id
                for feature_id in positive_ids
                if classify_transcript(feature_id) == label
            }
            result[f"{label}_features"] = len(ids)
            result[f"{label}_molecules"] = sum(
                feature_counts[feature_id] for feature_id in ids
            )
    return result


def build(contract_path: Path, output_dir: Path, selected_srr: list[str]) -> None:
    contract = pd.read_csv(contract_path, sep="\t", low_memory=False)
    eligible = contract.loc[contract["file_contract_status"].eq("isoquant_core_complete")].copy()
    if selected_srr:
        unknown = sorted(set(selected_srr) - set(eligible["srr"]))
        if unknown:
            raise ValueError(f"Requested runs are not feature-complete: {unknown}")
        eligible = eligible.loc[eligible["srr"].isin(selected_srr)].copy()
    elif len(eligible) != 155:
        raise RuntimeError(f"Expected 155 feature-complete runs, found {len(eligible)}")

    rows: list[dict[str, object]] = []
    reference_transcripts: defaultdict[str, set[str]] = defaultdict(set)
    observed_genes: defaultdict[str, set[str]] = defaultdict(set)
    reference_transcripts_by_stratum: defaultdict[
        tuple[str, str], set[str]
    ] = defaultdict(set)
    observed_genes_by_stratum: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for _, row in eligible.sort_values(["gse", "srr"]).iterrows():
        transcript = scan_counts(Path(row["transcript_counts_path"]), "transcript")
        discovered = scan_counts(
            Path(row["discovered_transcript_counts_path"]), "transcript"
        )
        gene = scan_counts(Path(row["gene_counts_path"]), "gene")
        expected = int(row["isoquant_transcripts_registry"])
        recounted = int(transcript["positive_features"])
        species = str(row["species"])
        platform = str(row["platform"])
        reference_ids = {
            feature_id
            for feature_id in discovered["positive_ids"]
            if classify_transcript(feature_id) == "reference_ensembl"
        }
        reference_transcripts[species].update(reference_ids)
        reference_transcripts_by_stratum[(species, platform)].update(reference_ids)
        observed_genes[species].update(gene["positive_ids"])
        observed_genes_by_stratum[(species, platform)].update(gene["positive_ids"])
        rows.append(
            {
                "srr": row["srr"],
                "gse": row["gse"],
                "species": species,
                "platform": row["platform"],
                "isoquant_cells": int(row["isoquant_cells"]),
                "registry_positive_reference_transcripts": expected,
                "transcript_rows_enumerated": transcript["enumerated_features"],
                "observed_run_transcript_detections": recounted,
                "registry_recount_delta": recounted - expected,
                "registry_recount_match": recounted == expected,
                "observed_reference_transcript_detections": transcript[
                    "reference_ensembl_features"
                ],
                "observed_isoquant_novel_transcript_detections": transcript[
                    "isoquant_novel_features"
                ],
                "observed_other_transcript_detections": transcript["other_features"],
                "transcript_assigned_molecules": transcript["positive_molecules"],
                "transcript_duplicate_headers_skipped": transcript[
                    "duplicate_headers_skipped"
                ],
                "discovered_transcript_rows_enumerated": discovered[
                    "enumerated_features"
                ],
                "observed_discovered_run_transcript_detections": discovered[
                    "positive_features"
                ],
                "observed_discovered_reference_transcript_detections": discovered[
                    "reference_ensembl_features"
                ],
                "observed_discovered_isoquant_novel_transcript_detections": discovered[
                    "isoquant_novel_features"
                ],
                "observed_discovered_other_transcript_detections": discovered[
                    "other_features"
                ],
                "discovered_transcript_assigned_molecules": discovered[
                    "positive_molecules"
                ],
                "discovered_duplicate_headers_skipped": discovered[
                    "duplicate_headers_skipped"
                ],
                "gene_rows_enumerated": gene["enumerated_features"],
                "observed_run_gene_detections": gene["positive_features"],
                "gene_assigned_molecules": gene["positive_molecules"],
                "gene_duplicate_headers_skipped": gene[
                    "duplicate_headers_skipped"
                ],
                "transcript_counts_path": row["transcript_counts_path"],
                "discovered_transcript_counts_path": row[
                    "discovered_transcript_counts_path"
                ],
                "gene_counts_path": row["gene_counts_path"],
            }
        )

    inventory = pd.DataFrame(rows)
    summary = (
        inventory.groupby(["species", "platform"], as_index=False)
        .agg(
            n_runs=("srr", "nunique"),
            n_studies=("gse", "nunique"),
            n_cells=("isoquant_cells", "sum"),
            run_transcript_detections=("observed_run_transcript_detections", "sum"),
            reference_transcript_detections=(
                "observed_reference_transcript_detections",
                "sum",
            ),
            isoquant_novel_transcript_detections=(
                "observed_isoquant_novel_transcript_detections",
                "sum",
            ),
            other_transcript_detections=("observed_other_transcript_detections", "sum"),
            transcript_assigned_molecules=("transcript_assigned_molecules", "sum"),
            discovered_run_transcript_detections=(
                "observed_discovered_run_transcript_detections",
                "sum",
            ),
            discovered_reference_transcript_detections=(
                "observed_discovered_reference_transcript_detections",
                "sum",
            ),
            discovered_isoquant_novel_transcript_detections=(
                "observed_discovered_isoquant_novel_transcript_detections",
                "sum",
            ),
            discovered_other_transcript_detections=(
                "observed_discovered_other_transcript_detections",
                "sum",
            ),
            discovered_transcript_assigned_molecules=(
                "discovered_transcript_assigned_molecules",
                "sum",
            ),
            run_gene_detections=("observed_run_gene_detections", "sum"),
            gene_assigned_molecules=("gene_assigned_molecules", "sum"),
        )
    )
    summary["unique_reference_transcripts_observed"] = summary.apply(
        lambda row: len(
            reference_transcripts_by_stratum[(row["species"], row["platform"])]
        ),
        axis=1,
    )
    summary["unique_gene_ids_observed"] = summary.apply(
        lambda row: len(observed_genes_by_stratum[(row["species"], row["platform"])]),
        axis=1,
    )

    release_rows = [
        (
            "feature_complete_runs",
            len(inventory),
            "runs",
            "Runs with all six contracted IsoQuant output types.",
        ),
        (
            "feature_complete_studies",
            inventory["gse"].nunique(),
            "studies",
            "Studies represented by feature-complete runs.",
        ),
        (
            "feature_complete_cells",
            int(inventory["isoquant_cells"].sum()),
            "called cells",
            "Sum of frozen isoquant_cells across feature-complete runs.",
        ),
        (
            "reference_run_transcript_detections",
            int(inventory["observed_run_transcript_detections"].sum()),
            "positive run-by-reference-transcript rows",
            "Positive rows in transcript_counts.tsv; repeated Ensembl identities in different runs remain separate detections.",
        ),
        (
            "discovered_run_transcript_detections",
            int(inventory["observed_discovered_run_transcript_detections"].sum()),
            "positive run-by-discovered-transcript rows",
            "Positive rows in discovered_transcript_counts.tsv; this is not a globally deduplicated isoform count.",
        ),
        (
            "discovered_isoquant_novel_run_transcript_detections",
            int(
                inventory[
                    "observed_discovered_isoquant_novel_transcript_detections"
                ].sum()
            ),
            "positive run-by-IsoQuant-novel-transcript rows",
            "Novel IDs are run-local and are not deduplicated across runs.",
        ),
        (
            "unique_reference_transcripts_observed",
            sum(len(values) for values in reference_transcripts.values()),
            "species-qualified Ensembl transcript IDs",
            "Union within species of positive Ensembl transcript IDs.",
        ),
        (
            "unique_gene_ids_observed",
            sum(len(values) for values in observed_genes.values()),
            "species-qualified gene IDs",
            "Union within species of positive gene-count feature IDs.",
        ),
    ]
    headline = pd.DataFrame(
        release_rows, columns=["metric", "value", "unit", "interpretation"]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = ".smoke" if selected_srr else ""
    inventory_path = output_dir / f"run_molecular_inventory{suffix}.tsv"
    summary_path = output_dir / f"molecular_inventory_by_stratum{suffix}.tsv"
    headline_path = output_dir / f"molecular_inventory_headline{suffix}.tsv"
    manifest_path = output_dir / f"molecular_inventory_manifest{suffix}.json"
    inventory.to_csv(inventory_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    headline.to_csv(headline_path, sep="\t", index=False)

    payload = {
        "status": "smoke" if selected_srr else "formal",
        "ase_policy": "not accessed",
        "source_contract": {
            "path": str(contract_path.resolve()),
            "sha256": sha256(contract_path),
        },
        "outputs": {
            path.name: {"rows": sum(1 for _ in path.open()) - 1, "sha256": sha256(path)}
            for path in (inventory_path, summary_path, headline_path)
        },
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--srr", action="append", default=[])
    args = parser.parse_args()
    build(args.contract.resolve(), args.output_dir.resolve(), args.srr)


if __name__ == "__main__":
    main()
