#!/usr/bin/env python3
"""Recount the PTPRC variable-exon choice with and without CB/UB deduplication.

The read-level classifier is shared with ``measure_cd45_exon_inclusion.py``.
UMI counts use one key per (corrected cell barcode, UB, splice-choice class);
keys whose reads support more than one class are reported as conflicts and are
excluded from the deduplicated class counts rather than resolved by dominance.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import pysam

from measure_cd45_exon_inclusion import SCAN, classify


ALIASES = {
    "Monocyte / myeloid": "Monocyte",
    "Monocyte/Myeloid": "Monocyte",
    "NK cell": "NK",
    "CD4 T": "T cell",
    "CD8 T": "T cell",
    "Malignant plasma cell": "Plasma cell",
}


def load_annotations(path: Path, run: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("run") != run:
                continue
            cell_type = ALIASES.get(row.get("cell_type", ""), row.get("cell_type", ""))
            if not cell_type or cell_type in {"Unassigned", "Unlabeled"}:
                continue
            barcode = (row.get("barcode") or "").strip()
            if not barcode:
                continue
            mapping[barcode] = cell_type
            mapping[barcode.removesuffix("-1")] = cell_type
    return mapping


def count_bam(bam: Path, annotations: dict[str, str]) -> tuple[dict[tuple[str, str], int], dict[tuple[str, str], set[str]], dict[str, int]]:
    reads: dict[tuple[str, str], int] = defaultdict(int)
    umi_labels: dict[tuple[str, str], set[str]] = defaultdict(set)
    missing_ub: dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(bam) as handle:
        for read in handle.fetch(*SCAN):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            try:
                barcode = str(read.get_tag("CB"))
            except KeyError:
                continue
            cell_type = annotations.get(barcode) or annotations.get(barcode.removesuffix("-1"))
            if cell_type is None:
                continue
            label = classify(read)
            if label is None:
                continue
            reads[(cell_type, label)] += 1
            try:
                umi = str(read.get_tag("UB"))
            except KeyError:
                missing_ub[cell_type] += 1
                continue
            umi_labels[(cell_type, f"{barcode}\t{umi}")].add(label)
    return dict(reads), dict(umi_labels), dict(missing_ub)


def process_run(bam: Path, annotation: Path, study: str) -> tuple[list[dict[str, object]], dict[str, object]]:
    run = bam.name.removesuffix(".tagged.bam")
    labels = load_annotations(annotation, run)
    if not labels:
        return [], {"study": study, "run": run, "status": "no_annotation_rows"}
    reads, umi_labels, missing_ub = count_bam(bam, labels)
    umi_counts: dict[tuple[str, str], int] = defaultdict(int)
    conflicts: dict[str, int] = defaultdict(int)
    for (cell_type, _key), classes in umi_labels.items():
        if len(classes) == 1:
            umi_counts[(cell_type, next(iter(classes)))] += 1
        else:
            conflicts[cell_type] += 1
    keys = sorted(set(reads) | set(umi_counts))
    rows = [
        {
            "study": study,
            "run": run,
            "cell_type": cell_type,
            "class": label,
            "reads": reads.get((cell_type, label), 0),
            "umi_molecules": umi_counts.get((cell_type, label), 0),
            "conflicted_umis": conflicts.get(cell_type, 0),
        }
        for cell_type, label in keys
    ]
    manifest = {
        "study": study,
        "run": run,
        "status": "ok",
        "annotated_barcodes": len(labels),
        "read_informative": sum(reads.values()),
        "umi_informative": sum(umi_counts.values()),
        "conflicted_umis": sum(conflicts.values()),
        "reads_without_ub": sum(missing_ub.values()),
    }
    return rows, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam-root", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--run", help="Process one run for a bounded smoke test")
    args = parser.parse_args()

    bams = sorted(args.bam_root.glob("*/*/*.tagged.bam"))
    if args.run:
        bams = [path for path in bams if path.name == f"{args.run}.tagged.bam"]
    if not bams:
        raise SystemExit("No tagged BAM matched the requested root/run")

    rows: list[dict[str, object]] = []
    manifests: list[dict[str, object]] = []
    for bam in bams:
        run_rows, run_manifest = process_run(bam, args.annotation, args.study)
        rows.extend(run_rows)
        manifests.append(run_manifest)
        print(json.dumps(run_manifest, sort_keys=True), flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["study", "run", "cell_type", "class", "reads", "umi_molecules", "conflicted_umis"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    args.manifest.write_text(json.dumps({"runs": manifests}, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
