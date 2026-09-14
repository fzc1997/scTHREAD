#!/usr/bin/env python3
"""Summarize the real CD45 read-vs-UMI sensitivity recount."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


SELECTED = ("B cell", "NK", "Monocyte", "T cell")


def read_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    return rows


def summarize(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    totals: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for row in rows:
        if row["cell_type"] not in SELECTED:
            continue
        key = (row["study"], row["run"], row["cell_type"])
        totals[key][0] += int(row["reads"])
        totals[key][1] += int(row["umi_molecules"])

    summary_rows: list[dict[str, object]] = []
    by_study: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for study in sorted({row["study"] for row in rows}):
        for cell_type in SELECTED:
            runs = {
                run
                for (row_study, run, row_cell_type), values in totals.items()
                if row_study == study and row_cell_type == cell_type and values[0] >= 100
            }
            read_total = sum(totals[(study, run, cell_type)][0] for run in runs)
            umi_total = sum(totals[(study, run, cell_type)][1] for run in runs)
            read_ro = sum(
                int(row["reads"])
                for row in rows
                if row["study"] == study and row["cell_type"] == cell_type
                and row["run"] in runs and row["class"] == "RO"
            )
            umi_ro = sum(
                int(row["umi_molecules"])
                for row in rows
                if row["study"] == study and row["cell_type"] == cell_type
                and row["run"] in runs and row["class"] == "RO"
            )
            record = {
                "study": study,
                "cell_type": cell_type,
                "eligible_runs": len(runs),
                "informative_reads": read_total,
                "informative_umi_molecules": umi_total,
                "ro_read_fraction": read_ro / read_total if read_total else None,
                "ro_umi_fraction": umi_ro / umi_total if umi_total else None,
                "absolute_delta_percentage_points": abs(read_ro / read_total - umi_ro / umi_total) * 100 if read_total and umi_total else None,
            }
            summary_rows.append(record)
            by_study[study][cell_type] = record

    pattern_checks = {}
    for study, records in by_study.items():
        read_order = [ct for ct in SELECTED if records[ct]["ro_read_fraction"] is not None]
        umi_order = [ct for ct in SELECTED if records[ct]["ro_umi_fraction"] is not None]
        read_order.sort(key=lambda ct: records[ct]["ro_read_fraction"])
        umi_order.sort(key=lambda ct: records[ct]["ro_umi_fraction"])
        pattern_checks[study] = {
            "read_ro_fraction_order": read_order,
            "umi_ro_fraction_order": umi_order,
            "order_preserved": read_order == umi_order,
        }
    deltas = [r["absolute_delta_percentage_points"] for r in summary_rows if r["absolute_delta_percentage_points"] is not None]
    manifest = {
        "schema": "scTHREAD.cd45_umi_sensitivity.summary.v1",
        "input_scope": "GSE276974 and GSE307660 indexed tagged BAM recounts",
        "selected_cell_types": list(SELECTED),
        "stratum_filter": "at least 100 informative read-level observations per study × run × cell-type stratum",
        "umi_key": "corrected CB + UB + splice-choice class",
        "conflict_rule": "CB/UB keys supporting multiple splice-choice classes are excluded from UMI class counts",
        "pattern_checks": pattern_checks,
        "max_absolute_delta_percentage_points": max(deltas) if deltas else None,
        "all_study_orders_preserved": all(item["order_preserved"] for item in pattern_checks.values()),
        "rows": len(summary_rows),
    }
    return summary_rows, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args()
    rows, manifest = summarize(read_rows(args.input))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
