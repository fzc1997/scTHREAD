#!/usr/bin/env python3
"""Compare the structure of assembled and annotated transcripts in the catalog.

IsoQuant's transcript models carry an `exons` attribute, and reference-derived
models additionally carry `transcript_name`, so each model can be assigned to
one of the two groups without re-deriving anything. Exon count comes from the
attribute; spliced length is summed over the exon lines.

Scope is the blood/marrow studies, the largest biological system in the
catalog, matching the exemplar used in the Results.
"""
from __future__ import annotations

import os

import argparse
import json
import re
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ISOQUANT = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/data")
STUDIES = ["GSE307660", "GSE276974", "GSE283658"]

TID = re.compile(r'transcript_id "([^"]+)"')
EXONS = re.compile(r'exons "(\d+)"')
NAMED = re.compile(r'transcript_name "')


def measure(gtf: Path) -> dict[str, list[tuple[int, int]]]:
    """Return {group: [(exon_count, spliced_length), ...]} for one run."""
    exon_count: dict[str, int] = {}
    group: dict[str, str] = {}
    length: dict[str, int] = {}
    with gtf.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature, start, end, attrs = parts[2], parts[3], parts[4], parts[8]
            match = TID.search(attrs)
            if not match:
                continue
            tid = match.group(1)
            if feature == "transcript":
                exons = EXONS.search(attrs)
                if exons:
                    exon_count[tid] = int(exons.group(1))
                group[tid] = "annotated" if NAMED.search(attrs) else "assembled"
            elif feature == "exon":
                length[tid] = length.get(tid, 0) + (int(end) - int(start) + 1)
    out: dict[str, list[tuple[int, int]]] = {"annotated": [], "assembled": []}
    for tid, grp in group.items():
        if tid in exon_count and tid in length:
            out[grp].append((exon_count[tid], length[tid]))
    return out


def summarize(values: list[int]) -> dict:
    values = sorted(values)
    if not values:
        return {}
    return {
        "n": len(values),
        "median": st.median(values),
        "q1": values[len(values) // 4],
        "q3": values[3 * len(values) // 4],
        "mono_exon_pct": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-runs", type=int, default=8)
    parser.add_argument("--out", type=Path, default=ROOT / "tables/transcript_structure.json")
    args = parser.parse_args()

    gtfs: list[Path] = []
    for study in STUDIES:
        for base in (ISOQUANT / "isoquant" / study, ISOQUANT / "isoquant_pacbio" / study):
            if not base.is_dir():
                continue
            for run_dir in sorted(base.iterdir()):
                hit = run_dir / run_dir.name / f"{run_dir.name}.transcript_models.gtf"
                if hit.is_file():
                    gtfs.append(hit)
    gtfs = gtfs[: args.max_runs]
    if not gtfs:
        raise SystemExit("no transcript_models.gtf found")

    pooled: dict[str, list[tuple[int, int]]] = {"annotated": [], "assembled": []}
    for gtf in gtfs:
        got = measure(gtf)
        for key, rows in got.items():
            pooled[key].extend(rows)
        print(f"{gtf.parent.name}: annotated={len(got['annotated']):,} "
              f"assembled={len(got['assembled']):,}", flush=True)

    result = {"runs": [g.parent.name for g in gtfs], "groups": {}}
    for key, rows in pooled.items():
        exons = [e for e, _ in rows]
        lengths = [l for _, l in rows]
        result["groups"][key] = {
            "n_models": len(rows),
            "exons": summarize(exons),
            "spliced_length_bp": summarize(lengths),
            "multi_exon_pct": round(100 * sum(1 for e in exons if e > 1) / len(exons), 1) if exons else None,
            "exons_ge5_pct": round(100 * sum(1 for e in exons if e >= 5) / len(exons), 1) if exons else None,
        }
    args.out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
