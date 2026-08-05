#!/usr/bin/env python3
"""Test whether the catalog recovers the classical CD45 RA/RO isoform split.

CD45 isoforms are defined by inclusion or exclusion of the three variable exons
4, 5 and 6 (A, B and C), which is the one fact about PTPRC that is fixed by
prior work rather than by this database. Two annotated models in the catalog
differ by exactly those three exons, so their relative usage across cell types
is a check that can only come out one way if the records are sound: lymphoid
cells should favour the exon-containing model and monocytes the shorter one.

The exon identities are derived from the annotation here rather than asserted,
so the check fails loudly if the transcript pair or the annotation build change.

Writes tables/cd45_ra_ro_recovery.{tsv,json}.
"""
from __future__ import annotations

import os

import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GTF = Path(os.environ.get("SCTHREAD_GENCODE_ROOT", "/gpfs/home/fuzc/project/ABE_transcriptomics_off_target/dep_files") + "/"
           "hg38/hg38_annotation/gencode.v41.annotation.gtf")
USAGE = ROOT / "tables/Table_S_PTPRC_isoform_usage.tsv"

# The long model carries the variable exons; the short model lacks all three.
LONG, SHORT = "ENST00000442510", "ENST00000348564"
# Canonical lengths of PTPRC exons 4, 5 and 6 (A, B, C).
VARIABLE_EXON_BP = (198, 141, 144)


def ptprc_exons() -> tuple[dict[str, list[tuple[int, int]]], dict[str, str]]:
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    names: dict[str, str] = {}
    with GTF.open() as handle:
        for line in handle:
            if 'gene_name "PTPRC"' not in line:
                continue
            fields = line.split("\t")
            match = re.search(r'transcript_id "([^.]+)', fields[8])
            if not match or match.group(1) not in (LONG, SHORT):
                continue
            tid = match.group(1)
            if fields[2] == "transcript":
                names[tid] = re.search(r'transcript_name "([^"]+)', fields[8]).group(1)
            elif fields[2] == "exon":
                exons[tid].append((int(fields[3]), int(fields[4])))
    return exons, names


def main() -> None:
    exons, names = ptprc_exons()
    for tid in (LONG, SHORT):
        if not exons.get(tid):
            raise SystemExit(f"{tid} not found in {GTF.name}")

    only_long = sorted(set(exons[LONG]) - set(exons[SHORT]))
    variable = [(s, e) for s, e in only_long if (e - s + 1) in VARIABLE_EXON_BP]
    if len(variable) != 3:
        raise SystemExit(
            f"expected the three variable exons {VARIABLE_EXON_BP} between "
            f"{names[LONG]} and {names[SHORT]}, found "
            f"{[e - s + 1 for s, e in only_long]}"
        )

    usage = pd.read_csv(USAGE, sep="\t")
    pair = usage[usage.transcript_id.isin([LONG, SHORT])]
    wide = pair.pivot(index="ct", columns="transcript_id", values="count").fillna(0)
    wide["molecules"] = wide[LONG] + wide[SHORT]
    wide["short_model_fraction"] = (wide[SHORT] / wide["molecules"]).round(4)
    wide = wide.rename(columns={LONG: f"{names[LONG]}_molecules",
                                SHORT: f"{names[SHORT]}_molecules"})
    wide = wide.sort_values("short_model_fraction")
    wide.to_csv(ROOT / "tables/cd45_ra_ro_recovery.tsv", sep="\t")

    result = {
        "long_model": {"transcript": LONG, "name": names[LONG],
                       "exons": len(exons[LONG])},
        "short_model": {"transcript": SHORT, "name": names[SHORT],
                        "exons": len(exons[SHORT])},
        "variable_exons_bp": [e - s + 1 for s, e in variable],
        "variable_exon_coords": [f"chr1:{s}-{e}" for s, e in variable],
        "molecules_total": int(wide["molecules"].sum()),
        "short_model_fraction": {ct: float(v) for ct, v
                                 in wide["short_model_fraction"].items()},
    }
    (ROOT / "tables/cd45_ra_ro_recovery.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print()
    print(wide.to_string())


if __name__ == "__main__":
    main()
