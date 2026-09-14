#!/usr/bin/env python3
"""Build the Figure 3c-e figure data from the raw CD45 exon-inclusion counts.

This step used to be an ad-hoc command, which is why the inclusion thresholds
appeared only in the figure legend and could not be found anywhere in the
repository. They are constants here, applied once, and written into a sidecar
manifest so the legend, the Methods and the QA gate all read the same rule.

The rule has two levels, and they are deliberately different things:

MIN_STRATUM_READS
    A run x cell-type stratum enters the analysis only with this many
    informative reads. Below it a per-run fraction is not worth reporting.
MIN_CELL_TYPE_READS / MIN_CELL_TYPE_RUNS
    A cell type is shown only if, in *both* studies, its surviving strata carry
    this many reads across this many runs. This is what keeps Dendritic cell,
    Erythroid, Plasma cell and Progenitor out of the figure.

Both levels are applied to the same filtered set, so every read and run count
the figure prints reconciles with every other. An earlier version tested the
cell-type thresholds against unfiltered totals, which made the printed
"informative reads" and "run x cell-type strata" numbers disagree by 876 reads
and by up to 7 runs for one cell type.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]

MIN_STRATUM_READS = 100
MIN_CELL_TYPE_READS = 1_000
MIN_CELL_TYPE_RUNS = 3

# Harmonised because the two studies label the same population differently.
# Left unharmonised, a pivot treats "Monocyte / myeloid" and "Monocyte" as
# different cell types and drops the most concordant population out of the
# shared set, which turns replication into apparent non-replication.
CELL_TYPE_ALIASES = {
    "Monocyte / myeloid": "Monocyte",
    "NK cell": "NK",
    "Malignant plasma cell": "Plasma cell",
}
UNTYPED = "Unassigned"

# The exon-3 donor's first acceptor, derived from which variable exons a read
# covers: a read reaching exon A spliced into A, one covering B but not A
# spliced into B, and so on.
ACCEPTOR_OF = {"RO": "exon 7 (skip all)"}
SKIP = "exon 7 (skip all)"


def first_acceptor(label: str) -> str | None:
    if label in ACCEPTOR_OF:
        return ACCEPTOR_OF[label]
    if label == "other":
        return None
    for exon in "ABC":
        if exon in label:
            return f"exon {exon}"
    return None


def load_counts(inputs: list[Path]) -> pd.DataFrame:
    frame = pd.concat([pd.read_csv(path, sep="\t") for path in inputs])
    frame["cell_type"] = frame.cell_type.replace(CELL_TYPE_ALIASES)
    frame = frame[frame.cell_type.ne(UNTYPED)]
    frame["acceptor"] = frame["class"].map(first_acceptor)
    frame = frame.dropna(subset=["acceptor"]).reset_index(drop=True)
    frame["skip_reads"] = frame.reads.where(frame.acceptor.eq(SKIP), 0)
    return frame


def surviving_strata(frame: pd.DataFrame) -> pd.DataFrame:
    """Reads per run x cell-type stratum, keeping only well-covered strata."""
    strata = (
        frame.groupby(["study", "run", "cell_type"])
        .agg(n=("reads", "sum"), skip=("skip_reads", "sum"))
        .reset_index()
    )
    return strata[strata.n >= MIN_STRATUM_READS].copy()


def select_cell_types(strata: pd.DataFrame, studies: list[str]) -> tuple[list[str], list[dict]]:
    audit = []
    selected = []
    for cell_type in sorted(strata.cell_type.unique()):
        rows = {}
        for study in studies:
            group = strata[strata.study.eq(study) & strata.cell_type.eq(cell_type)]
            rows[study] = {"reads": int(group.n.sum()), "runs": int(group.run.nunique())}
        passes = all(
            rows[s]["reads"] >= MIN_CELL_TYPE_READS and rows[s]["runs"] >= MIN_CELL_TYPE_RUNS
            for s in studies
        )
        audit.append({"cell_type": cell_type, "selected": passes, **{
            f"{s}_{k}": v for s, d in rows.items() for k, v in d.items()
        }})
        if passes:
            selected.append(cell_type)
    return selected, audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs="+", required=True)
    parser.add_argument("--out-dir", type=Path, default=PROJECT / "tables")
    args = parser.parse_args()

    frame = load_counts(args.inputs)
    studies = sorted(frame.study.unique())
    strata = surviving_strata(frame)
    selected, audit = select_cell_types(strata, studies)
    if not selected:
        raise SystemExit("No cell type met the inclusion thresholds")

    keep = frame.merge(
        strata[["study", "run", "cell_type"]], on=["study", "run", "cell_type"]
    )
    keep = keep[keep.cell_type.isin(selected)]

    rows = []
    for (study, cell_type), group in keep.groupby(["study", "cell_type"]):
        total = int(group.reads.sum())
        runs = int(group.run.nunique())
        for acceptor, sub in group.groupby("acceptor"):
            rows.append(
                {
                    "study": study,
                    "cell_type": cell_type,
                    "acceptor": acceptor,
                    "reads": int(sub.reads.sum()),
                    "frac": sub.reads.sum() / total,
                    "total_reads": total,
                    "runs": runs,
                }
            )
    acceptor_table = pd.DataFrame(rows)

    per_run = strata[strata.cell_type.isin(selected)].copy()
    per_run["RO_frac"] = per_run.skip / per_run.n
    per_run = per_run.rename(columns={"skip": "RO"})[
        ["study", "run", "cell_type", "RO", "n", "RO_frac"]
    ]

    # Row order for panel e: pooled skip-all fraction, so the ordering is a
    # property of the data rather than of whichever study is listed first.
    skip_rows = acceptor_table[acceptor_table.acceptor.eq(SKIP)]
    totals = acceptor_table.groupby(["study", "cell_type"]).total_reads.first()
    pooled = {
        cell_type: skip_rows[skip_rows.cell_type.eq(cell_type)].reads.sum()
        / totals[totals.index.get_level_values("cell_type") == cell_type].sum()
        for cell_type in selected
    }
    order = sorted(selected, key=lambda c: pooled[c])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    acceptor_path = args.out_dir / "figdata_cd45_acceptor_usage.tsv"
    per_run_path = args.out_dir / "figdata_cd45_ro_per_run.tsv"
    manifest_path = args.out_dir / "figdata_cd45_inclusion_manifest.json"
    acceptor_table.to_csv(acceptor_path, sep="\t", index=False)
    per_run.to_csv(per_run_path, sep="\t", index=False)
    manifest = {
        "thresholds": {
            "min_stratum_reads": MIN_STRATUM_READS,
            "min_cell_type_reads": MIN_CELL_TYPE_READS,
            "min_cell_type_runs": MIN_CELL_TYPE_RUNS,
            "applied_to": "strata surviving min_stratum_reads, in every study",
        },
        "cell_type_aliases": CELL_TYPE_ALIASES,
        "studies": studies,
        "selected_cell_types": order,
        "strata": int(len(per_run)),
        "informative_reads": int(acceptor_table.groupby(
            ["study", "cell_type"]).total_reads.first().sum()),
        "audit": audit,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    print(json.dumps({k: v for k, v in manifest.items() if k != "audit"}, indent=2))
    print("\nper cell type (reads / runs per study, after the stratum filter):")
    for row in audit:
        mark = "SHOWN " if row["selected"] else "  --  "
        detail = "  ".join(
            f"{s}: {row[f'{s}_reads']:>7,d}r/{row[f'{s}_runs']:>2d}runs" for s in studies
        )
        print(f"  {mark} {row['cell_type']:14s} {detail}")
    for path in (acceptor_path, per_run_path, manifest_path):
        print("wrote", path)


if __name__ == "__main__":
    main()
