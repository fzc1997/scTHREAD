#!/usr/bin/env python3
"""Summarize evidence-completeness tiers in the frozen NAR run manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = (
    ROOT
    / "tables/release_20260803/release_manifest.tsv"
)
DEFAULT_OUTPUT = ROOT / "tables/Table_S_release_completeness_tiers.tsv"


def summarize(frame: pd.DataFrame, tier: str) -> dict[str, object]:
    transcript_documented = frame["isoquant_transcripts"].notna()
    if tier == "transcript_count_documented":
        subset = frame.loc[transcript_documented]
        definition = "isoquant_transcripts is recorded in the frozen release manifest"
        interpretation = (
            "Run membership and a run-level transcript count are documented; "
            "a complete immutable feature bundle is still pending."
        )
    else:
        subset = frame.loc[~transcript_documented]
        definition = "isoquant_transcripts is missing in the frozen release manifest"
        interpretation = (
            "One-well run membership and one called cell per row are documented; "
            "these rows must not be described as transcript-feature-complete."
        )
    return {
        "tier": tier,
        "definition": definition,
        "n_runs": subset["srr"].nunique(),
        "n_studies": subset["gse"].nunique(),
        "n_pipeline_called_cells": int(subset["isoquant_cells"].sum()),
        "human_runs": int((subset["species"] == "human").sum()),
        "mouse_runs": int((subset["species"] == "mouse").sum()),
        "ont_runs": int((subset["platform"] == "ONT").sum()),
        "pacbio_runs": int((subset["platform"] == "PacBio").sum()),
        "pipeline_values": "|".join(sorted(subset["pipeline"].dropna().unique())),
        "cell_count_methods": "|".join(
            sorted(subset["isoquant_cells_method"].dropna().unique())
        ),
        "interpretation": interpretation,
    }


def build(manifest: Path, output: Path) -> None:
    frame = pd.read_csv(manifest, sep="\t", keep_default_na=True)
    if (
        len(frame) != 453
        or frame["srr"].nunique() != 453
        or frame["gse"].nunique() != 34
    ):
        raise RuntimeError(
            f"frozen release is 453 runs / 34 studies, manifest gives "
            f"{frame['srr'].nunique()} / {frame['gse'].nunique()}"
        )

    missing = frame["isoquant_transcripts"].isna()
    if missing.sum() != 298:
        raise RuntimeError("Expected 298 membership-only rows")
    membership_only = frame.loc[missing]
    if not (
        membership_only["pipeline"].value_counts().to_dict()
        == {"EXCLUDE": 283, "unclassified": 15}
        and (membership_only["isoquant_cells_method"] == "smartseq2_one_well").all()
        and (membership_only["isoquant_cells"] == 1).all()
    ):
        raise RuntimeError("Membership-only tier no longer matches its audited definition")

    result = pd.DataFrame(
        [
            summarize(frame, "transcript_count_documented"),
            summarize(frame, "one_well_membership_only"),
        ]
    )
    if int(result["n_runs"].sum()) != len(frame) or int(
        result["n_pipeline_called_cells"].sum()
    ) != int(frame["isoquant_cells"].sum()):
        raise RuntimeError("Completeness tiers do not partition the frozen release")
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    print(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.manifest.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
