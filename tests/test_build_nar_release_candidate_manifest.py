from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))
import build_nar_release_candidate_manifest as builder  # noqa: E402


def _row(run: str, study: str, species: str, cells: int) -> dict:
    return {
        "srr": run,
        "gse": study,
        "species": species,
        "platform": "ONT",
        "biology_group": "test",
        "description": "historical",
        "pipeline": "test",
        "isoquant_status": "done",
        "isoquant_cells": cells,
        "isoquant_cells_method": "test",
        "isoquant_transcripts": 10,
        "annotation_status": "done",
        "annotation_method": "test",
        "annotation_cells": cells,
    }


def test_exact_run_metadata_resolution(tmp_path: Path, monkeypatch) -> None:
    frozen = pd.DataFrame(
        [_row("RUN1", "STUDY1", "human", 3), _row("RUN2", "STUDY2", "?", 4)]
    )
    current = frozen.copy()
    current.loc[current["srr"] == "RUN2", ["species", "description"]] = [
        "mouse",
        "resolved",
    ]
    frozen_path = tmp_path / "frozen.tsv"
    current_path = tmp_path / "current.tsv"
    out = tmp_path / "out"
    frozen.to_csv(frozen_path, sep="\t", index=False)
    current.to_csv(current_path, sep="\t", index=False)

    monkeypatch.setattr(
        builder,
        "EXPECTED",
        {"runs": 2, "studies": 2, "cells": 7, "historical_unresolved_species": 1},
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build",
            "--frozen-registry",
            str(frozen_path),
            "--current-registry",
            str(current_path),
            "--output-dir",
            str(out),
        ],
    )
    builder.main()

    manifest = pd.read_csv(out / "release_candidate_manifest.tsv", sep="\t")
    assert manifest.set_index("srr").loc["RUN2", "species"] == "mouse"
    assert manifest["isoquant_cells"].sum() == 7
    summary = json.loads((out / "release_candidate_summary.json").read_text())
    assert summary["counts"]["runs"] == 2
    assert summary["metadata_resolution_registry"]["resolved_run_count"] == 1


def test_resolution_rejects_cell_count_change(tmp_path: Path, monkeypatch) -> None:
    frozen = pd.DataFrame([_row("RUN1", "STUDY1", "?", 3)])
    current = pd.DataFrame([_row("RUN1", "STUDY1", "mouse", 4)])
    frozen_path = tmp_path / "frozen.tsv"
    current_path = tmp_path / "current.tsv"
    frozen.to_csv(frozen_path, sep="\t", index=False)
    current.to_csv(current_path, sep="\t", index=False)
    monkeypatch.setattr(
        builder,
        "EXPECTED",
        {"runs": 1, "studies": 1, "cells": 3, "historical_unresolved_species": 1},
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "build",
            "--frozen-registry",
            str(frozen_path),
            "--current-registry",
            str(current_path),
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )
    with pytest.raises(ValueError, match="cell count changed"):
        builder.main()
