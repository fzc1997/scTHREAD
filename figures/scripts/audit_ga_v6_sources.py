#!/usr/bin/env python3
"""Audit registry breadth, PTPRC isoform UMAP APIs and ASE summary for GA v6."""

from __future__ import annotations

import os

import json
from pathlib import Path

import pandas as pd

import render_ga_gpt2_real_components as components


ROOT = Path(__file__).resolve().parents[2]
REGISTRY = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv")
REGISTRY_FROZEN = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
RAW_API = ROOT / "figures/ga_gpt2_vector_v6/raw_api"
ASE = Path(
    os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD") + "/results/paper1/f2_grammar/"
    "figdata/ase_interaction.tsv"
)


def registry_summary(path: Path, label: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str)
    table["isoquant_cells_num"] = pd.to_numeric(
        table.get("isoquant_cells"), errors="coerce"
    ).fillna(0)
    done = table[
        table.get("isoquant_status", "").fillna("").str.lower().eq("done")
    ].copy()
    current = table[
        table.get("in_scthread_current", "").fillna("").str.lower().eq("yes")
    ].copy()
    print(
        "REGISTRY",
        label,
        {
            "rows": len(table),
            "done_rows": len(done),
            "done_cells": int(done["isoquant_cells_num"].sum()),
            "done_studies": int(done["gse"].nunique()),
            "current_rows": len(current),
            "current_studies": int(current["gse"].nunique()),
        },
    )
    return done


def audit_registry() -> None:
    frozen = registry_summary(REGISTRY_FROZEN, "frozen")
    live = registry_summary(REGISTRY, "live")
    if len(frozen) != 469 or int(frozen["isoquant_cells_num"].sum()) != 845_781:
        raise AssertionError("Frozen registry no longer matches 469 / 845,781")
    frozen["system"] = frozen.apply(components.classify_system, axis=1)
    studies = (
        frozen.groupby(["gse", "system"], as_index=False)
        .agg(
            samples=("srr", "nunique"),
            cells=("isoquant_cells_num", "sum"),
            species=("species", lambda x: "|".join(sorted(set(x.dropna())))),
            biology_group=(
                "biology_group",
                lambda x: "|".join(sorted(v for v in set(x.dropna()) if v)),
            ),
            descriptions=(
                "description",
                lambda x: "|".join(sorted(v for v in set(x.dropna()) if v)),
            ),
            notes=(
                "note",
                lambda x: "|".join(sorted(v for v in set(x.dropna()) if v)),
            ),
        )
        .sort_values(["system", "cells"], ascending=[True, False])
    )
    print("FROZEN STUDY COVERAGE")
    print(studies.to_csv(sep="\t", index=False))

    frozen_keys = set(frozen["srr"].astype(str))
    additions = live[~live["srr"].astype(str).isin(frozen_keys)].copy()
    if not additions.empty:
        columns = [
            "srr",
            "gse",
            "species",
            "platform",
            "biology_group",
            "description",
            "note",
            "isoquant_status",
            "isoquant_cells_num",
            "in_scthread_current",
        ]
        print("LIVE ADDITIONS SINCE FREEZE")
        print(additions[columns].to_csv(sep="\t", index=False))


def load_api(transcript: str) -> pd.DataFrame:
    path = RAW_API / f"{transcript}.json"
    payload = json.loads(path.read_text())
    assert payload["species"] == "human"
    assert payload["gene"]["gname"] == "PTPRC"
    assert payload["signal"] == "isoform"
    assert payload["transcript"] == transcript
    assert payload["count"] == len(payload["items"]) == 71_913
    frame = pd.DataFrame(payload["items"])
    required = {"cell_id", "cell_type", "umap1", "umap2", "expression"}
    if not required.issubset(frame.columns):
        raise AssertionError(f"Missing API columns: {sorted(required - set(frame))}")
    return frame[list(required)]


def audit_umap() -> None:
    first = load_api("ENST00000367364").rename(
        columns={"expression": "expr_67364"}
    )
    second = load_api("ENST00000697630").rename(
        columns={"expression": "expr_97630"}
    )
    merged = first.merge(
        second,
        on=["cell_id", "cell_type", "umap1", "umap2"],
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != 71_913:
        raise AssertionError(f"Matched UMAP coordinates changed: {len(merged)}")
    for column in ("expr_67364", "expr_97630"):
        values = merged[column]
        print(
            "UMAP EXPRESSION",
            column,
            {
                "cells": len(values),
                "positive_cells": int((values > 0).sum()),
                "sum_molecules": int(values.sum()),
                "max": float(values.max()),
                "q95_positive": float(values[values > 0].quantile(0.95)),
                "q99_positive": float(values[values > 0].quantile(0.99)),
            },
        )
    union = (merged["expr_67364"] > 0) | (merged["expr_97630"] > 0)
    high = (merged[["expr_67364", "expr_97630"]].max(axis=1) >= 2)
    print(
        "UMAP UNION",
        {
            "positive_either": int(union.sum()),
            "both_positive": int(
                ((merged["expr_67364"] > 0) & (merged["expr_97630"] > 0)).sum()
            ),
            "max_expression_ge_2": int(high.sum()),
            "zero_both": int((~union).sum()),
        },
    )


def audit_ase() -> None:
    table = pd.read_csv(ASE, sep="\t")
    print(
        "ASE",
        {
            "genes_tested": len(table),
            "sig_true": int(table["sig"].sum()),
            "raw_p_lt_0.05": int((table["pval"] < 0.05).sum()),
            "raw_p_fraction": float((table["pval"] < 0.05).mean()),
            "q_lt_0.05": int((table["qval"] < 0.05).sum()),
            "min_q": float(table["qval"].min()),
        },
    )


def main() -> None:
    audit_registry()
    audit_umap()
    audit_ase()
    print("GA V6 SOURCE AUDIT PASS")


if __name__ == "__main__":
    main()
