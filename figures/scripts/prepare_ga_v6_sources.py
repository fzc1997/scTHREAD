#!/usr/bin/env python3
"""Freeze registry tissue examples and real PTPRC isoform UMAP points for GA v6."""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW = ROOT / "figures/ga_gpt2_vector_v6/raw_api"
DEFAULT_OUT = ROOT / "figures/ga_gpt2_vector_v6"
REGISTRY_FROZEN = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
RAW_HASHES = {
    "ENST00000367364.json":
        "fbf4086d97602b183b96119e6657598f2391a98f26ae96501ac3b8160a08a086",
    "ENST00000697630.json":
        "e6c60b580817af64057f8b9db811b017c4c4e4ba6f5526c54d4356a4c3c8bdfb",
}
TISSUE_EVIDENCE = (
    ("Marrow", "Blood/immune", "GSE276974", ("marrow", "骨髓")),
    ("Myeloma", "Blood/immune", "GSE307660", ("myeloma", "骨髓瘤")),
    (
        "Frontal cortex",
        "Neural/sensory",
        "GSE178175",
        ("cortex", "皮层"),
    ),
    ("Retina", "Neural/sensory", "GSE255520", ("retina", "视网膜")),
    ("Heart", "Heart/vascular", "GSE288222", ("heart", "心脏")),
    ("Prostate", "Cancer", "GSE289790", ("prostate", "前列腺")),
    ("Glioma", "Cancer", "GSE301658", ("glioma", "胶质瘤")),
    ("Ovary", "Reproductive", "GSE248118", ("ovarian", "ovary", "卵巢")),
    ("Brain", "Neural/sensory", "GSE314176", ("brain", "脑")),
    (
        "Fibroblast→HSC",
        "Development/embryo",
        "GSE283658",
        ("reprogramming", "重编程", "hsc"),
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_rank(value: str) -> str:
    return hashlib.blake2b(value.encode(), digest_size=12).hexdigest()


def load_api(raw_dir: Path, transcript: str) -> tuple[pd.DataFrame, dict]:
    path = raw_dir / f"{transcript}.json"
    expected = RAW_HASHES[path.name]
    observed = sha256(path)
    if observed != expected:
        raise AssertionError(f"Raw API changed: {path} {observed} != {expected}")
    payload = json.loads(path.read_text())
    assert payload["species"] == "human"
    assert payload["gene"]["gname"] == "PTPRC"
    assert payload["signal"] == "isoform"
    assert payload["transcript"] == transcript
    assert payload["count"] == len(payload["items"]) == 71_913
    frame = pd.DataFrame(payload["items"])[
        ["cell_id", "cell_type", "umap1", "umap2", "expression"]
    ]
    return frame, payload


def freeze_umap(raw_dir: Path, outdir: Path) -> dict[str, object]:
    first, _ = load_api(raw_dir, "ENST00000367364")
    second, _ = load_api(raw_dir, "ENST00000697630")
    first = first.rename(columns={"expression": "expr_67364"})
    second = second.rename(columns={"expression": "expr_97630"})
    merged = first.merge(
        second,
        on=["cell_id", "cell_type", "umap1", "umap2"],
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != 71_913:
        raise AssertionError(f"Coordinate-matched UMAP count changed: {len(merged)}")

    positive = (merged["expr_67364"] > 0) | (merged["expr_97630"] > 0)
    selected_positive = merged[positive].copy()
    zero = merged[~positive].copy()
    zero["stable_rank"] = zero["cell_id"].map(stable_rank)
    selected_zero = (
        zero.sort_values(["cell_type", "stable_rank"])
        .groupby("cell_type", observed=True, group_keys=False)
        .head(80)
        .drop(columns="stable_rank")
        .copy()
    )
    selected_positive["sample_reason"] = "positive_either_isoform"
    selected_zero["sample_reason"] = "double_zero_background"
    sample = pd.concat([selected_zero, selected_positive], ignore_index=True)
    sample = sample.sort_values(
        ["sample_reason", "cell_type", "cell_id"],
        kind="mergesort",
    ).reset_index(drop=True)
    if sample["cell_id"].duplicated().any():
        raise AssertionError("Duplicate UMAP cell in frozen sample")
    if int((sample["sample_reason"] == "positive_either_isoform").sum()) != 8_432:
        raise AssertionError("Positive-cell union changed")

    sample_path = outdir / "ptprc_two_isoform_umap_points.tsv"
    sample.to_csv(sample_path, sep="\t", index=False)
    stats = {
        "embedding_cells": len(merged),
        "retained_cells": len(sample),
        "positive_either_isoform": int(positive.sum()),
        "double_positive": int(
            ((merged["expr_67364"] > 0) & (merged["expr_97630"] > 0)).sum()
        ),
        "background_cells_sampled": len(selected_zero),
        "ENST00000367364": {
            "positive_cells": int((merged["expr_67364"] > 0).sum()),
            "molecules_in_embedding": int(merged["expr_67364"].sum()),
            "display_cap": 5,
        },
        "ENST00000697630": {
            "positive_cells": int((merged["expr_97630"] > 0).sum()),
            "molecules_in_embedding": int(merged["expr_97630"].sum()),
            "display_cap": 3,
        },
        "raw_api_sha256": RAW_HASHES,
        "sample_tsv_sha256": sha256(sample_path),
    }
    (outdir / "ptprc_two_isoform_umap_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True) + "\n"
    )
    return stats


def freeze_tissues(outdir: Path) -> pd.DataFrame:
    registry = pd.read_csv(REGISTRY_FROZEN, sep="\t", dtype=str)
    registry["isoquant_cells_num"] = pd.to_numeric(
        registry["isoquant_cells"], errors="coerce"
    ).fillna(0)
    done = registry[
        registry["isoquant_status"].fillna("").str.lower().eq("done")
    ].copy()
    if len(done) != 469 or int(done["isoquant_cells_num"].sum()) != 845_781:
        raise AssertionError("Frozen registry no longer matches 469 / 845,781")

    rows: list[dict[str, object]] = []
    candidate_columns = [
        "description",
        "note",
        "annotation_source",
        "biology_group",
    ]
    searchable_columns = [
        column for column in candidate_columns if column in done.columns
    ]
    if not searchable_columns:
        raise AssertionError("Frozen registry has no searchable evidence columns")
    for display, system, gse, terms in TISSUE_EVIDENCE:
        study = done[done["gse"].eq(gse)].copy()
        if study.empty:
            raise AssertionError(f"Registry study missing for tissue {display}: {gse}")
        evidence = " | ".join(
            sorted(
                {
                    str(value)
                    for column in searchable_columns
                    for value in study[column].dropna()
                    if str(value).strip()
                }
            )
        )
        if not any(term.lower() in evidence.lower() for term in terms):
            raise AssertionError(
                f"Registry does not support tissue label {display}: {evidence}"
            )
        rows.append(
            {
                "display": display,
                "system": system,
                "source_gse": gse,
                "samples": int(study["srr"].nunique()),
                "cells": int(study["isoquant_cells_num"].sum()),
                "registry_evidence": evidence,
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(outdir / "registry_tissue_examples.tsv", sep="\t", index=False)
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    tissues = freeze_tissues(args.outdir)
    stats = freeze_umap(args.raw_dir, args.outdir)
    print("REGISTRY TISSUES")
    print(tissues.to_csv(sep="\t", index=False))
    print("UMAP STATS", stats)
    print("GA V6 SOURCES PASS")


if __name__ == "__main__":
    main()
