#!/usr/bin/env python3
"""Reduce one ANCHOR signed h5ad to per-cell-type allelic gene counts."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import duckdb
import numpy as np
import pandas as pd
import scipy.sparse as sp

from species_resources import ANNOTATION, annotation_pattern, run_config

SIGNED_ROOT = Path(os.environ.get("SCTHREAD_SIGNED_ASE_ROOT", "."))
DEFAULT_OUTPUT = Path(os.environ.get("SCTHREAD_ASE_OUTPUT", "outputs/f2_grammar/agg_ase_ct"))


def normalize_barcode(values) -> pd.Series:
    """Use the release annotation namespace: remove a suffix starting at '-'."""
    return pd.Series(np.asarray(values).astype(str), dtype="string").str.split(
        "-", n=1
    ).str[0]


def vector_sum(matrix, mask: np.ndarray) -> np.ndarray:
    return np.asarray(matrix[mask].sum(axis=0)).ravel()


def expressed_cells(hap_a, hap_b, mask: np.ndarray) -> np.ndarray:
    subset = hap_a[mask] + hap_b[mask]
    if sp.issparse(subset):
        return np.asarray(subset.getnnz(axis=0)).ravel()
    return np.count_nonzero(np.asarray(subset), axis=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    registry_row, _ = run_config(args.run)
    gse = registry_row["gse"]
    signed_h5ad = SIGNED_ROOT / args.run / "cellgene" / f"{args.run}_signed.h5ad"
    if not signed_h5ad.is_file():
        raise SystemExit(f"MISSING {signed_h5ad}")

    pattern = annotation_pattern(args.run, registry_row)
    labels = duckdb.connect().execute(
        """
        SELECT barcode_norm, active_label, evidence_tier
        FROM read_parquet(?)
        WHERE gse = ? AND source_cell_id LIKE ?
          AND active_label NOT IN ('Unlabeled', 'Unassigned', '')
        """,
        [str(ANNOTATION), gse, pattern],
    ).df()
    if labels.empty:
        raise SystemExit(
            f"NO_ANNOTATION_ROWS {args.run} {gse} pattern={pattern}"
        )
    labels["barcode_norm"] = normalize_barcode(labels["barcode_norm"]).values
    labels["evidence_tier"] = pd.to_numeric(
        labels["evidence_tier"], errors="coerce"
    ).fillna(-1)
    labels = labels.sort_values(
        "evidence_tier", ascending=False
    ).drop_duplicates("barcode_norm")
    barcode_to_type = dict(zip(labels.barcode_norm, labels.active_label))

    atlas = ad.read_h5ad(signed_h5ad)
    obs_barcode = normalize_barcode(atlas.obs_names)
    cell_types = obs_barcode.map(barcode_to_type).fillna("Unlabeled").to_numpy()
    labeled_cells = int(np.count_nonzero(cell_types != "Unlabeled"))
    if labeled_cells < 5:
        raise SystemExit(
            f"INSUFFICIENT_BARCODE_OVERLAP {args.run} "
            f"labeled={labeled_cells}/{atlas.n_obs} pattern={pattern}"
        )

    hap_a = atlas.layers["count_b6"]
    hap_b = atlas.layers["count_dba"]
    genes = np.asarray(atlas.var_names.astype(str))
    frames: list[pd.DataFrame] = []
    for cell_type in pd.unique(cell_types):
        if cell_type == "Unlabeled":
            continue
        mask = cell_types == cell_type
        n_cells = int(mask.sum())
        if n_cells < 5:
            continue
        count_a = vector_sum(hap_a, mask)
        count_b = vector_sum(hap_b, mask)
        n_cells_expr = expressed_cells(hap_a, hap_b, mask)
        keep = (count_a + count_b) > 0
        if np.any(keep):
            frames.append(
                pd.DataFrame(
                    {
                        "gene": genes[keep],
                        "cell_type": cell_type,
                        "hapA": count_a[keep].astype(float),
                        "hapB": count_b[keep].astype(float),
                        "n_cells": n_cells,
                        "n_cells_expr": n_cells_expr[keep].astype(int),
                    }
                )
            )
    if not frames:
        raise SystemExit(f"NO_ASE_OUTPUT_ROWS {args.run}")

    result = pd.concat(frames, ignore_index=True)
    result["run"] = args.run
    result["gse"] = gse
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / f"{args.run}.ase_ct.parquet"
    partial = target.with_suffix(target.suffix + f".partial.{os.getpid()}")
    result.to_parquet(partial, index=False)
    os.replace(partial, target)
    genes_two_types = (
        result[result.n_cells_expr >= 10]
        .groupby("gene").cell_type.nunique().ge(2).sum()
    )
    print(
        f"[{args.run}] {gse} cells={atlas.n_obs} labeled={labeled_cells} "
        f"genes={atlas.n_vars} rows={len(result):,} "
        f"cell_types={result.cell_type.nunique()} genes_in_2ct={genes_two_types:,} "
        f"pattern={pattern} output={target}"
    )


if __name__ == "__main__":
    main()
