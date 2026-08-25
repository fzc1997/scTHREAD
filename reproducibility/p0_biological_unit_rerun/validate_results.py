#!/usr/bin/env python3
"""Independent validation and compact summary of the P0 rerun outputs."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "."))
TABLES = ROOT / "tables/p0_biological_unit_rerun"
ANALYSES = ("diu", "apa", "ase")
NULL_SEEDS = (101, 202, 303)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def independent_bh(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg from its rank definition, without statsmodels."""
    count = len(pvalues)
    order = np.argsort(pvalues, kind="mergesort")
    ranked = pvalues[order] * count / np.arange(1, count + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(count, dtype=float)
    result[order] = np.minimum(ranked, 1.0)
    return result


def validate_one(path: Path, expected_null: bool) -> tuple[pd.DataFrame, dict]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    manifest_path = path.with_suffix(path.suffix + ".manifest.json")
    manifest = json.loads(manifest_path.read_text())
    if bool(manifest["negative_control"]) != expected_null:
        raise AssertionError(f"negative-control identity mismatch: {path}")
    if manifest["output_sha256"] != sha256(path):
        raise AssertionError(f"output hash mismatch: {path}")
    data = pd.read_csv(path, sep="\t")
    required = {
        "gene",
        "pval",
        "qval",
        "effect_equal_donor",
        "effect_pooled_counts",
        "n_celltypes",
        "n_donors",
        "min_donors_per_celltype",
        "sig",
    }
    if required - set(data.columns):
        raise AssertionError(f"missing result columns: {path}")
    if data.gene.duplicated().any():
        raise AssertionError(f"duplicate genes: {path}")
    if not data.pval.between(0, 1).all() or not data.qval.between(0, 1).all():
        raise AssertionError(f"invalid P/Q range: {path}")
    q_independent = independent_bh(data.pval.to_numpy(float))
    max_q_error = float(np.max(np.abs(q_independent - data.qval.to_numpy(float))))
    if max_q_error > 1e-12:
        raise AssertionError(f"independent BH mismatch {max_q_error}: {path}")
    independent_sig = (q_independent < 0.05) & (
        data.effect_equal_donor.to_numpy(float) >= 0.20
    )
    parsed_sig = data.sig.astype(str).str.lower().eq("true").to_numpy()
    if not np.array_equal(independent_sig, parsed_sig):
        raise AssertionError(f"significance flag mismatch: {path}")
    audit = manifest["unit_audit"]
    if audit["expected_registry_runs"] != 25:
        raise AssertionError(f"registry run closure mismatch: {path}")
    if audit["available_registry_runs"] != 25:
        raise AssertionError(f"input coverage mismatch: {path}")
    if audit["missing_registry_runs"]:
        raise AssertionError(f"missing registry runs: {path}")
    if manifest["aggregation_unit"] != "study_id::donor_or_source_id":
        raise AssertionError(f"aggregation-unit mismatch: {path}")
    return data, {
        "path": str(path),
        "sha256": sha256(path),
        "rows": int(len(data)),
        "significant": int(independent_sig.sum()),
        "raw_p_lt_0_05_fraction": float(data.pval.lt(0.05).mean()),
        "max_independent_bh_error": max_q_error,
        "inner_permutations": int(manifest["permutations"]),
        "permutation_seed": int(manifest["seed"]),
    }


def main() -> None:
    rows = []
    detailed = {}
    for analysis in ANALYSES:
        observed_path = TABLES / f"{analysis}_observed_9999.tsv"
        observed, audit = validate_one(observed_path, expected_null=False)
        audit.update({"analysis": analysis, "dataset": "observed", "seed": ""})
        rows.append(audit)
        requested = [
            "PTPRC",
            "MS4A1",
            "ENSG00000081237",
            "ENSG00000156738",
        ]
        detailed[analysis] = {
            "top10": observed.head(10).to_dict(orient="records"),
            "requested_genes": observed[observed.gene.isin(requested)].to_dict(
                orient="records"
            ),
            "requested_gene_identifiers": {
                "PTPRC": ["PTPRC", "ENSG00000081237"],
                "MS4A1": ["MS4A1", "ENSG00000156738"],
            },
        }
        for seed in NULL_SEEDS:
            null_path = TABLES / f"{analysis}_null_seed{seed}.tsv"
            _, null_audit = validate_one(null_path, expected_null=True)
            null_audit.update(
                {"analysis": analysis, "dataset": "null", "seed": seed}
            )
            rows.append(null_audit)
    summary = pd.DataFrame(rows)
    null = summary[summary.dataset.eq("null")]
    failures = null[
        ~null.raw_p_lt_0_05_fraction.between(0.02, 0.08)
        | null.significant.ne(0)
    ]
    if not failures.empty:
        raise AssertionError(
            "null diagnostic outside the recorded 0.02–0.08 raw-P range "
            "or contains FDR discoveries:\n"
            + failures.to_string(index=False)
        )
    summary_path = TABLES / "validation_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    report = {
        "schema_version": "scTHREAD.biological_unit_rerun.validation.v2",
        "status": "PASS",
        "independence": (
            "BH and significance flags recomputed without statsmodels; "
            "file hashes and unit/input contracts replayed from manifests"
        ),
        "null_recorded_range_raw_p_lt_0_05": [0.02, 0.08],
        "summary_path": str(summary_path),
        "summary_sha256": sha256(summary_path),
        "details": detailed,
    }
    report_path = TABLES / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
