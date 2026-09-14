#!/usr/bin/env python3
"""Targeted, study-stratified replication audit for featured marrow genes.

Each study is processed independently with the same biological-unit collapse,
depth filter, shared-lineage restriction and within-source permutation used by
the corrected pooled analysis.  In addition to targeted empirical P values,
the script compares signed cell-type contrasts across studies.  These targeted
P values are not genome-wide FDR values.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import run_celltype_permutation as core


STUDIES = ("GSE276974", "GSE307660")
DEFAULT_TARGETS = ("PTPRC", "MS4A1")
GENE_IDS = {
    "PTPRC": {"diu": "ENSG00000081237", "apa": "ENSG00000081237", "ase": "PTPRC"},
    "MS4A1": {"diu": "ENSG00000156738", "apa": "ENSG00000156738", "ase": "MS4A1"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--registry", type=Path, default=Path(core.DEFAULT_REGISTRY))
    parser.add_argument("--studies", nargs="+", default=list(STUDIES))
    parser.add_argument("--targets", nargs="+", default=list(DEFAULT_TARGETS))
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--min-donors", type=int, default=3)
    parser.add_argument("--contrast-threshold", type=float, default=0.02)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def analysis_args(args: argparse.Namespace, analysis: str, seed: int) -> SimpleNamespace:
    return SimpleNamespace(
        analysis=analysis,
        min_depth=args.min_depth,
        min_donors=args.min_donors,
        negative_control=False,
        include_private_celltypes=False,
        seed=seed,
    )


def usage_profiles(
    sub: pd.DataFrame,
    feature_col: str,
    study: str,
    analysis: str,
    gene_name: str,
) -> pd.DataFrame:
    pivot = sub.pivot_table(
        index=["donor", "cell_type"],
        columns=feature_col,
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    usage = pivot.div(pivot.sum(axis=1), axis=0)
    equal_source = usage.groupby(level="cell_type").mean()
    long = (
        equal_source.rename_axis(index="cell_type", columns="feature")
        .stack(future_stack=True)
        .rename("mean_usage")
        .reset_index()
    )
    long.insert(0, "gene", gene_name)
    long.insert(0, "analysis", analysis)
    long.insert(0, "study", study)
    return long


def ase_profiles(
    sub: pd.DataFrame,
    study: str,
    gene_name: str,
) -> pd.DataFrame:
    current = sub.copy()
    current["fraction"] = current.hapA / current["count"]
    equal_source = current.groupby("cell_type")["fraction"].mean()
    long = equal_source.rename("mean_usage").reset_index()
    long["feature"] = "hapA_fraction"
    long.insert(0, "gene", gene_name)
    long.insert(0, "analysis", "ase")
    long.insert(0, "study", study)
    return long


def signed_contrast_vector(
    profile: pd.DataFrame,
    cell_types: list[str],
    features: list[str],
) -> np.ndarray:
    matrix = (
        profile.pivot(index="cell_type", columns="feature", values="mean_usage")
        .reindex(index=cell_types, columns=features, fill_value=0.0)
        .fillna(0.0)
    )
    rows: list[np.ndarray] = []
    for left_index, left in enumerate(cell_types):
        for right in cell_types[left_index + 1 :]:
            rows.append(
                matrix.loc[left].to_numpy(dtype=float)
                - matrix.loc[right].to_numpy(dtype=float)
            )
    return np.concatenate(rows) if rows else np.array([], dtype=float)


def correlation(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    if len(left) < 2 or np.std(left) == 0 or np.std(right) == 0:
        return np.nan, np.nan
    pearson = float(np.corrcoef(left, right)[0, 1])
    left_rank = pd.Series(left).rank(method="average").to_numpy()
    right_rank = pd.Series(right).rank(method="average").to_numpy()
    spearman = float(np.corrcoef(left_rank, right_rank)[0, 1])
    return pearson, spearman


def concordance_rows(
    profiles: pd.DataFrame,
    studies: list[str],
    threshold: float,
) -> list[dict]:
    if len(studies) != 2:
        raise ValueError("cross-study concordance currently requires exactly two studies")
    rows: list[dict] = []
    for (analysis, gene), group in profiles.groupby(["analysis", "gene"]):
        by_study = {study: group[group.study.eq(study)] for study in studies}
        common_ct = sorted(
            set(by_study[studies[0]].cell_type)
            & set(by_study[studies[1]].cell_type)
        )
        features = sorted(set(group.feature))
        left = signed_contrast_vector(by_study[studies[0]], common_ct, features)
        right = signed_contrast_vector(by_study[studies[1]], common_ct, features)
        pearson, spearman = correlation(left, right)
        informative = (np.abs(left) >= threshold) | (np.abs(right) >= threshold)
        sign_agreement = (
            float(np.mean(np.sign(left[informative]) == np.sign(right[informative])))
            if informative.any()
            else np.nan
        )
        cosine_denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        cosine = (
            float(np.dot(left, right) / cosine_denominator)
            if cosine_denominator > 0
            else np.nan
        )
        rows.append(
            {
                "analysis": analysis,
                "gene": gene,
                "study_left": studies[0],
                "study_right": studies[1],
                "n_common_celltypes": len(common_ct),
                "common_celltypes": ";".join(common_ct),
                "n_union_features": len(features),
                "n_signed_contrasts": len(left),
                "contrast_threshold": threshold,
                "n_informative_contrasts": int(informative.sum()),
                "pearson_signed_contrasts": pearson,
                "spearman_signed_contrasts": spearman,
                "cosine_signed_contrasts": cosine,
                "sign_agreement_fraction": sign_agreement,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    if args.permutations < 1:
        raise ValueError("--permutations must be positive")
    unknown = sorted(set(args.targets) - set(GENE_IDS))
    if unknown:
        raise ValueError(f"unknown targets: {unknown}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tests: list[dict] = []
    profile_frames: list[pd.DataFrame] = []
    input_paths: set[Path] = set()
    unit_audits: dict[str, dict] = {}

    for analysis in ("diu", "apa", "ase"):
        for study_index, study in enumerate(args.studies):
            current_args = analysis_args(
                args,
                analysis,
                args.seed + 1000 * study_index + 100 * ("diu", "apa", "ase").index(analysis),
            )
            registry = core.load_registry(args.registry, [study])
            raw, paths = core.read_inputs(core.DEFAULT_INPUTS[analysis], analysis, [study])
            input_paths.update(paths)
            bound, audit = core.bind_units(raw, registry, allow_missing=False)
            unit_audits[f"{analysis}:{study}"] = audit
            if analysis == "ase":
                prepared, gene_col, feature_col, _ = core.prepare_ase(bound, current_args)
            else:
                prepared, gene_col, feature_col, _ = core.prepare_usage(
                    bound, analysis, current_args
                )
            eligible = set(
                core.eligible_genes(
                    prepared, gene_col, feature_col, current_args.min_donors
                )
            )
            grouped = {str(key): frame for key, frame in prepared.groupby(gene_col)}
            for target in args.targets:
                gene = GENE_IDS[target][analysis]
                base = {
                    "analysis": analysis,
                    "study": study,
                    "gene_name": target,
                    "gene": gene,
                    "eligible": gene in eligible,
                    "targeted_permutations": args.permutations,
                    "genome_wide_fdr_claim": False,
                }
                if gene not in eligible:
                    tests.append(base)
                    continue
                sub = grouped[gene]
                rng = core.stable_rng(
                    current_args.seed,
                    f"targeted:{analysis}:{study}:{gene}",
                )
                if analysis == "ase":
                    result = core.test_ase_gene(sub, args.permutations, rng)
                    profile_frames.append(ase_profiles(sub, study, target))
                else:
                    result = core.test_usage_gene(
                        sub, str(feature_col), args.permutations, rng
                    )
                    profile_frames.append(
                        usage_profiles(
                            sub, str(feature_col), study, analysis, target
                        )
                    )
                tests.append({**base, **result})

    tests_frame = pd.DataFrame(tests)
    profiles = pd.concat(profile_frames, ignore_index=True)
    concordance = pd.DataFrame(
        concordance_rows(profiles, list(args.studies), args.contrast_threshold)
    )
    outputs = {
        "tests": args.output_dir / "study_stratified_targeted_tests.tsv",
        "profiles": args.output_dir / "study_stratified_profiles.tsv",
        "concordance": args.output_dir / "cross_study_concordance.tsv",
    }
    tests_frame.to_csv(outputs["tests"], sep="\t", index=False)
    profiles.to_csv(outputs["profiles"], sep="\t", index=False)
    concordance.to_csv(outputs["concordance"], sep="\t", index=False)
    manifest = {
        "schema_version": "scTHREAD.study_stratified_replication.v1",
        "studies": list(args.studies),
        "targets": list(args.targets),
        "analyses": ["diu", "apa", "ase"],
        "aggregation_unit": "study_id::donor_or_source_id",
        "targeted_permutations": args.permutations,
        "targeted_p_values_are_genome_wide_fdr": False,
        "min_depth": args.min_depth,
        "min_donors": args.min_donors,
        "shared_celltypes_only": True,
        "contrast_definition": (
            "signed differences between equal-source cell-type mean usage "
            "profiles for every cell-type pair and union feature"
        ),
        "contrast_threshold": args.contrast_threshold,
        "unit_audits": unit_audits,
        "registry": str(args.registry.resolve()),
        "registry_sha256": sha256(args.registry),
        "input_file_count": len(input_paths),
        "outputs": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in outputs.items()
        },
    }
    manifest_path = args.output_dir / "study_stratified_replication_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
