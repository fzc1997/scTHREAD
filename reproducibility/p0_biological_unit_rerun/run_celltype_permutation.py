#!/usr/bin/env python3
"""Biological-unit-aware DIU, APA, or ASE cell-type permutation analysis.

The publication analysis unit is ``study_id::donor_or_source_id`` from the
validated scLong analysis-unit registry. Technical runs and repeated states
from the same source are summed before testing. Cell-type labels are permuted
within source, and only cell types observed in at least ``min_donors`` sources
are retained for a gene.

This script deliberately writes versioned outputs and never overwrites the
legacy Fig. 2 tables.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


DEFAULT_REGISTRY = os.environ.get(
    "SCTHREAD_ANALYSIS_UNIT_REGISTRY",
    "metadata/analysis_unit_manifest.tsv",
)
DEFAULT_INPUTS = {
    "diu": os.environ.get("SCTHREAD_DIU_INPUT", "results/paper1/f2_grammar/agg_diu/*.diu.parquet"),
    "apa": os.environ.get("SCTHREAD_APA_INPUT", "results/paper1/f2_grammar/agg_apa/*.apa.parquet"),
    "ase": os.environ.get("SCTHREAD_ASE_INPUT", "results/paper1/f2_grammar/agg_ase_ct/*.ase_ct.parquet"),
}
HARMONIZER_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(HARMONIZER_DIR))
from ct_harmonize import SHARED, harmonize  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", choices=("diu", "apa", "ase"), required=True)
    parser.add_argument("--registry", default=DEFAULT_REGISTRY)
    parser.add_argument("--input-glob")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--studies", nargs="+", default=("GSE276974", "GSE307660")
    )
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--min-donors", type=int, default=3)
    parser.add_argument("--effect-threshold", type=float, default=0.20)
    parser.add_argument("--max-genes", type=int)
    parser.add_argument("--negative-control", action="store_true")
    parser.add_argument("--include-private-celltypes", action="store_true")
    parser.add_argument("--allow-missing-registry-runs", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_rng(seed: int, token: str) -> np.random.RandomState:
    payload = f"{seed}:{token}".encode("utf-8")
    value = int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")
    return np.random.RandomState(value)


def load_registry(path: Path, studies: list[str]) -> pd.DataFrame:
    registry = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "study_id",
        "run_id",
        "analysis_unit_id",
        "donor_or_source_id",
        "included",
    }
    missing = required - set(registry.columns)
    if missing:
        raise ValueError(f"registry lacks columns: {sorted(missing)}")
    registry = registry[
        registry.study_id.isin(studies)
        & registry.included.str.lower().eq("true")
    ].copy()
    if registry.run_id.duplicated().any():
        dup = registry.loc[registry.run_id.duplicated(False), "run_id"].tolist()
        raise ValueError(f"duplicate registry run IDs: {dup}")
    if registry.donor_or_source_id.eq("").any():
        raise ValueError("blank donor_or_source_id in target registry rows")
    registry["donor"] = (
        registry.study_id.astype(str)
        + "::"
        + registry.donor_or_source_id.astype(str)
    )
    return registry


def read_inputs(pattern: str, analysis: str, studies: list[str]) -> tuple[pd.DataFrame, list[Path]]:
    files = [Path(x) for x in sorted(glob.glob(pattern))]
    if not files:
        raise FileNotFoundError(f"no inputs match {pattern}")
    columns = {
        "diu": ["gene_id", "transcript_id", "cell_type", "run", "gse", "count"],
        "apa": ["gene_id", "pas_id", "cell_type", "run", "gse", "n_molecules"],
        "ase": ["gene", "cell_type", "run", "gse", "hapA", "hapB"],
    }[analysis]
    frames = []
    for path in files:
        frame = pd.read_parquet(path, columns=columns)
        frame = frame[frame.gse.isin(studies)]
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise ValueError("target studies have no rows in the matched inputs")
    return pd.concat(frames, ignore_index=True), files


def bind_units(
    data: pd.DataFrame,
    registry: pd.DataFrame,
    allow_missing: bool,
) -> tuple[pd.DataFrame, dict]:
    available = set(data.run.astype(str))
    expected = set(registry.run_id.astype(str))
    excluded_present = sorted(available - expected)
    missing_expected = sorted(expected - available)
    if missing_expected and not allow_missing:
        raise ValueError(f"registry runs absent from inputs: {missing_expected}")
    data = data[data.run.astype(str).isin(expected)].copy()
    data = data.merge(
        registry[
            [
                "run_id",
                "study_id",
                "analysis_unit_id",
                "donor_or_source_id",
                "donor",
            ]
        ],
        left_on="run",
        right_on="run_id",
        how="left",
        validate="many_to_one",
    )
    if data.donor.isna().any():
        raise ValueError("unmapped target rows remain after registry join")
    audit = {
        "expected_registry_runs": len(expected),
        "available_registry_runs": len(expected & available),
        "excluded_input_runs_ignored": excluded_present,
        "missing_registry_runs": missing_expected,
        "analysis_units": int(registry.analysis_unit_id.nunique()),
        "independent_donors": int(registry.donor.nunique()),
    }
    return data, audit


def outer_null_labels(data: pd.DataFrame, seed: int) -> pd.Series:
    mapping: dict[tuple[str, str], str] = {}
    for donor, group in data[["donor", "cell_type"]].drop_duplicates().groupby("donor"):
        labels = sorted(group.cell_type.astype(str))
        shuffled = stable_rng(seed, f"outer:{donor}").permutation(labels)
        mapping.update({(donor, old): new for old, new in zip(labels, shuffled)})
    return pd.Series(
        [mapping[(d, c)] for d, c in zip(data.donor, data.cell_type)],
        index=data.index,
    )


def fixed_point_filter(
    sub: pd.DataFrame,
    gene_col: str,
    min_donors: int,
) -> pd.DataFrame:
    """Keep cell types replicated across donors and donors informative for CT."""
    current = sub
    while not current.empty:
        before = len(current)
        ct_n = current.groupby([gene_col, "cell_type"]).donor.nunique()
        keep_ct = ct_n[ct_n >= min_donors].index
        key = pd.MultiIndex.from_frame(current[[gene_col, "cell_type"]])
        current = current[key.isin(keep_ct)]
        donor_n = current.groupby([gene_col, "donor"]).cell_type.nunique()
        keep_donor = donor_n[donor_n >= 2].index
        key = pd.MultiIndex.from_frame(current[[gene_col, "donor"]])
        current = current[key.isin(keep_donor)]
        if len(current) == before:
            break
    return current.copy()


def prepare_usage(
    data: pd.DataFrame,
    analysis: str,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, str, str, str]:
    gene_col = "gene_id"
    feature_col = "transcript_id" if analysis == "diu" else "pas_id"
    if analysis == "apa":
        data = data.rename(columns={"n_molecules": "count"})
    data["cell_type"] = harmonize(data.cell_type)
    data = data[data.cell_type.ne("Unassigned")].copy()
    if not args.include_private_celltypes:
        data = data[data.cell_type.isin(SHARED)].copy()
    data = data.groupby(
        [gene_col, feature_col, "cell_type", "donor"], as_index=False
    )["count"].sum()
    if args.negative_control:
        data["cell_type"] = outer_null_labels(data, args.seed)
        data = data.groupby(
            [gene_col, feature_col, "cell_type", "donor"], as_index=False
        )["count"].sum()
    depth = data.groupby([gene_col, "cell_type", "donor"])["count"].transform("sum")
    data = data[depth >= args.min_depth].copy()
    data = fixed_point_filter(data, gene_col, args.min_donors)
    return data, gene_col, feature_col, "count"


def prepare_ase(
    data: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, str, None, str]:
    gene_col = "gene"
    data["cell_type"] = harmonize(data.cell_type)
    data = data[data.cell_type.ne("Unassigned")].copy()
    if not args.include_private_celltypes:
        data = data[data.cell_type.isin(SHARED)].copy()
    data = data.groupby(
        [gene_col, "cell_type", "donor"], as_index=False
    )[["hapA", "hapB"]].sum()
    if args.negative_control:
        data["cell_type"] = outer_null_labels(data, args.seed)
        data = data.groupby(
            [gene_col, "cell_type", "donor"], as_index=False
        )[["hapA", "hapB"]].sum()
    data["count"] = data.hapA + data.hapB
    data = data[data["count"] >= args.min_depth].copy()
    data = fixed_point_filter(data, gene_col, args.min_donors)
    return data, gene_col, None, "count"


def eligible_genes(
    data: pd.DataFrame,
    gene_col: str,
    feature_col: str | None,
    min_donors: int,
) -> list[str]:
    summary = data.groupby(gene_col).agg(
        n_celltypes=("cell_type", "nunique"),
        n_donors=("donor", "nunique"),
    )
    keep = summary[(summary.n_celltypes >= 2) & (summary.n_donors >= min_donors)]
    genes = set(keep.index.astype(str))
    if feature_col:
        n_feature = data.groupby(gene_col)[feature_col].nunique()
        genes &= set(n_feature[n_feature >= 2].index.astype(str))
    return sorted(genes)


def permutation_stat(
    values: np.ndarray,
    depth: np.ndarray,
    donor: np.ndarray,
    cell_type: np.ndarray,
    permutations: int,
    rng: np.random.RandomState,
) -> tuple[float, float]:
    """Depth-weighted donor-centred statistic and restricted empirical P."""
    cell_levels = sorted(set(cell_type))
    code_map = {label: index for index, label in enumerate(cell_levels)}
    codes = np.array([code_map[x] for x in cell_type], dtype=int)
    residual = values.copy()
    donor_groups = []
    for label in sorted(set(donor)):
        index = np.where(donor == label)[0]
        residual[index] = values[index] - np.average(
            values[index], axis=0, weights=depth[index]
        )
        donor_groups.append(index)
    weighted = residual * depth[:, None]
    total_depth = float(depth.sum())

    def statistic(current_codes: np.ndarray) -> float:
        numerator = np.zeros((len(cell_levels), values.shape[1]))
        np.add.at(numerator, current_codes, weighted)
        denominator = np.bincount(
            current_codes, weights=depth, minlength=len(cell_levels)
        )
        means = np.divide(
            numerator,
            denominator[:, None],
            out=np.zeros_like(numerator),
            where=denominator[:, None] > 0,
        )
        return float(
            np.sum((denominator / total_depth) * np.sum(means**2, axis=1))
        )

    observed = statistic(codes)
    greater_equal = 0
    for _ in range(permutations):
        permuted = codes.copy()
        for index in donor_groups:
            permuted[index] = rng.permutation(codes[index])
        greater_equal += statistic(permuted) >= observed
    return observed, (1 + greater_equal) / (permutations + 1)


def usage_effect(
    pivot: pd.DataFrame,
    donor: np.ndarray,
    cell_type: np.ndarray,
) -> tuple[float, float]:
    counts = pivot.to_numpy(dtype=float)
    usage = counts / counts.sum(axis=1, keepdims=True)
    levels = sorted(set(cell_type))
    equal_means = {
        label: usage[cell_type == label].mean(axis=0) for label in levels
    }
    pooled_means = {
        label: counts[cell_type == label].sum(axis=0)
        / counts[cell_type == label].sum()
        for label in levels
    }
    equal_tv = pooled_tv = 0.0
    for left_index, left in enumerate(levels):
        for right in levels[left_index + 1 :]:
            equal_tv = max(
                equal_tv, 0.5 * np.abs(equal_means[left] - equal_means[right]).sum()
            )
            pooled_tv = max(
                pooled_tv,
                0.5 * np.abs(pooled_means[left] - pooled_means[right]).sum(),
            )
    return float(equal_tv), float(pooled_tv)


def test_usage_gene(
    sub: pd.DataFrame,
    feature_col: str,
    permutations: int,
    rng: np.random.RandomState,
) -> dict:
    pivot = sub.pivot_table(
        index=["donor", "cell_type"],
        columns=feature_col,
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    meta = pivot.index.to_frame(index=False)
    donor = meta.donor.to_numpy()
    cell_type = meta.cell_type.to_numpy()
    depth = pivot.sum(axis=1).to_numpy(dtype=float)
    values = pivot.to_numpy(dtype=float) / depth[:, None]
    statistic, pvalue = permutation_stat(
        values, depth, donor, cell_type, permutations, rng
    )
    effect, pooled_effect = usage_effect(pivot, donor, cell_type)
    donor_counts = meta.groupby("cell_type").donor.nunique()
    return {
        "statistic": statistic,
        "pval": pvalue,
        "effect_equal_donor": effect,
        "effect_pooled_counts": pooled_effect,
        "n_features": int(pivot.shape[1]),
        "n_celltypes": int(meta.cell_type.nunique()),
        "n_donors": int(meta.donor.nunique()),
        "min_donors_per_celltype": int(donor_counts.min()),
        "max_donors_per_celltype": int(donor_counts.max()),
    }


def test_ase_gene(
    sub: pd.DataFrame,
    permutations: int,
    rng: np.random.RandomState,
) -> dict:
    donor = sub.donor.to_numpy()
    cell_type = sub.cell_type.to_numpy()
    depth = sub["count"].to_numpy(dtype=float)
    fraction = (sub.hapA / sub["count"]).to_numpy(dtype=float)[:, None]
    statistic, pvalue = permutation_stat(
        fraction, depth, donor, cell_type, permutations, rng
    )
    levels = sorted(set(cell_type))
    equal = {
        label: float(fraction[cell_type == label].mean()) for label in levels
    }
    pooled = {
        label: float(
            sub.loc[sub.cell_type.eq(label), "hapA"].sum()
            / sub.loc[sub.cell_type.eq(label), "count"].sum()
        )
        for label in levels
    }
    equal_hi = max(equal, key=equal.get)
    equal_lo = min(equal, key=equal.get)
    donor_counts = sub.groupby("cell_type").donor.nunique()
    return {
        "statistic": statistic,
        "pval": pvalue,
        "effect_equal_donor": equal[equal_hi] - equal[equal_lo],
        "effect_pooled_counts": max(pooled.values()) - min(pooled.values()),
        "hi_celltype": equal_hi,
        "lo_celltype": equal_lo,
        "n_features": 2,
        "n_celltypes": int(sub.cell_type.nunique()),
        "n_donors": int(sub.donor.nunique()),
        "min_donors_per_celltype": int(donor_counts.min()),
        "max_donors_per_celltype": int(donor_counts.max()),
    }


def main() -> None:
    args = parse_args()
    if args.permutations < 1:
        raise ValueError("--permutations must be positive")
    registry_path = Path(args.registry)
    pattern = args.input_glob or DEFAULT_INPUTS[args.analysis]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    registry = load_registry(registry_path, list(args.studies))
    raw, files = read_inputs(pattern, args.analysis, list(args.studies))
    bound, unit_audit = bind_units(
        raw, registry, args.allow_missing_registry_runs
    )
    if args.analysis == "ase":
        data, gene_col, feature_col, _ = prepare_ase(bound, args)
    else:
        data, gene_col, feature_col, _ = prepare_usage(bound, args.analysis, args)
    genes = eligible_genes(data, gene_col, feature_col, args.min_donors)
    if args.max_genes is not None:
        genes = genes[: args.max_genes]
    grouped = {str(key): frame for key, frame in data.groupby(gene_col)}
    rows = []
    for gene in genes:
        rng = stable_rng(args.seed, f"inner:{args.analysis}:{gene}")
        if args.analysis == "ase":
            result = test_ase_gene(grouped[gene], args.permutations, rng)
        else:
            result = test_usage_gene(
                grouped[gene], str(feature_col), args.permutations, rng
            )
        result["gene"] = gene
        rows.append(result)
    result = pd.DataFrame(rows)
    if result.empty:
        raise ValueError("no eligible genes after filtering")
    result["qval"] = multipletests(result.pval, method="fdr_bh")[1]
    result["sig"] = (
        result.qval.lt(0.05)
        & result.effect_equal_donor.ge(args.effect_threshold)
    )
    result = result.sort_values(["pval", "gene"]).reset_index(drop=True)
    result.to_csv(output, sep="\t", index=False)
    manifest = {
        "schema_version": "scTHREAD.biological_unit_celltype_permutation.v1",
        "analysis": args.analysis,
        "negative_control": args.negative_control,
        "studies": list(args.studies),
        "shared_celltypes_only": not args.include_private_celltypes,
        "shared_celltypes": SHARED,
        "aggregation_unit": "study_id::donor_or_source_id",
        "technical_runs_and_repeated_states": "summed_before_testing",
        "permutations": args.permutations,
        "seed": args.seed,
        "min_depth": args.min_depth,
        "min_donors": args.min_donors,
        "effect_threshold": args.effect_threshold,
        "unit_audit": unit_audit,
        "input_glob": pattern,
        "input_file_count": len(files),
        "registry_path": str(registry_path.resolve()),
        "registry_sha256": sha256(registry_path),
        "genes_tested": int(len(result)),
        "significant_genes": int(result.sig.sum()),
        "raw_p_lt_0_05_fraction": float(result.pval.lt(0.05).mean()),
        "output_path": str(output.resolve()),
        "output_sha256": sha256(output),
    }
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
