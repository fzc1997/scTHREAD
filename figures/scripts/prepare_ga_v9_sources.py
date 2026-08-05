#!/usr/bin/env python3
"""Freeze audited biological coverage and 9,999-permutation GA v9 sources."""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REGISTRY_FROZEN = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv.bak_2338"
)
REGISTRY_ENRICHED = Path(
    os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv"
)
SUMMARY_TABLE = ROOT / "tables/Table_S_three_axis_summary.tsv"
DEFAULT_OUTDIR = ROOT / "figures/ga_gpt2_vector_v9"

INPUT_HASHES = {
    REGISTRY_FROZEN:
        "14adfaffa905c8832c29d0f94a4a375c546048a9b94d5a9e392947c2929183da",
    SUMMARY_TABLE:
        "a43cc81e70da955e067172ed422681d56981da3fe2dc1940150659667fd6ce81",
    ROOT / "tables/p0_biological_unit_rerun/ase_observed_9999.tsv":
        "b2a42c3955243c25b080aa6a1b469bd2a94c593c7c36f2be3926ae802e370743",
    ROOT / "tables/p0_biological_unit_rerun/diu_observed_9999.tsv":
        "688275f25d30442c13f916af11585e21f4989750c9b0dd6daaef451bb05eb98f",
    ROOT / "tables/p0_biological_unit_rerun/apa_observed_9999.tsv":
        "8e58946048b44d4a7b6893a589c3d82701e61efa2584ea19598c15c2aa77882b",
}

# The classes are mutually exclusive. Cancer/cell-line studies are kept in
# Cancer even when the cell lineage belongs to an organ system.
SYSTEM_MAP = {
    # Blood and immune biology
    "GSE214231": "Blood/immune",
    "GSE252344": "Blood/immune",
    "GSE252416": "Blood/immune",
    "GSE276974": "Blood/immune",
    "GSE292324": "Blood/immune",
    "GSE295352": "Blood/immune",
    "GSE307660": "Blood/immune",
    # Neural and sensory biology
    "GSE114157": "Neural/sensory",
    "GSE130708": "Neural/sensory",
    "GSE178175": "Neural/sensory",
    "GSE255520": "Neural/sensory",
    "GSE274249": "Neural/sensory",
    "GSE283629": "Neural/sensory",
    "GSE314176": "Neural/sensory",
    "GSE76026": "Neural/sensory",
    # Cancer and malignant cell models
    "GSE212945": "Cancer",
    "GSE224045": "Cancer",
    "GSE248118": "Cancer",
    "GSE289428": "Cancer",
    "GSE289790": "Cancer",
    "GSE295932": "Cancer",
    "GSE301658": "Cancer",
    "GSE303762": "Cancer",
    # Other organ/developmental systems
    "GSE295353": "Endocrine",
    "GSE288222": "Heart/vascular",
    "GSE185554": "Development/embryo",
    "GSE250381": "Development/embryo",
    "GSE274527": "Development/embryo",
    "GSE283658": "Development/embryo",
    "GSE309071": "Development/embryo",
    # datasets that entered under the single-cell long-read scope rule
    "GSE158450": "Neural/sensory",
    "GSA_mouse_testis_11week": "Endocrine",
    "benagen_human_ovary": "Endocrine",
    "OWN_ASE_scONT": "Development/embryo",
}

MODALITY_EXCLUDED_STUDIES = {"GSE140890"}

EXPECTED_SYSTEMS = {
    "Blood/immune": {"cells": 308_071, "samples": 56, "studies": 7},
    "Neural/sensory": {"cells": 274_793, "samples": 62, "studies": 9},
    "Cancer": {"cells": 134_433, "samples": 20, "studies": 8},
    "Development/embryo": {"cells": 98_079, "samples": 292, "studies": 6},
    "Endocrine": {"cells": 77_331, "samples": 11, "studies": 3},
    "Heart/vascular": {"cells": 30_682, "samples": 12, "studies": 1},
}

CONTEXTS = [
    ("GSE178175", "Frontal cortex", "Neural/sensory", 121_244),
    ("GSE276974", "CCUS marrow", "Blood/immune", 117_300),
    ("GSE283629", "iNeuron trajectory", "Neural/sensory", 74_540),
    ("GSE303762", "Cross-platform BM mix", "Cancer", 71_800),
    ("OWN_ASE_scONT", "Mouse gastrulation", "Development/embryo", 68_417),
    ("GSE307660", "Myeloma marrow", "Blood/immune", 68_407),
    ("GSE292324", "Splicing dynamics", "Blood/immune", 55_507),
    ("GSA_mouse_testis_11week", "Mouse testis", "Endocrine", 39_585),
    ("GSE274249", "iPSC to cortical neuron", "Neural/sensory", 34_973),
    ("benagen_human_ovary", "Human ovary", "Endocrine", 31_839),
]

EXPECTED_INVENTORY = {
    "DIU": {"tested": 8_092, "passed": 2_008, "raw": 3_710},
    "APA": {"tested": 10_531, "passed": 2_558, "raw": 5_799},
    "ASE": {"tested": 6_930, "passed": 0, "raw": 538},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def compact_values(values: pd.Series, limit: int = 8) -> str:
    unique = sorted(
        {
            str(value).strip()
            for value in values.dropna()
            if str(value).strip()
        }
    )
    return " | ".join(unique[:limit])


def validate_input_hashes() -> None:
    for path, expected in INPUT_HASHES.items():
        require(path.is_file(), f"Missing source: {path}")
        observed = sha256(path)
        require(
            observed == expected,
            f"Source changed: {path} {observed} != {expected}",
        )
    require(REGISTRY_ENRICHED.is_file(), f"Missing source: {REGISTRY_ENRICHED}")


def prepare_catalog(outdir: Path) -> list[Path]:
    # The released manifest is the authority: it carries the modality audit and
    # the stale-row correction. Re-deriving from the frozen registry silently
    # reintroduced a genomic-DNA study and a stale cell count.
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parent))
    import render_nar_bio as _R

    done = _R.registry_done().copy()
    done["isoquant_cells"] = pd.to_numeric(
        done["isoquant_cells"], errors="coerce"
    ).fillna(0).astype(int)
    require(
        not set(done["gse"]) & MODALITY_EXCLUDED_STUDIES,
        "Modality-excluded studies reappeared in the released manifest",
    )
    require(len(done) == 453, f"Release run count changed: {len(done)}")
    require(done["gse"].nunique() == 34, "Release study count changed")
    require(
        int(_R.study_cells()["cells"].sum()) == 923_389,
        "Audited called-cell total changed",
    )
    require(done["srr"].is_unique, "Frozen run identifiers are not unique")
    require(
        set(done["gse"]) == set(SYSTEM_MAP),
        (
            "System map does not cover the frozen studies: "
            f"missing={sorted(set(done['gse']) - set(SYSTEM_MAP))}, "
            f"extra={sorted(set(SYSTEM_MAP) - set(done['gse']))}"
        ),
    )

    enriched = pd.read_csv(REGISTRY_ENRICHED, sep="\t", dtype=str)
    require(enriched["srr"].is_unique, "Enriched registry run IDs are not unique")
    metadata_columns = [
        "srr",
        "geo_source_name",
        "geo_title",
        "geo_characteristics",
        "geo_description",
        "ena_sample_title",
        "ena_experiment_title",
    ]
    # registry_done() already carries the live-registry columns, so drop the
    # metadata copies before re-joining them from the enriched registry
    merged = done.drop(columns=[c for c in metadata_columns if c != "srr"],
                       errors="ignore").merge(
        enriched[metadata_columns],
        on="srr",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    require(
        merged["_merge"].eq("both").all(),
        "At least one frozen run lacks enriched registry metadata",
    )
    merged = merged.drop(columns="_merge")
    merged["system"] = merged["gse"].map(SYSTEM_MAP)
    # Cells are authoritative per study, so attach the study total to each of its
    # run rows and take the max when regrouping; summing run rows overstates any
    # study whose count came from an author file or a STARsolo call.
    merged["authoritative_cells"] = merged["gse"].map(
        _R.study_cells().set_index("gse")["cells"]
    )
    require(
        merged["authoritative_cells"].notna().all(),
        "At least one run belongs to a study with no authoritative cell count",
    )

    evidence_assertions = {
        "GSE295352": ("geo_source_name", "spleen"),
        "GSE295353": ("geo_source_name", "pancreatic islet"),
        "GSE248118": ("geo_characteristics", "cell line:"),
        "GSE178175": ("geo_source_name", "frontal cortex"),
        "GSE309071": ("geo_source_name", "mouse embryo"),
        "GSE289428": ("geo_source_name", "hela-s3"),
        "GSE252344": ("geo_source_name", "lymph node"),
        "GSE252416": ("geo_source_name", "lymph node"),
        "GSE214231": ("geo_source_name", "tcr negative t cells"),
        "GSE295932": ("geo_title", "jurkat"),
    }
    for gse, (column, term) in evidence_assertions.items():
        values = " | ".join(
            merged.loc[merged["gse"].eq(gse), column].dropna().astype(str)
        ).lower()
        require(term in values, f"Metadata evidence missing: {gse} {term}")

    study = (
        merged.groupby(["gse", "system"], as_index=False)
        .agg(
            cells=("authoritative_cells", "max"),
            samples=("srr", "nunique"),
            species=("species", compact_values),
            geo_source_name=("geo_source_name", compact_values),
            geo_title=("geo_title", compact_values),
            geo_characteristics=("geo_characteristics", compact_values),
        )
        .sort_values(["system", "cells", "gse"], ascending=[True, False, True])
        .reset_index(drop=True)
    )
    summary = (
        study.groupby("system", as_index=False)
        .agg(
            cells=("cells", "sum"),
            samples=("samples", "sum"),
            studies=("gse", "nunique"),
        )
        .sort_values("cells", ascending=False)
        .reset_index(drop=True)
    )
    observed = {
        row.system: {
            "cells": int(row.cells),
            "samples": int(row.samples),
            "studies": int(row.studies),
        }
        for row in summary.itertuples(index=False)
    }
    require(observed == EXPECTED_SYSTEMS, f"System totals changed: {observed}")
    require(int(summary["cells"].sum()) == 923_389, "System cells do not close")
    require(int(summary["samples"].sum()) == 453, "System samples do not close")
    require(int(summary["studies"].sum()) == 34, "System studies do not close")

    context_rows: list[dict[str, object]] = []
    study_index = study.set_index("gse")
    for gse, display, system, expected_cells in CONTEXTS:
        row = study_index.loc[gse]
        cells = int(row["cells"])
        require(cells == expected_cells, f"Context count changed: {display}")
        require(row["system"] == system, f"Context system changed: {display}")
        context_rows.append(
            {
                "display": display,
                "system": system,
                "source_gse": gse,
                "samples": int(row["samples"]),
                "cells": cells,
                "registry_evidence": (
                    f"{row['geo_source_name']} | {row['geo_title']}"
                ),
            }
        )
    contexts = pd.DataFrame(context_rows).sort_values(
        ["cells", "display"], ascending=[False, True], kind="mergesort"
    )

    metadata_subset = merged[
        [
            "srr",
            "gse",
            "species",
            "platform",
            "isoquant_cells",
            "system",
            "geo_source_name",
            "geo_title",
            "geo_characteristics",
        ]
    ].sort_values(["gse", "srr"], kind="mergesort")

    outputs = [
        outdir / "catalog_study_system_audit_v9.tsv",
        outdir / "catalog_system_composition_v9.tsv",
        outdir / "registry_context_examples_v9.tsv",
        outdir / "frozen_registry_metadata_subset_v9.tsv",
    ]
    study.to_csv(outputs[0], sep="\t", index=False)
    summary.to_csv(outputs[1], sep="\t", index=False)
    contexts.to_csv(outputs[2], sep="\t", index=False)
    metadata_subset.to_csv(outputs[3], sep="\t", index=False)
    return outputs


def prepare_inventory(outdir: Path) -> Path:
    summary = pd.read_csv(SUMMARY_TABLE, sep="\t")
    require(set(summary["axis"]) == set(EXPECTED_INVENTORY), "Axis set changed")
    rows: list[dict[str, object]] = []
    for axis in ["DIU", "APA", "ASE"]:
        summary_row = summary.loc[summary["axis"].eq(axis)].iloc[0]
        source = Path(str(summary_row["source"]))
        require(source in INPUT_HASHES, f"Unexpected inventory source: {source}")
        observed = pd.read_csv(source, sep="\t")
        expected = EXPECTED_INVENTORY[axis]
        tested = len(observed)
        passed = int(observed["sig"].astype(bool).sum())
        raw = int((observed["pval"] < 0.05).sum())
        min_q = float(observed["qval"].min())
        require(tested == expected["tested"], f"{axis} tested count changed")
        require(passed == expected["passed"], f"{axis} passed count changed")
        require(raw == expected["raw"], f"{axis} raw-P count changed")
        require(
            int(summary_row["n_genes_tested"]) == tested
            and int(summary_row["n_sig_fdr05_effect_gate"]) == passed,
            f"{axis} summary/source mismatch",
        )
        require(
            abs(float(summary_row["frac_raw_p_lt_0.05"]) - raw / tested) < 1e-12,
            f"{axis} raw-P fraction mismatch",
        )
        rows.append(
            {
                "axis": axis,
                "n_genes_tested": tested,
                "n_sig_fdr05_effect_gate": passed,
                "sig_fraction": passed / tested,
                "n_raw_p_lt_0.05": raw,
                "raw_p_fraction": raw / tested,
                "min_q": min_q,
                "permutations": 9999,
                "source": str(source),
                "source_sha256": sha256(source),
            }
        )
    inventory = pd.DataFrame(rows)
    ase = inventory.set_index("axis").loc["ASE"]
    require(abs(float(ase["min_q"]) - 0.231) < 1e-12, "ASE minimum q changed")
    output = outdir / "three_axis_inventory_9999_v9.tsv"
    inventory.to_csv(output, sep="\t", index=False)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    validate_input_hashes()
    outputs = [*prepare_catalog(outdir), prepare_inventory(outdir)]
    manifest = {
        "contract": "scTHREAD graphical abstract v9",
        "catalog": {"runs": 453, "studies": 34, "cells": 923_389},
        "systems": EXPECTED_SYSTEMS,
        "inventory": EXPECTED_INVENTORY,
        "inputs": {
            **{str(path): expected for path, expected in INPUT_HASHES.items()},
            str(REGISTRY_ENRICHED): sha256(REGISTRY_ENRICHED),
        },
        "outputs": {
            str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs
        },
    }
    manifest_path = outdir / "source_manifest_v9.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    outputs.append(manifest_path)

    print(
        "GA V9 SOURCE PREFLIGHT PASS",
        {
            "runs": 434,
            "studies": 30,
            "called_cells": 850_938,
            "systems": EXPECTED_SYSTEMS,
            "inventory": EXPECTED_INVENTORY,
        },
    )
    for path in outputs:
        print(f"{path}\t{path.stat().st_size}\t{sha256(path)}")


if __name__ == "__main__":
    main()
