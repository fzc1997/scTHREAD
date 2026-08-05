#!/usr/bin/env python3
"""Build a traceable OWN_ASE_scONT per-cell UMAP overlay for the web portal.

The output contains real UMAP coordinates and long-read isoform counts for two
curated DTU examples (Tra2a and Tnrc6c). It does not infer DTU from the UMAP:
the UMAP localizes expression, while the existing DRIMSeq tables provide the
cell-type/stage-specific statistical evidence.
"""

from __future__ import annotations

import os

import argparse
import hashlib
import json
from pathlib import Path

import duckdb
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
TABLES = PROJECT / "tables"
WEB_DATA = PROJECT.parent / "web" / "data"
ASE_ROOT = Path(
    os.environ.get("SCTHREAD_ISOFORM_ROOT", "/gpfs/home/fuzc/project/ASE/cc/config/ISOform") + "/"
    "06_isoform_identification/dtu_analysis"
)

UMAP_TSV = TABLES / "OWN_ASE_scONT_all_cells_umap.tsv"
ANNOTATION_TSV = ASE_ROOT / "results/pseudotime_switching/slingshot_pseudotime.tsv"
USAGE_STAGE_TSV = {
    "E65": ASE_ROOT / "results/pseudotime_switching/E65_isoform_usage.tsv",
    "E75": ASE_ROOT / "results/pseudotime_switching/E75_isoform_usage.tsv",
    "E85": ASE_ROOT / "results/pseudotime_switching/E85_isoform_usage.tsv",
}
CELLTYPE_DTU = ASE_ROOT / "results/nature_cases_celltype.tsv"
DEVELOPMENTAL_DTU = ASE_ROOT / "results/nature_cases_developmental.tsv"
ALL_CELLTYPE_DTU_GENES = ASE_ROOT / "tables/ST2_all_significant_DTU_genes.tsv"
ALL_CELLTYPE_DTU_ISOFORMS = ASE_ROOT / "tables/ST3_all_significant_DTU_isoforms.tsv"

POINTS_OUT = WEB_DATA / "mouse_scont_umap_points.parquet"
ISOFORMS_OUT = WEB_DATA / "mouse_scont_umap_isoform_expression.parquet"
AUDIT_OUT = WEB_DATA / "mouse_scont_umap_audit.json"
DTU_OUT = TABLES / "Fig3_mouse_scONT_portal_DTU_examples.tsv"

TARGETS = {
    "ENSMUSG00000092341": "Malat1",
    "ENSMUSG00000029817": "Tra2a",
    "ENSMUSG00000025571": "Tnrc6c",
}
TARGET_ISOFORMS = {
    "ENSMUSG00000092341": [
        "ENSMUST00000245150",
        "ENSMUST00000172812",
    ],
    "ENSMUSG00000029817": [
        "ENSMUST00000031841",
        "ENSMUST00000204013",
    ],
    "ENSMUSG00000025571": [
        "ENSMUST00000026658",
        "ENSMUST00000138299",
    ],
}
STAGE_SUFFIX = {"_1": "E6.5", "_2": "E7.5", "_3": "E8.5"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_points(umap: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    required_umap = {"cell_barcode", "UMAP_1", "UMAP_2"}
    required_annotation = {"cell_id", "celltype"}
    if not required_umap.issubset(umap.columns):
        raise ValueError(f"UMAP columns missing: {required_umap - set(umap.columns)}")
    if not required_annotation.issubset(annotation.columns):
        raise ValueError(
            f"Annotation columns missing: {required_annotation - set(annotation.columns)}"
        )
    points = umap.merge(
        annotation[["cell_id", "celltype"]].drop_duplicates("cell_id"),
        left_on="cell_barcode",
        right_on="cell_id",
        how="left",
        validate="one_to_one",
    )
    if points["celltype"].isna().any():
        raise ValueError(
            f"{int(points['celltype'].isna().sum())} UMAP cells lack annotations"
        )
    points = points.drop(columns=["cell_id"])
    suffix = points["cell_barcode"].str.extract(r"(_[123])$", expand=False)
    points["stage"] = suffix.map(STAGE_SUFFIX)
    if points["stage"].isna().any():
        raise ValueError("Unexpected cell-id stage suffix")
    points = points.rename(
        columns={
            "cell_barcode": "cell_id",
            "UMAP_1": "umap1",
            "UMAP_2": "umap2",
        }
    )
    points["gse"] = points["stage"]
    points["run"] = "OWN_ASE_scONT_" + points["stage"].str.replace(".", "", regex=False)
    points["source_study"] = "OWN_ASE_scONT"
    return points[
        [
            "cell_id",
            "gse",
            "run",
            "celltype",
            "umap1",
            "umap2",
            "stage",
            "source_study",
        ]
    ].rename(columns={"celltype": "cell_type"})


def self_test() -> None:
    umap = pd.DataFrame(
        {
            "cell_barcode": ["AAA-1_1", "BBB-1_3"],
            "UMAP_1": [1.0, 2.0],
            "UMAP_2": [-1.0, -2.0],
        }
    )
    annotation = pd.DataFrame(
        {
            "cell_id": ["AAA-1_1", "BBB-1_3"],
            "celltype": ["Neuroectoderm", "Embryonic_Epiblast"],
        }
    )
    points = build_points(umap, annotation)
    assert points["stage"].tolist() == ["E6.5", "E8.5"]
    assert points["gse"].tolist() == ["E6.5", "E8.5"]
    assert points["source_study"].eq("OWN_ASE_scONT").all()
    print("SELF-TEST PASS")


def build() -> None:
    for path in (
        UMAP_TSV,
        ANNOTATION_TSV,
        CELLTYPE_DTU,
        DEVELOPMENTAL_DTU,
        ALL_CELLTYPE_DTU_GENES,
        ALL_CELLTYPE_DTU_ISOFORMS,
        *USAGE_STAGE_TSV.values(),
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    umap = pd.read_csv(UMAP_TSV, sep="\t")
    annotation = pd.read_csv(
        ANNOTATION_TSV,
        sep="\t",
        usecols=["cell_id", "celltype"],
    )
    points = build_points(umap, annotation)
    if len(points) != 25_621 or points["cell_id"].nunique() != 25_621:
        raise ValueError("Expected exactly 25,621 unique OWN_ASE_scONT cells")

    WEB_DATA.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(database=":memory:")
    con.register("points_df", points)
    con.execute(
        f"COPY points_df TO '{POINTS_OUT.as_posix()}' "
        "(FORMAT PARQUET, COMPRESSION ZSTD)"
    )

    celltype = pd.read_csv(CELLTYPE_DTU, sep="\t")
    developmental = pd.read_csv(DEVELOPMENTAL_DTU, sep="\t")
    transcript_to_gene = {
        transcript_id: gene_id
        for gene_id, transcript_ids in TARGET_ISOFORMS.items()
        for transcript_id in transcript_ids
    }
    if set(transcript_to_gene.values()) != set(TARGETS):
        raise ValueError("Both curated genes must have target ENSMUST isoforms")
    transcript_ids = ",".join(f"'{item}'" for item in transcript_to_gene)
    gene_case = " ".join(
        f"WHEN gene_id='{transcript}' THEN '{gene_id}'"
        for transcript, gene_id in transcript_to_gene.items()
    )
    stage_queries = []
    for stage, path in USAGE_STAGE_TSV.items():
        suffix = {"E65": "_1", "E75": "_2", "E85": "_3"}[stage]
        usage_path = path.as_posix().replace("'", "''")
        stage_queries.append(
            f"""
            SELECT
              cell_id || '{suffix}' AS cell_id,
              CASE {gene_case} END AS gene_id,
              gene_id AS transcript_id,
              CAST(count AS DOUBLE) AS count,
              '{stage}' AS stage
            FROM read_csv(
              '{usage_path}',
              delim='\\t',
              header=true,
              columns={{
                'cell_id':'VARCHAR',
                'gene_id':'VARCHAR',
                'isoform_id':'VARCHAR',
                'count':'DOUBLE',
                'gene_total':'DOUBLE',
                'usage':'DOUBLE'
              }}
            )
            WHERE gene_id IN ({transcript_ids}) AND count > 0
            """
        )
    con.execute("CREATE TABLE target_usage AS " + " UNION ALL ".join(stage_queries))
    missing_cells = con.execute(
        """
        SELECT count(*) FROM target_usage u
        LEFT JOIN points_df p USING(cell_id)
        WHERE p.cell_id IS NULL
        """
    ).fetchone()[0]
    target_records = con.execute("SELECT count(*) FROM target_usage").fetchone()[0]
    missing_breakdown = con.execute(
        """
        SELECT u.stage, count(*)::INTEGER AS records,
               count(DISTINCT u.cell_id)::INTEGER AS cells
        FROM target_usage u LEFT JOIN points_df p USING(cell_id)
        WHERE p.cell_id IS NULL
        GROUP BY u.stage ORDER BY u.stage
        """
    ).fetchdf()
    if missing_cells / max(target_records, 1) > 0.35:
        raise ValueError(
            f"{missing_cells}/{target_records} target-expression rows lack UMAP cells"
        )
    con.execute(
        """
        CREATE TABLE matched_target_usage AS
        SELECT u.cell_id,u.gene_id,u.transcript_id,u.count,u.stage
        FROM target_usage u JOIN points_df p USING(cell_id)
        """
    )
    con.execute(
        f"COPY (SELECT cell_id,gene_id,transcript_id,count FROM matched_target_usage) "
        f"TO '{ISOFORMS_OUT.as_posix()}' (FORMAT PARQUET, COMPRESSION ZSTD)"
    )

    overlay_summary = con.execute(
        """
        SELECT gene_id, transcript_id, count(DISTINCT cell_id)::INTEGER AS cells,
               sum(count)::DOUBLE AS count_units
        FROM matched_target_usage
        GROUP BY gene_id, transcript_id
        ORDER BY gene_id, count_units DESC
        """
    ).fetchdf()
    if set(overlay_summary["gene_id"]) != set(TARGETS):
        raise ValueError("Both curated genes must have per-cell overlay records")

    celltype = celltype.loc[celltype["gene_id"].isin(TARGETS)].copy()
    celltype["contrast"] = "E6.5 Neuroectoderm vs Embryonic epiblast"
    celltype["evidence_type"] = "cell-type DTU"
    developmental = developmental.loc[developmental["gene_id"].isin(TARGETS)].copy()
    developmental["contrast"] = "E6.5/E7.5/E8.5 within lineage"
    developmental["evidence_type"] = "stage-associated DTU"
    common = [
        "gene_id",
        "gene_name",
        "isoform_id",
        "contrast",
        "evidence_type",
    ]
    all_isoform_dtu = pd.read_csv(ALL_CELLTYPE_DTU_ISOFORMS, sep="\t")
    malat1 = all_isoform_dtu.loc[
        all_isoform_dtu["Isoform_ID"].isin(TARGET_ISOFORMS["ENSMUSG00000092341"])
    ].rename(
        columns={
            "Gene_ID": "gene_id",
            "Isoform_ID": "isoform_id",
            "FDR": "feature_fdr",
            "Prop_Neuroectoderm": "prop_Neuroectoderm",
            "Prop_Epiblast": "prop_Epiblast",
            "Prop_Diff": "prop_diff",
        }
    )
    malat1["gene_name"] = "Malat1"
    malat1["contrast"] = "E6.5 Neuroectoderm vs Embryonic epiblast"
    malat1["evidence_type"] = "cell-type DTU"
    dtu_examples = pd.concat(
        [
            celltype[common + ["prop_Neuroectoderm", "prop_Epiblast", "prop_diff"]],
            malat1[
                common
                + [
                    "prop_Neuroectoderm",
                    "prop_Epiblast",
                    "prop_diff",
                    "feature_fdr",
                ]
            ],
            developmental[
                common + ["prop_E65", "prop_E75", "prop_E85", "prop_diff_E65_E85"]
            ],
        ],
        ignore_index=True,
        sort=False,
    )
    dtu_examples.to_csv(DTU_OUT, sep="\t", index=False)

    audit = {
        "dataset": "OWN_ASE_scONT mouse gastrulation",
        "method": "project Seurat UMAP; long-read isoform counts overlaid without coordinate refitting",
        "expression_unit": "IsoQuant-derived long-read count units (fractional values retained)",
        "cells": int(len(points)),
        "cell_types": int(points["cell_type"].nunique()),
        "stages": points["stage"].value_counts().sort_index().to_dict(),
        "genes": TARGETS,
        "overlay_records": int(
            con.execute("SELECT count(*) FROM matched_target_usage").fetchone()[0]
        ),
        "overlay_records_without_umap": int(missing_cells),
        "overlay_records_without_umap_by_stage": missing_breakdown.to_dict(
            orient="records"
        ),
        "overlay_summary": overlay_summary.to_dict(orient="records"),
        "guardrails": [
            "UMAP localizes isoform expression and does not itself test DTU",
            "Malat1 cell-type DTU statistics come from the E6.5 DRIMSeq comparison",
            "Tra2a cell-type DTU statistics come from the E6.5 DRIMSeq comparison",
            "Tnrc6c stage-associated statistics use cell-bootstrap pseudo-replicates, not independent embryos",
        ],
        "inputs": {
            str(path): {"sha256": sha256(path)}
            for path in (
                UMAP_TSV,
                ANNOTATION_TSV,
                CELLTYPE_DTU,
                DEVELOPMENTAL_DTU,
                ALL_CELLTYPE_DTU_GENES,
                ALL_CELLTYPE_DTU_ISOFORMS,
                *USAGE_STAGE_TSV.values(),
            )
        },
        "outputs": [str(POINTS_OUT), str(ISOFORMS_OUT), str(DTU_OUT)],
    }
    AUDIT_OUT.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        build()


if __name__ == "__main__":
    main()
