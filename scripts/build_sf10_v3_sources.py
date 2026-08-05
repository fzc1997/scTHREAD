#!/usr/bin/env python3
"""Build SF10 v3 source tables from the current portal mouse embedding.

SF10 previously asserted a 25,621-cell / 24-cell-type portal cell map. The portal
was rebuilt on 30 Jul 2026 onto the full six-run embedding (Forward and
Reciprocal crosses, 55,729 cells, 36 cell types), so the figure could no longer
be reproduced. These tables are derived directly from the files the portal
serves, so the figure and the site cannot drift apart again.

Note the two denominators this dataset carries, which must never be mixed:
  * portal cell map  - 55,729 cells / 36 cell types / 6 runs (this script)
  * ANCHOR DTU analysis - 25,621 cells / 24 cell types, forward cross only
"""
from __future__ import annotations

import os

import hashlib
import json
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WEB = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD") + "/web/data")
POINTS = WEB / "mouse_scont_umap_points.parquet"
ISOFORMS = WEB / "mouse_scont_umap_isoform_expression.parquet"
AUDIT = WEB / "mouse_scont_umap_audit.json"
OUT = ROOT / "tables/sf10_v3"

MALAT1 = "ENSMUSG00000092341"
ISO1 = "ENSMUST00000245150"
ISO2 = "ENSMUST00000172812"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("PRAGMA threads=4")
    pts, iso = str(POINTS), str(ISOFORMS)

    composition = con.execute(
        f"""SELECT "cross", stage, run, count(*) AS cells,
                   count(DISTINCT cell_type) AS cell_types
            FROM read_parquet('{pts}') GROUP BY 1,2,3 ORDER BY 1,2"""
    ).df()
    composition.to_csv(OUT / "embedding_composition.tsv", sep="\t", index=False)

    celltypes = con.execute(
        f"""SELECT cell_type,
                   count(*) AS cells,
                   sum(CASE WHEN stage='E6.5' THEN 1 ELSE 0 END) AS E6_5,
                   sum(CASE WHEN stage='E7.5' THEN 1 ELSE 0 END) AS E7_5,
                   sum(CASE WHEN stage='E8.5' THEN 1 ELSE 0 END) AS E8_5,
                   sum(CASE WHEN "cross"='Forward' THEN 1 ELSE 0 END) AS forward_cells,
                   sum(CASE WHEN "cross"='Reciprocal' THEN 1 ELSE 0 END) AS reciprocal_cells
            FROM read_parquet('{pts}') GROUP BY 1 ORDER BY cells DESC"""
    ).df()
    celltypes.to_csv(OUT / "celltype_composition.tsv", sep="\t", index=False)

    # per-cell coordinates with the two displayed Malat1 isoforms
    umap = con.execute(
        f"""WITH sel AS (
              SELECT cell_id,
                     sum(CASE WHEN transcript_id='{ISO1}' THEN count ELSE 0 END) AS iso1,
                     sum(CASE WHEN transcript_id='{ISO2}' THEN count ELSE 0 END) AS iso2,
                     sum(count) AS malat1_total
              FROM read_parquet('{iso}')
              WHERE gene_id='{MALAT1}'
              GROUP BY 1)
            SELECT p.cell_id, p.run, p."cross", p.stage, p.cell_type,
                   p.umap1, p.umap2,
                   coalesce(s.iso1, 0) AS iso1_count,
                   coalesce(s.iso2, 0) AS iso2_count,
                   coalesce(s.malat1_total, 0) AS malat1_total
            FROM read_parquet('{pts}') p LEFT JOIN sel s USING (cell_id)"""
    ).df()
    umap.to_parquet(OUT / "malat1_umap.parquet", index=False)

    # within-gene isoform usage per cell type, on molecule counts
    usage = con.execute(
        f"""WITH j AS (
              SELECT p.cell_type, e.transcript_id, e.count
              FROM read_parquet('{iso}') e
              JOIN read_parquet('{pts}') p USING (cell_id)
              WHERE e.gene_id='{MALAT1}')
            SELECT cell_type, transcript_id, sum(count) AS molecules
            FROM j GROUP BY 1,2"""
    ).df()
    totals = usage.groupby("cell_type")["molecules"].sum().rename("gene_molecules")
    usage = usage.join(totals, on="cell_type")
    usage["usage_fraction"] = usage["molecules"] / usage["gene_molecules"]
    usage.sort_values(["cell_type", "molecules"], ascending=[True, False]).to_csv(
        OUT / "malat1_usage_by_celltype.tsv", sep="\t", index=False
    )

    # headline contrast for the two displayed isoforms, eligible cell types only
    eligible = usage[usage["gene_molecules"] >= 100]
    pair = eligible[eligible["transcript_id"].isin({ISO1, ISO2})].pivot(
        index="cell_type", columns="transcript_id", values="usage_fraction"
    ).fillna(0.0)
    pair["difference"] = pair[ISO1] - pair[ISO2]
    pair = pair.sort_values("difference")
    pair.to_csv(OUT / "malat1_iso_pair_by_celltype.tsv", sep="\t")

    audit = json.loads(AUDIT.read_text())
    summary = {
        "portal_embedding": {
            "cells": int(len(umap)),
            "runs": int(umap["run"].nunique()),
            "cell_types": int(umap["cell_type"].nunique()),
            "crosses": sorted(umap["cross"].unique()),
            "stages": sorted(umap["stage"].unique()),
            "cells_by_cross_stage": {
                f"{row.cross}:{row.stage}": int(row.cells)
                for row in composition.itertuples(index=False)
            },
        },
        "audit_json_agrees": {
            "cells": int(audit["cells"]) == int(len(umap)),
            "runs": int(audit["runs"]) == int(umap["run"].nunique()),
            "cell_types": int(audit["cell_types"]) == int(umap["cell_type"].nunique()),
        },
        "malat1": {
            "gene_id": MALAT1,
            "cells_with_signal": int((umap["malat1_total"] > 0).sum()),
            "displayed_isoforms": {
                ISO1: {
                    "cells": int((umap["iso1_count"] > 0).sum()),
                    "molecules": float(umap["iso1_count"].sum()),
                },
                ISO2: {
                    "cells": int((umap["iso2_count"] > 0).sum()),
                    "molecules": float(umap["iso2_count"].sum()),
                },
            },
            "eligible_cell_types_min100_molecules": int(len(pair)),
            "max_usage_difference": float(pair["difference"].abs().max()),
            "most_iso1_biased_cell_type": str(pair["difference"].idxmax()),
            "most_iso2_biased_cell_type": str(pair["difference"].idxmin()),
        },
        "inputs": {str(p): sha256(p) for p in (POINTS, ISOFORMS, AUDIT)},
    }
    (OUT / "sf10_v3_manifest.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
