#!/usr/bin/env python3
"""Rank real E6.5 cell-type DTU genes for readable isoform-expression UMAPs."""

from __future__ import annotations

import os

import argparse
import json
from pathlib import Path

import duckdb
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
ASE_ROOT = Path(
    os.environ.get("SCTHREAD_ISOFORM_ROOT", "/gpfs/home/fuzc/project/ASE/cc/config/ISOform") + "/"
    "06_isoform_identification/dtu_analysis"
)
ST3 = ASE_ROOT / "tables/ST3_all_significant_DTU_isoforms.tsv"
USAGE = ASE_ROOT / "results/pseudotime_switching/E65_isoform_usage.tsv"
MOUSE_GENES = PROJECT.parent / "web/data/mouse_gene_intervals_grcm39.tsv"
OUTPUT = PROJECT / "tables/Fig3_mouse_scONT_DTU_UMAP_candidates.tsv"
AUDIT = PROJECT / "tables/Fig3_mouse_scONT_DTU_UMAP_candidates.json"


def rank_candidates(coverage: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    ranked_isoforms = coverage.sort_values(
        ["Gene_ID", "expressing_cells", "molecules"],
        ascending=[True, False, False],
    )
    top2 = ranked_isoforms.groupby("Gene_ID", as_index=False).head(2)
    ranked = (
        top2.groupby("Gene_ID")
        .agg(
            n_visible_isoforms=("Isoform_ID", "size"),
            top2_expressing_cells=("expressing_cells", "sum"),
            top2_molecules=("molecules", "sum"),
            max_abs_prop_diff=("abs_prop_diff", "max"),
            min_fdr=("FDR", "min"),
        )
        .reset_index()
    )
    ranked = ranked.loc[ranked["n_visible_isoforms"] >= 2].merge(
        genes, on="Gene_ID", how="left"
    )
    ranked = ranked.sort_values(
        ["top2_expressing_cells", "max_abs_prop_diff"],
        ascending=[False, False],
    )
    ranked["top_isoforms"] = ranked["Gene_ID"].map(
        top2.groupby("Gene_ID")["Isoform_ID"].apply(lambda x: ";".join(x))
    )
    return ranked


def self_test() -> None:
    coverage = pd.DataFrame(
        [
            ["g1", "i1", 0.01, 0.4, 10, 12.0],
            ["g1", "i2", 0.02, 0.3, 8, 9.0],
            ["g2", "i3", 0.01, 0.5, 30, 31.0],
        ],
        columns=[
            "Gene_ID",
            "Isoform_ID",
            "FDR",
            "abs_prop_diff",
            "expressing_cells",
            "molecules",
        ],
    )
    genes = pd.DataFrame({"Gene_ID": ["g1", "g2"], "gene_name": ["A", "B"]})
    ranked = rank_candidates(coverage, genes)
    assert ranked["Gene_ID"].tolist() == ["g1"]
    assert ranked.iloc[0]["top_isoforms"] == "i1;i2"
    print("SELF-TEST PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return

    significant = pd.read_csv(ST3, sep="\t")
    genes = (
        pd.read_csv(MOUSE_GENES, sep="\t", usecols=["gid", "gname"])
        .drop_duplicates("gid")
        .rename(columns={"gid": "Gene_ID", "gname": "gene_name"})
    )
    con = duckdb.connect(database=":memory:")
    con.register("significant", significant)
    usage = USAGE.as_posix().replace("'", "''")
    coverage = con.execute(
        f"""
        SELECT
          s.Gene_ID,
          s.Isoform_ID,
          min(s.FDR)::DOUBLE AS FDR,
          max(abs(s.Prop_Diff))::DOUBLE AS abs_prop_diff,
          count(DISTINCT u.cell_id)::INTEGER AS expressing_cells,
          sum(CAST(u.count AS DOUBLE))::DOUBLE AS molecules
        FROM significant s
        JOIN read_csv(
          '{usage}',
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
        ) u ON u.gene_id=s.Isoform_ID
        WHERE u.count > 0
        GROUP BY s.Gene_ID,s.Isoform_ID
        """
    ).fetchdf()
    ranked = rank_candidates(coverage, genes)
    ranked.to_csv(OUTPUT, sep="\t", index=False)
    audit = {
        "purpose": "coverage-aware selection of a real E6.5 cell-type DTU web-UMAP example",
        "selection_rule": "at least two significant isoforms; rank by summed non-zero cells of the two most-covered isoforms, then maximum absolute usage difference",
        "candidate_genes": int(len(ranked)),
        "top10": ranked.head(10).to_dict(orient="records"),
        "guardrail": "UMAP coverage ranks visualization suitability; DRIMSeq FDR and usage difference remain the inferential evidence",
    }
    AUDIT.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(ranked.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
