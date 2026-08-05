#!/usr/bin/env python3
"""Ask whether the precomputed cell-context annotations land on disease genes.

A catalog of candidate genes is only useful if it points somewhere. This tests
the simplest version of that: are genes flagged for cell-type-associated
isoform or poly(A)-site usage more likely than expected to carry a pathogenic
ClinVar variant?

The background is the set of genes that were actually **tested**, not the whole
genome. That matters: testability tracks expression and coverage, and so does
the chance of a gene having been sequenced enough to accumulate ClinVar
submissions, so a genome-wide background would manufacture enrichment out of
study effort alone.

Writes tables/clinvar_enrichment.{tsv,json}.
"""
from __future__ import annotations

import os

import gzip
import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
GTF = Path(os.environ.get("SCTHREAD_GENCODE_ROOT", "/gpfs/home/fuzc/project/ABE_transcriptomics_off_target/dep_files") + "/"
           "hg38/hg38_annotation/gencode.v41.annotation.gtf")
CLINVAR = Path(os.environ.get("SCTHREAD_CISREG_ROOT", "/gpfs/home/fuzc/project/ASE/papers/paper2_cis_regulation/data") + "/"
               "external/human_relevance/variant_summary.txt.gz")
AXES = {
    "isoform usage": ROOT / "tables/p0_biological_unit_rerun/diu_observed_9999.tsv",
    "poly(A)-site usage": ROOT / "tables/p0_biological_unit_rerun/apa_observed_9999.tsv",
    "allele-specific expression": ROOT / "tables/p0_biological_unit_rerun/ase_observed_9999.tsv",
}
EFFECT_GATE = 0.20
QVAL_GATE = 0.05


def ensembl_to_symbol() -> dict[str, str]:
    """Map versionless ENSG to gene symbol from the annotation build in use."""
    mapping: dict[str, str] = {}
    gene_id = re.compile(r'gene_id "([^."]+)')
    gene_name = re.compile(r'gene_name "([^"]+)')
    with GTF.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gid, name = gene_id.search(fields[8]), gene_name.search(fields[8])
            if gid and name:
                mapping[gid.group(1)] = name.group(1)
    return mapping


def pathogenic_genes() -> set[str]:
    """Gene symbols carrying at least one GRCh38 pathogenic ClinVar variant.

    "Conflicting" and "Pathogenic/Likely pathogenic, low penetrance" style
    records still count as pathogenic evidence only when the field starts with
    Pathogenic or Likely pathogenic; anything conflicting is dropped.
    """
    genes: set[str] = set()
    with gzip.open(CLINVAR, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        sig_i = idx["ClinicalSignificance"]
        sym_i = idx["GeneSymbol"]
        asm_i = idx["Assembly"]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sig_i, sym_i, asm_i):
                continue
            if fields[asm_i] != "GRCh38":
                continue
            sig = fields[sig_i]
            if "Conflicting" in sig:
                continue
            if not sig.startswith(("Pathogenic", "Likely pathogenic")):
                continue
            symbol = fields[sym_i]
            if symbol and symbol not in ("-", "not provided"):
                genes.update(s.strip() for s in symbol.split(";") if s.strip())
    return genes


def main() -> None:
    symbols = ensembl_to_symbol()
    disease = pathogenic_genes()
    print(f"annotation genes: {len(symbols):,}   ClinVar pathogenic genes: {len(disease):,}")

    rows = []
    for axis, path in AXES.items():
        frame = pd.read_csv(path, sep="\t")
        frame["symbol"] = frame["gene"].map(symbols)
        frame = frame.dropna(subset=["symbol"])
        flagged = (
            frame["qval"].lt(QVAL_GATE)
            & frame["effect_equal_donor"].ge(EFFECT_GATE)
        )
        has_disease = frame["symbol"].isin(disease)
        table = [
            [int((flagged & has_disease).sum()), int((flagged & ~has_disease).sum())],
            [int((~flagged & has_disease).sum()), int((~flagged & ~has_disease).sum())],
        ]
        odds, pval = fisher_exact(table, alternative="greater")
        n_flagged = int(flagged.sum())
        rows.append({
            "axis": axis,
            "genes_tested": int(len(frame)),
            "genes_flagged": n_flagged,
            "flagged_with_pathogenic_variant": table[0][0],
            "flagged_pct": round(100 * table[0][0] / n_flagged, 1) if n_flagged else None,
            "background_with_pathogenic_variant": table[1][0],
            "background_pct": round(100 * table[1][0] / max(sum(table[1]), 1), 1),
            "odds_ratio": round(float(odds), 3),
            "p_value": float(pval),
        })
        print(f"{axis:28s} {n_flagged:>5,} flagged, "
              f"{table[0][0]:>4,} with a pathogenic variant "
              f"({rows[-1]['flagged_pct']}% vs {rows[-1]['background_pct']}% background), "
              f"OR={odds:.2f}, P={pval:.3g}")

    # Sensitivity: within the tested set, a gene with more features has more
    # power to be flagged AND more chance of being well enough studied to carry
    # ClinVar submissions. Stratify on that and combine with Mantel-Haenszel;
    # if the association is a power artefact it collapses here.
    strat_rows = []
    for axis, path in AXES.items():
        frame = pd.read_csv(path, sep="\t")
        frame["symbol"] = frame["gene"].map(symbols)
        frame = frame.dropna(subset=["symbol"])
        flagged = (frame["qval"].lt(QVAL_GATE)
                   & frame["effect_equal_donor"].ge(EFFECT_GATE))
        if not flagged.any():
            continue
        has_disease = frame["symbol"].isin(disease)
        strata = pd.qcut(frame["n_features"].rank(method="first"), 5, labels=False)
        num = den = 0.0
        for s in sorted(strata.unique()):
            m = strata == s
            a = int((m & flagged & has_disease).sum())
            b = int((m & flagged & ~has_disease).sum())
            c = int((m & ~flagged & has_disease).sum())
            d = int((m & ~flagged & ~has_disease).sum())
            n = a + b + c + d
            if n == 0:
                continue
            num += a * d / n
            den += b * c / n
        mh = num / den if den else float("nan")
        strat_rows.append({"axis": axis, "mh_odds_ratio": round(mh, 3), "strata": 5})
        print(f"{axis:28s} Mantel-Haenszel OR (5 coverage strata) = {mh:.2f}")
    pd.DataFrame(strat_rows).to_csv(
        ROOT / "tables/clinvar_enrichment_stratified.tsv", sep="\t", index=False)

    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "tables/clinvar_enrichment.tsv", sep="\t", index=False)
    (ROOT / "tables/clinvar_enrichment.json").write_text(
        json.dumps({"clinvar_pathogenic_genes": len(disease),
                    "effect_gate": EFFECT_GATE, "qval_gate": QVAL_GATE,
                    "results": rows}, indent=2) + "\n"
    )
    print(f"\nwrote {ROOT / 'tables/clinvar_enrichment.tsv'}")


if __name__ == "__main__":
    main()
