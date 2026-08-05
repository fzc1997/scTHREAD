#!/usr/bin/env python3
"""Render corrected scTHREAD NAR Figure 3 v7.

Core conclusion
---------------
The PTPRC portal workflow links queryable multi-layer evidence to corrected
biological-unit-aware DIU/APA/ASE statistics and an explicit study-stratified
replication audit.

Figure archetype
----------------
Asymmetric mixed-modality figure: portal localization and inspectability are
supporting evidence; corrected PTPRC inference and study stratification are
the quantitative evidence.  All rendering and export use Python/matplotlib.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import nar_style as S
import render_nar_fig3_v3 as V3
import render_nar_fig3_v4 as V4
import render_nar_fig3_v5 as V5
import render_nar_fig3_v6 as V6


PROJECT = Path(__file__).resolve().parents[2]
FIGURES = PROJECT / "figures"
TABLES = PROJECT / "tables"
P0 = TABLES / "p0_biological_unit_rerun"
REPLICATION = P0 / "study_stratified_replication"
LEGACY_API = TABLES / "PTPRC_api_overview.json"
CORRECTED_SNAPSHOT = TABLES / "PTPRC_api_overview_corrected_20260727.json"

WIDTH_MM = 183.0
HEIGHT_MM = 180.0


def _one(frame: pd.DataFrame, gene: str, label: str) -> pd.Series:
    rows = frame[frame.gene.astype(str).eq(gene)]
    if len(rows) != 1:
        raise ValueError(f"Expected one corrected {label} row for {gene}; found {len(rows)}")
    return rows.iloc[0]


def _normalized_usage_row(row: pd.Series, feature_name: str) -> pd.Series:
    return pd.Series(
        {
            "gene": row.gene,
            "pval": float(row.pval),
            "effect": float(row.effect_equal_donor),
            feature_name: int(row.n_features),
            "qval": float(row.qval),
            "sig": bool(row.sig),
        }
    )


def _normalized_ase_row(row: pd.Series) -> pd.Series:
    return pd.Series(
        {
            "gene": row.gene,
            "pval": float(row.pval),
            "effect": float(row.effect_equal_donor),
            "hi_ct": row.hi_celltype,
            "lo_ct": row.lo_celltype,
            "qval": float(row.qval),
            "sig": bool(row.sig),
        }
    )


def load_corrected_inputs(
    observed_suffix: str = "", base_api: Path | None = None
) -> V3.Inputs:
    diu_raw = _one(
        pd.read_csv(P0 / f"diu_observed{observed_suffix}.tsv", sep="\t"),
        V3.GID,
        "DIU",
    )
    apa_raw = _one(
        pd.read_csv(P0 / f"apa_observed{observed_suffix}.tsv", sep="\t"),
        V3.GID,
        "APA",
    )
    ase_raw = _one(
        pd.read_csv(P0 / f"ase_observed{observed_suffix}.tsv", sep="\t"),
        V3.GENE,
        "ASE",
    )
    diu = _normalized_usage_row(diu_raw, "n_iso")
    apa = _normalized_usage_row(apa_raw, "n_pas")
    ase = _normalized_ase_row(ase_raw)

    api_source = base_api or LEGACY_API
    api = copy.deepcopy(json.loads(api_source.read_text()))
    api["analyses"] = {
        "diu": diu.to_dict(),
        "apa": apa.to_dict(),
        "ase": ase.to_dict(),
    }
    api["analysis_provenance"] = {
        "schema_version": "scTHREAD.biological_unit_celltype_permutation.v1",
        "aggregation_unit": "study_id::donor_or_source_id",
        "eligible_runs": 25,
        "independent_sources": 19,
        "shared_celltypes": 7,
        "multiple_testing": "Benjamini-Hochberg",
        "effect": "maximum equal-source cell-type usage contrast",
        "effect_threshold": 0.20,
    }
    inputs = V3.Inputs(
        diu=diu,
        apa=apa,
        ase=ase,
        usage=pd.read_csv(V3.FIGDATA / "ptprc_isoform_usage.tsv", sep="\t"),
        api=api,
        capture=json.loads((V3.ASSETS / "ptprc_live_capture_v3.json").read_text()),
        portal_image=V3.ASSETS / "ptprc_gene_card_live_v3.png",
    )
    V3.validate_inputs(inputs)
    return inputs


def replication_text() -> str:
    tests = pd.read_csv(
        REPLICATION / "study_stratified_targeted_tests.tsv",
        sep="\t",
    )
    concordance = pd.read_csv(
        REPLICATION / "cross_study_concordance.tsv",
        sep="\t",
    )
    selected = tests[
        tests.analysis.eq("diu")
        & tests.gene_name.eq("PTPRC")
        & tests.eligible.astype(str).str.lower().eq("true")
    ].set_index("study")
    if set(selected.index) != {"GSE276974", "GSE307660"}:
        raise ValueError("PTPRC DIU is not eligible in both study-specific audits")
    agreement_rows = concordance[
        concordance.analysis.eq("diu") & concordance.gene.eq("PTPRC")
    ]
    if len(agreement_rows) != 1:
        raise ValueError("Expected one PTPRC DIU cross-study concordance row")
    agreement = agreement_rows.iloc[0]
    ccus = selected.loc["GSE276974"]
    myeloma = selected.loc["GSE307660"]
    return (
        f"Targeted DIU: CCUS P={ccus.pval:.3g}, Δ={ccus.effect_equal_donor:.2f} · "
        f"myeloma P={myeloma.pval:.3g}, Δ={myeloma.effect_equal_donor:.2f}\n"
        f"Signed cell-type contrasts: ρ={agreement.spearman_signed_contrasts:.2f} · "
        f"sign agreement={agreement.sign_agreement_fraction * 100:.0f}%"
    )


def render(inputs: V3.Inputs, junctions: pd.DataFrame, stem: Path, dpi: int) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(
        {
            "font.family": S._FAM,
            "font.size": 6.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": None,
            "savefig.pad_inches": 0,
            "axes.unicode_minus": False,
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
        }
    ):
        fig = plt.figure(
            figsize=(WIDTH_MM * S.MM, HEIGHT_MM * S.MM),
            facecolor="none",
        )
        V4._draw_query_panel(fig, inputs)
        V6._draw_cell_type_umap(fig)
        V6._draw_isoform_umap(
            fig,
            inputs,
            letter="c",
            x=0.375,
            image=V5.PRIMARY_CELL_MAP,
            isoform=V5.PRIMARY_ISOFORM,
            molecules="57,185",
            scale="0–5",
            usage_parts=(
                ("Monocyte", "Monocyte"),
                ("Plasma", "Plasma cell"),
                ("Progenitor", "Progenitor"),
            ),
            color=V6.TEAL,
        )
        V6._draw_isoform_umap(
            fig,
            inputs,
            letter="d",
            x=0.695,
            image=V5.SECONDARY_CELL_MAP,
            isoform=V5.SECONDARY_ISOFORM,
            molecules="25,697",
            scale="0–3",
            usage_parts=(
                ("Progenitor", "Progenitor"),
                ("B cell", "B cell"),
                ("Monocyte", "Monocyte"),
            ),
            color=V6.BLUE,
        )
        V6._draw_umap_guardrail(fig)
        V6._draw_dtu_quantification(
            fig,
            inputs,
            replication_text=replication_text(),
        )
        V6._draw_junction_access(fig, junctions)

        outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
        for output in outputs:
            fig.savefig(
                output,
                dpi=dpi,
                transparent=True,
                bbox_inches=None,
                pad_inches=0,
            )
        plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stem",
        type=Path,
        default=FIGURES / "NAR_Fig3_v7",
    )
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--observed-suffix",
        default="",
        help="Optional suffix for observed tables, for example _9999",
    )
    parser.add_argument(
        "--base-api",
        type=Path,
        default=None,
        help="Optional validated live API snapshot for coverage and gene metadata",
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=CORRECTED_SNAPSHOT,
    )
    args = parser.parse_args()

    inputs = load_corrected_inputs(args.observed_suffix, args.base_api)
    junctions = V6._load_junction_table()
    replication = replication_text()
    print(
        "CORRECTED INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "diu_q": inputs.diu["qval"],
            "diu_effect": inputs.diu["effect"],
            "apa_q": inputs.apa["qval"],
            "ase_q": inputs.ase["qval"],
            "replication": replication,
            "canvas_mm": [WIDTH_MM, HEIGHT_MM],
        },
    )
    if args.validate_only:
        return
    args.snapshot.parent.mkdir(parents=True, exist_ok=True)
    args.snapshot.write_text(json.dumps(inputs.api, indent=2) + "\n")
    outputs = render(inputs, junctions, args.stem.resolve(), args.dpi)
    for output in [args.snapshot, *outputs]:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
