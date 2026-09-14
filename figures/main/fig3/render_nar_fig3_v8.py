#!/usr/bin/env python3
"""Render scTHREAD NAR Figure 3 v9 — CD45 variable-exon evidence.

What changed from v8 and why
----------------------------
v8 panels c/d/e were built on ENST00000367364 and ENST00000697630, chosen
because they are the two most abundant PTPRC models. Auditing that choice
showed it does not survive:

* neither model carries any of the three CD45 variable exons, so the panels
  were not showing the CD45RA/RO switch the text claims;
* ENST00000367364 is a 3-exon, 1,017 bp fragment whose last exon reads 860 bp
  into intron 3 and terminates 64 kb short of the gene's 3' end, and
  ENST00000697630 is a `retained_intron` model terminating 14 kb short;
* only 17-27% of PTPRC long reads reach the canonical 3' end, and read 3' ends
  in the proximal zone show no poly(A) support above background (6.2-7.7%
  within 50 bp of a PolyASite 2.0 cluster against a 4.8% background), so those
  two models act as truncation sinks and their relative usage is confounded
  with read completeness rather than being transcript biology.

c/d/e now measure the splice choice directly. The three variable exons sit at
chr1:198.696-198.703 Mb while the median read 3' end is ~198.707 Mb, so nearly
every read spans the decision point regardless of where it stops — the
measurement is both better powered (465,308 informative reads over 31 runs and
2 studies, against 6,561 pair-assigned molecules at the isoform level) and
immune to the truncation artefact.

Figure contract
---------------
a, query entry and analysis status;
b, shared cell-type UMAP reference;
c-d, matched sashimi views of the exon-3 splice choice in B cells versus
     monocytes, drawn on identical axes;
e, that choice quantified across the four adequately powered cell types in two
   independent studies, with the gene-level statistics and replication audit;
f, junction support plus export/API access.

Guardrails
    Panels c-e report reads, not UMI-collapsed molecules. Only cell types with
    >=1,000 informative reads and >=3 runs in both studies are shown; Dendritic,
    Erythroid, Plasma cell and Progenitor are excluded as single-run or thin.
    The B-cell-versus-monocyte contrast replicates across both studies; the
    finer ordering among NK, T cell and monocyte does not, and four cell types
    cannot support a concordance statistic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "style"))
import nar_style as S
import render_nar_fig3_v3 as V3
import render_nar_fig3_v4 as V4
import render_nar_fig3_v5 as V5
import render_nar_fig3_v6 as V6
import render_nar_fig3_v7 as V7


REPO = Path(__file__).resolve().parents[3]
PROJECT = REPO
FIGURES = REPO / "figures"
TABLES = REPO / "source_data" / "tables"
ACCEPTOR_TABLE = TABLES / "figdata_cd45_acceptor_usage.tsv"
PER_RUN_TABLE = TABLES / "figdata_cd45_ro_per_run.tsv"
INCLUSION_MANIFEST = TABLES / "figdata_cd45_inclusion_manifest.json"

# Which cell types appear, and in what order, is a property of the data: it is
# decided by scripts/build_cd45_figdata.py from the read/run thresholds and
# recorded in the manifest. Hard-coding the list here would let the legend's
# "the four cell types carrying at least ..." drift out of being exhaustive.
MANIFEST = json.loads(INCLUSION_MANIFEST.read_text())
THRESHOLDS = MANIFEST["thresholds"]

WIDTH_MM = 183.0
HEIGHT_MM = 180.0

INK, SLATE, TEAL, BLUE = V3.INK, V3.SLATE, V3.TEAL, V3.BLUE
WHITE = V3.WHITE

# GENCODE v44, GRCh38. PTPRC + strand.
EXON3 = (198_692_347, 198_692_373)
EXON7 = (198_703_298, 198_703_372)
VARIABLE = {
    "A": (198_696_712, 198_696_909),
    "B": (198_699_564, 198_699_704),
    "C": (198_702_387, 198_702_530),
}
ACCEPTOR_X = {
    "exon A": VARIABLE["A"][0],
    "exon B": VARIABLE["B"][0],
    "exon C": VARIABLE["C"][0],
    "exon 7 (skip all)": EXON7[0],
}
ACCEPTOR_ORDER = ["exon A", "exon B", "exon C", "exon 7 (skip all)"]
ACCEPTOR_COLOR = {
    "exon A": "#2C6E9B",
    "exon B": "#6BA292",
    "exon C": "#C9A525",
    "exon 7 (skip all)": "#C1553B",
}
STUDIES = ("GSE307660", "GSE276974")
# The author asked for accessions rather than disease shorthand: a reader
# can look up GSE307660 and cannot look up "myeloma". The disease is kept
# in the caption below the arcs, so the mapping stays on the page.
STUDY_LABEL = {"GSE307660": "GSE307660", "GSE276974": "GSE276974"}
CELL_TYPE_ORDER = list(MANIFEST["selected_cell_types"])
VIEW = (198_691_600, 198_704_200)


def load_acceptors() -> pd.DataFrame:
    frame = pd.read_csv(ACCEPTOR_TABLE, sep="\t")
    expected = {"study", "cell_type", "acceptor", "reads", "frac", "total_reads", "runs"}
    if set(frame.columns) != expected:
        raise ValueError(f"Unexpected acceptor columns: {sorted(frame.columns)}")
    if set(frame.study) != set(STUDIES):
        raise ValueError(f"Unexpected studies: {sorted(set(frame.study))}")
    if set(frame.cell_type) != set(CELL_TYPE_ORDER):
        raise ValueError(f"Unexpected cell types: {sorted(set(frame.cell_type))}")
    return frame


def _acceptor_fracs(frame: pd.DataFrame, study: str, cell_type: str) -> dict[str, float]:
    rows = frame[frame.study.eq(study) & frame.cell_type.eq(cell_type)]
    return {row.acceptor: float(row.frac) for row in rows.itertuples()}


def _reads(frame: pd.DataFrame, study: str, cell_type: str) -> int:
    rows = frame[frame.study.eq(study) & frame.cell_type.eq(cell_type)]
    return int(rows.total_reads.iloc[0])


def _draw_sashimi(
    fig: plt.Figure,
    frame: pd.DataFrame,
    *,
    letter: str,
    x: float,
    cell_type: str,
) -> None:
    """One cell type's exon-3 splice choice, as arcs over the exon diagram."""
    V4._panel_heading(
        fig,
        letter,
        f"{cell_type}: exon-3 splice choice",
        x=x,
        y=0.705,
        title_size=5.95,
    )
    myeloma = _acceptor_fracs(frame, "GSE307660", cell_type)
    ccus = _acceptor_fracs(frame, "GSE276974", cell_type)
    fig.text(
        x,
        0.674,
        f"{_reads(frame,'GSE307660',cell_type):,} + "
        f"{_reads(frame,'GSE276974',cell_type):,} informative reads",
        fontsize=4.45,
        fontweight="bold",
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    ax = fig.add_axes([x, 0.445, 0.285, 0.215])
    ax.set_xlim(*VIEW)
    ax.set_ylim(-0.30, 1.06)
    ax.axis("off")

    # intron line + exon boxes
    ax.plot(VIEW, [0, 0], color="#9AA4AC", linewidth=0.5, zorder=1)
    boxes = [(EXON3, "3", "#8A949C"), (EXON7, "7", "#8A949C")]
    boxes += [(VARIABLE[k], k, ACCEPTOR_COLOR[f"exon {k}"]) for k in "ABC"]
    for (start, end), label, color in boxes:
        width = max(end - start, 90)          # keep 27 bp exons visible
        ax.add_patch(
            Rectangle(
                (start - (width - (end - start)) / 2, -0.085),
                width,
                0.17,
                facecolor=color,
                edgecolor="none",
                zorder=3,
            )
        )
        ax.text(
            (start + end) / 2,
            -0.175,
            label,
            ha="center",
            va="top",
            fontsize=4.1,
            fontweight="bold",
            color=color if label in "ABC" else INK,
            fontfamily=S._FAM,
        )

    donor = EXON3[1]
    for acceptor in ACCEPTOR_ORDER:
        frac_m = myeloma.get(acceptor, 0.0)
        frac_c = ccus.get(acceptor, 0.0)
        if max(frac_m, frac_c) < 0.01:
            continue
        target = ACCEPTOR_X[acceptor]
        height = 0.28 + 0.70 * max(frac_m, frac_c)
        mid = (donor + target) / 2
        path = MPath(
            [(donor, 0.02), (mid, height * 1.30), (target, 0.02)],
            [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3],
        )
        ax.add_patch(
            PathPatch(
                path,
                facecolor="none",
                edgecolor=ACCEPTOR_COLOR[acceptor],
                linewidth=0.45 + 3.6 * max(frac_m, frac_c),
                alpha=0.9,
                zorder=4,
            )
        )
        # Labels sit above their acceptor exon, not at the arc apex: apexes of
        # a dominant and a minor arc land close together and the bigger arc
        # draws straight through the smaller one's label.
        ax.text(
            target,
            0.125,
            f"{frac_m*100:.0f} / {frac_c*100:.0f}%",
            ha="center",
            va="bottom",
            fontsize=4.0,
            fontweight="bold",
            color=ACCEPTOR_COLOR[acceptor],
            fontfamily=S._FAM,
            zorder=6,
            bbox={"boxstyle": "square,pad=0.10", "facecolor": "white",
                  "edgecolor": "none", "alpha": 0.82},
        )
    fig.text(
        x,
        0.430,
        "Arc width = share of reads. Labels GSE307660 (myeloma) / GSE276974 (CCUS)",
        fontsize=4.1,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_sashimi_guardrail(fig: plt.Figure) -> None:
    fig.text(
        0.375,
        0.410,
        "Exons 3 and 7 are constitutive; A, B and C (198, 141 and 144 bp) are the exons whose "
        "inclusion defines CD45RA/RB/RC versus CD45RO.",
        fontsize=4.15,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.375,
        0.392,
        "Reads, not molecules. Measured at the splice junction, so unaffected by the 3'-truncation "
        "that confounds isoform-level PTPRC usage.",
        fontsize=4.15,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def _draw_exon_quantification(
    fig: plt.Figure,
    inputs: V3.Inputs,
    frame: pd.DataFrame,
    *,
    replication_text: str,
) -> None:
    V4._panel_heading(
        fig,
        "E",
        "Variable-exon usage across cell types and studies",
        x=0.055,
        y=0.355,
        title_size=6.25,
    )
    fig.text(
        0.055,
        0.324,
        f"CELL-TYPE DTU · q = {float(inputs.diu['qval']):.5f} "
        f"· effect = {float(inputs.diu['effect']):.3f} · 23 isoforms",
        fontsize=4.25,
        fontweight="bold",
        color=TEAL,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.305,
        replication_text,
        fontsize=4.2,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )

    ax = fig.add_axes([0.115, 0.100, 0.215, 0.180])
    labels: list[str] = []
    y = 0.0
    yticks: list[float] = []
    for cell_type in CELL_TYPE_ORDER:
        for study in STUDIES:
            fracs = _acceptor_fracs(frame, study, cell_type)
            left = 0.0
            for acceptor in ACCEPTOR_ORDER:
                value = fracs.get(acceptor, 0.0) * 100
                if value <= 0:
                    continue
                ax.barh(
                    y,
                    value,
                    left=left,
                    height=0.62,
                    color=ACCEPTOR_COLOR[acceptor],
                    edgecolor=WHITE,
                    linewidth=0.3,
                )
                if value >= 12:
                    ax.text(
                        left + value / 2,
                        y,
                        f"{value:.0f}",
                        ha="center",
                        va="center",
                        fontsize=3.9,
                        fontweight="bold",
                        color=WHITE,
                        fontfamily=S._FAM,
                    )
                left += value
            labels.append(f"{cell_type} · {STUDY_LABEL[study]}")
            yticks.append(y)
            y += 1.0
        y += 0.45
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels, fontsize=4.15)
    ax.tick_params(axis="y", length=0, pad=1.8)
    ax.set_ylim(y - 0.45, -0.7)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=4.0)
    ax.tick_params(axis="x", length=1.6, pad=1.2, width=0.4)
    ax.set_xlabel("% of informative reads", fontsize=4.2, labelpad=1.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.45)
    ax.spines["bottom"].set_color("#BBC3C8")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=ACCEPTOR_COLOR[a]) for a in ACCEPTOR_ORDER
    ]
    legend = ax.legend(
        handles,
        ["exon A", "exon B", "exon C", "skip all (RO)"],
        loc="lower left",
        bbox_to_anchor=(-0.02, 1.01),
        ncol=4,
        frameon=False,
        fontsize=4.0,
        handlelength=0.9,
        handleheight=0.75,
        columnspacing=0.9,
        handletextpad=0.35,
    )
    for text in legend.get_texts():
        text.set_fontfamily(S._FAM)

    per_run = pd.read_csv(PER_RUN_TABLE, sep="\t")
    strata = len(per_run)
    fig.text(
        0.055,
        0.072,
        f"B cells keep the variable exons in both studies (skip-all 0.4% and 7.3%);",
        fontsize=4.1,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.056,
        f"monocytes skip in ~30%. {strata} run × cell-type strata, each "
        f"≥{THRESHOLDS['min_stratum_reads']} reads.",
        fontsize=4.1,
        color=INK,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )
    fig.text(
        0.055,
        0.040,
        f"Shown: ≥{THRESHOLDS['min_cell_type_reads']:,} reads over "
        f"≥{THRESHOLDS['min_cell_type_runs']} runs in both studies. Four cannot "
        f"support a concordance test.",
        fontsize=4.1,
        color=SLATE,
        ha="left",
        va="center",
        fontfamily=S._FAM,
    )


def render(
    inputs: V3.Inputs,
    junctions: pd.DataFrame,
    frame: pd.DataFrame,
    stem: Path,
    dpi: int,
) -> list[Path]:
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
        fig = plt.figure(figsize=(WIDTH_MM * S.MM, HEIGHT_MM * S.MM), facecolor="none")
        V4._draw_query_panel(fig, inputs)
        V6._draw_cell_type_umap(
            fig,
            caption=f"{V5.DISPLAYED_CELLS:,} displayed cells · portal outlier filter on",
        )
        _draw_sashimi(fig, frame, letter="C", x=0.375, cell_type="B cell")
        _draw_sashimi(fig, frame, letter="D", x=0.695, cell_type="Monocyte")
        _draw_sashimi_guardrail(fig)
        _draw_exon_quantification(
            fig,
            inputs,
            frame,
            replication_text=V7.replication_text(),
        )
        V6._draw_junction_access(fig, junctions)

        outputs = [stem.with_suffix(ext) for ext in (".pdf", ".svg", ".png")]
        for output in outputs:
            fig.savefig(output, dpi=dpi, transparent=True, bbox_inches=None, pad_inches=0)
        plt.close(fig)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stem", type=Path, default=FIGURES / "NAR_Fig3_v9")
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--observed-suffix", default="")
    parser.add_argument("--base-api", type=Path, default=None)
    parser.add_argument("--snapshot", type=Path, default=V7.CORRECTED_SNAPSHOT)
    args = parser.parse_args()

    inputs = V7.load_corrected_inputs(args.observed_suffix, args.base_api)
    junctions = V6._load_junction_table()
    frame = load_acceptors()
    print(
        "V9 INPUT VALIDATION PASS",
        {
            "gene": inputs.api["gene"]["gid"],
            "diu_q": inputs.diu["qval"],
            "informative_reads": int(
                frame.groupby(["study", "cell_type"]).total_reads.first().sum()
            ),
            "cell_types": CELL_TYPE_ORDER,
            "studies": list(STUDIES),
        },
    )
    if args.validate_only:
        return
    args.snapshot.parent.mkdir(parents=True, exist_ok=True)
    args.snapshot.write_text(json.dumps(inputs.api, indent=2) + "\n")
    outputs = render(inputs, junctions, frame, args.stem.resolve(), args.dpi)
    for output in [args.snapshot, *outputs]:
        print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
