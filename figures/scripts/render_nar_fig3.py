#!/usr/bin/env python3
"""NAR Fig. 3 — Reproducible query-to-export utility (product use case, not mechanism)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

OUT = "NAR_Fig3"
GENE = "PTPRC"
GID = "ENSG00000081237"


def _box(ax, x, y, w, h, title, body, fc, ec, title_fs=5.6, body_fs=5.6):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.01,rounding_size=0.03",
            facecolor=fc, edgecolor=ec, linewidth=0.85,
            transform=ax.transAxes, clip_on=False,
        )
    )
    if title:
        ax.text(x + w / 2, y + h - 0.06, title, transform=ax.transAxes,
                ha="center", va="top", fontsize=title_fs, fontweight="bold",
                color=ec, linespacing=1.15, clip_on=False)
        body_y = y + h * 0.38
    else:
        body_y = y + h * 0.50
    ax.text(x + w / 2, body_y, body, transform=ax.transAxes,
            ha="center", va="center", fontsize=body_fs, color=S.INK, linespacing=1.25)


def panel_a_query_flow(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    S.panel_label(ax, "a", x=-0.02, y=1.05)
    ax.set_title("Query-to-export path (external user, no reanalysis)", pad=6)

    steps = [
        (0.02, "1. Search\ngene symbol\nor ID"),
        (0.22, "2. Browser\nmulti-layer\nevidence card"),
        (0.42, "3. Inspect\nisoform / PAS\njunction / ASE"),
        (0.62, "4. Export\nCSV / JSON\nvia UI or API"),
        (0.82, "5. Cite\nrelease + URL\n+ accession"),
    ]
    for x, t in steps:
        _box(ax, x, 0.28, 0.16, 0.50, "", t, S.CREAM if x < 0.7 else "#F8EFE8",
             S.TEAL if x < 0.7 else S.CORAL)
    for x0 in (0.18, 0.38, 0.58, 0.78):
        ax.annotate(
            "", xy=(x0 + 0.04, 0.53), xytext=(x0, 0.53),
            xycoords=ax.transAxes, textcoords=ax.transAxes,
            arrowprops=dict(arrowstyle="-|>", color=S.SLATE, lw=0.85, mutation_scale=7),
        )
    ax.text(
        0.5, 0.10,
        f"Worked example: {GENE} · https://scthread.ai4sc.ac.cn/browse?query={GENE}",
        transform=ax.transAxes, ha="center", fontsize=5.8, color=S.SLATE, fontweight="bold",
    )


def panel_b_gene_card(ax):
    """Badge-style multi-evidence summary from figdata (matches live API)."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    S.panel_label(ax, "b", x=-0.02, y=1.05)
    ax.set_title(f"Gene card for {GENE}: multi-layer summaries returned together", pad=6)

    diu = pd.read_csv(S.F2DATA / "diu_celltype.tsv", sep="\t")
    apa = pd.read_csv(S.F2DATA / "apa_celltype.tsv", sep="\t")
    ase = pd.read_csv(S.F2DATA / "ase_interaction.tsv", sep="\t")
    d = diu[diu.gene == GID].iloc[0]
    a = apa[apa.gene == GID].iloc[0]
    s = ase[ase.gene.astype(str).str.upper() == GENE].iloc[0]

    # header
    ax.text(0.5, 0.90, f"{GENE}  ·  {GID}  ·  chr1:198,638,457–198,757,476",
            transform=ax.transAxes, ha="center", fontsize=6.5, fontweight="bold")

    cards = [
        (0.02, "DIU",
         f"{'Significant' if bool(d.sig) else 'Not sig.'}\neffect {d.effect:.3f}\nq = {d.qval:.3g}\n{int(d.n_iso)} isoforms",
         S.TEAL if bool(d.sig) else S.GREY),
        (0.27, "APA",
         f"{'Significant' if bool(a.sig) else 'Not sig.'}\neffect {a.effect:.3f}\nq = {a.qval:.3g}\n{int(a.n_pas)} PAS",
         S.GOLD if bool(a.sig) else S.GREY),
        (0.52, "ASE",
         f"{'Significant' if bool(s.sig) else 'Not sig.'}\neffect {s.effect:.3f}\nq = {s.qval:.3g}\nCT interaction",
         S.GREY),
        (0.77, "Junctions",
         "8,994 jct.\n2.34M mol.\n13 studies\n(API)",
         S.BLUE),
    ]
    for x, title, body, ec in cards:
        _box(ax, x, 0.10, 0.21, 0.70, title, body, "#FFFFFF", ec, title_fs=6.5, body_fs=5.5)

    ax.text(
        0.5, 0.02,
        "Same values exposed by GET /api/gene/PTPRC/overview (product reproducibility, not a mechanism claim).",
        transform=ax.transAxes, ha="center", fontsize=5.3, color=S.SLATE, style="italic",
    )


def panel_c_isoform_usage(ax):
    pt = pd.read_csv(S.FIGDATA / "ptprc_isoform_usage.tsv", sep="\t")
    # keep top 4 isoforms by total count
    top = pt.groupby("transcript_id")["count"].sum().nlargest(4).index.tolist()
    pt = pt[pt.transcript_id.isin(top)].copy()
    # short labels (last 6 digits)
    pt["iso"] = pt.transcript_id.str[-6:]
    ct_order = ["Progenitor", "B cell", "T cell", "NK", "Monocyte",
                "Dendritic cell", "Plasma cell", "Erythroid"]
    ct_order = [c for c in ct_order if c in set(pt.ct)]
    mat = pt.pivot_table(index="ct", columns="iso", values="frac", aggfunc="sum").reindex(ct_order).fillna(0)
    # order isoforms by monocyte preference for visual switch
    if "Monocyte" in mat.index:
        mat = mat[mat.loc["Monocyte"].sort_values(ascending=False).index]

    im = ax.imshow(mat.values, aspect="auto", cmap="YlGnBu", vmin=0, vmax=max(0.6, mat.values.max()))
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right", fontsize=5.5)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.set_title(f"{GENE} isoform usage by cell type (database matrix)", pad=4)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Usage fraction", fontsize=5.8)
    cbar.ax.tick_params(labelsize=5.2)
    # highlight two dominant isoforms note
    ax.text(
        0.0, -0.22,
        "Top 4 isoforms by molecule count; values from uniform-cohort isoform tables.",
        transform=ax.transAxes, fontsize=5.3, color=S.SLATE, ha="left",
    )
    S.panel_label(ax, "c")
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
        sp.set_color(S.GREY)


def panel_d_export_contract(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    S.panel_label(ax, "d", x=-0.02, y=1.05)
    ax.set_title("Export & API contract for the same query", pad=6)

    _box(
        ax, 0.04, 0.08, 0.92, 0.82, "",
        "",
        S.CREAM, S.SLATE,
    )
    lines = [
        ("Portal", f"https://scthread.ai4sc.ac.cn/browse?query={GENE}"),
        ("API", f"GET /api/gene/{GENE}/overview"),
        ("Layers", "expression · isoforms · poly(A) · ASE · junctions"),
        ("Download", "Browser Export CSV · /api/download/table/{{asset}}"),
        ("Scope", "Catalog ≈450 runs; cohort 15 / 67 / 452,213 cells"),
        ("Replay", "API fields match frozen DIU/APA/ASE figdata rows"),
    ]
    y = 0.78
    for k, v in lines:
        ax.text(0.09, y, k, transform=ax.transAxes, ha="left", va="center",
                fontsize=6.0, fontweight="bold", color=S.TEAL)
        ax.text(0.28, y, v, transform=ax.transAxes, ha="left", va="center",
                fontsize=5.6, color=S.INK)
        y -= 0.105
    ax.text(
        0.5, 0.12,
        "Purpose: multi-layer evidence is queryable and exportable.",
        transform=ax.transAxes, ha="center", fontsize=5.4, color=S.SLATE, style="italic",
    )


def main():
    fig = plt.figure(figsize=(183 * S.MM, 150 * S.MM))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[0.85, 1.15],
        hspace=0.38, wspace=0.28,
        left=0.07, right=0.98, top=0.93, bottom=0.08,
    )
    panel_a_query_flow(fig.add_subplot(gs[0, 0]))
    panel_b_gene_card(fig.add_subplot(gs[0, 1]))
    panel_c_isoform_usage(fig.add_subplot(gs[1, 0]))
    panel_d_export_contract(fig.add_subplot(gs[1, 1]))

    fig.suptitle(
        f"Database utility: reproducible multi-evidence query ({GENE})",
        x=0.07, y=0.99, ha="left", fontsize=9, fontweight="bold", color=S.INK,
    )
    path = S.save(fig, OUT)
    print("wrote", path)


if __name__ == "__main__":
    main()
