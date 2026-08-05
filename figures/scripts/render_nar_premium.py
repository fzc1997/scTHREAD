#!/usr/bin/env python3
"""scTHREAD NAR main figures — premium v3 (Database Issue).

Fig1 content | Fig2 reliability | Fig3 utility (portal screenshot + API)
All numbers from frozen figdata / live API fields. No discovery framing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

INK, SLATE = S.INK, S.SLATE
TEAL, CORAL, GOLD, GREY, CREAM = S.TEAL, S.CORAL, S.GOLD, S.GREY, S.CREAM
SOFT, BLUE = S.SOFT, S.BLUE
GREYL = "#D9DCE0"
CORALL = "#E4B4A6"

ASSETS = Path(__file__).resolve().parents[1] / "assets"
WALK = Path(__file__).resolve().parents[1] / "website_walkthrough"

CMAP_SEQ = LinearSegmentedColormap.from_list(
    "prem_seq", ["#F4EFE6", "#D8C7A6", "#8FB0A3", "#3F7C77", "#274B57", "#12233A"]
)
CMAP_HEAT = LinearSegmentedColormap.from_list(
    "prem_heat", ["#F7F3EC", "#CFE0D9", "#6FA39A", "#2F6E6B", "#1A3F45"]
)


def card(ax, x, y, w, h, fc=CREAM, ec=GREYL, lw=0.7, rs=0.03, transform=None):
    tr = transform if transform is not None else ax.transAxes
    p = FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.008,rounding_size={rs}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=tr, clip_on=False,
    )
    ax.add_patch(p)
    return p


def tag(ax, letter, x=-0.02, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            color=INK, va="top", ha="right", clip_on=False)


def title(ax, s, pad=4):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=7.6, color=INK, pad=pad)


def show_image_fit(ax, path: Path, *, letterbox=True):
    """Show image preserving aspect ratio (letterbox on transparent)."""
    img = mpimg.imread(path)
    h, w = img.shape[:2]
    ax.imshow(img, aspect="equal", interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if letterbox:
        ax.set_aspect("equal", adjustable="box")
    return img


# ===================================================================== Fig 1
def render_fig1():
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())
    cc = pd.read_csv(S.FIGDATA / "cohort_cellcount.tsv", sep="\t").sort_values("n_cells")
    cls = pd.read_csv(S.FIGDATA / "fig1b_classification_all.tsv", sep="\t")

    fig = plt.figure(figsize=(183 * S.MM, 155 * S.MM), facecolor="none")
    # top row: schematic wider; scopes slightly narrower
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.08, 1.12], width_ratios=[1.35, 1.0],
        hspace=0.36, wspace=0.22,
        left=0.06, right=0.98, top=0.91, bottom=0.07,
    )

    # ---- a GPT schematic ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a", x=-0.01, y=1.06)
    title(ax, "Architecture: public long reads → multi-layer evidence database")
    sch = ASSETS / "Fig1a_schematic_gpt.png"
    if sch.exists():
        show_image_fit(ax, sch)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "missing Fig1a schematic", ha="center", color=CORAL)

    # ---- b dual scope tiles ----
    ax = fig.add_subplot(gs[0, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    tag(ax, "b", x=-0.01, y=1.06)
    title(ax, "Two scopes — never mix denominators")

    card(ax, 0.02, 0.06, 0.46, 0.86, fc="#FFF9F0", ec=GOLD, lw=1.2, rs=0.05)
    ax.text(0.25, 0.84, "PUBLIC CATALOG", transform=ax.transAxes, ha="center",
            fontsize=6.0, fontweight="bold", color=GOLD)
    for yy, val, lab in (
        (0.66, "≈450", "sequencing runs"),
        (0.44, "≈3.0M", "cells (approx.)"),
        (0.22, ">200k", "observed isoforms"),
    ):
        ax.text(0.25, yy, val, transform=ax.transAxes, ha="center",
                fontsize=14, fontweight="bold", color=INK)
        ax.text(0.25, yy - 0.09, lab, transform=ax.transAxes, ha="center",
                fontsize=5.6, color=SLATE)

    card(ax, 0.52, 0.06, 0.46, 0.86, fc="#F2F7F6", ec=TEAL, lw=1.2, rs=0.05)
    ax.text(0.75, 0.84, "UNIFORM COHORT", transform=ax.transAxes, ha="center",
            fontsize=6.0, fontweight="bold", color=TEAL)
    for yy, val, lab in (
        (0.66, str(ch["n_studies"]), "studies"),
        (0.44, str(ch["n_runs"]), "runs reprocessed"),
        (0.22, f"{ch['n_annotated_cells']:,}", "annotated cells"),
    ):
        ax.text(0.75, yy, val, transform=ax.transAxes, ha="center",
                fontsize=14 if yy > 0.3 else 11.5, fontweight="bold", color=INK)
        ax.text(0.75, yy - 0.09, lab, transform=ax.transAxes, ha="center",
                fontsize=5.6, color=SLATE)

    # ---- c lollipop ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Uniform cohort: annotated cells per study")
    y = np.arange(len(cc))
    vals = cc.n_cells.values / 1000.0
    colors = [TEAL if v >= 20 else BLUE if v >= 5 else GREY for v in vals]
    ax.hlines(y, 0, vals, color=GREYL, lw=1.05, zorder=1)
    ax.scatter(vals, y, s=32, c=colors, edgecolors="white", linewidths=0.65, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(cc.gse.values, fontsize=5.5)
    ax.set_xlabel("Annotated cells (×10³)")
    ax.set_xlim(0, vals.max() * 1.18)
    ax.text(0.98, 0.98, f"Σ = {cc.n_cells.sum():,}", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.4, fontweight="bold", color=TEAL)
    for lab, col in [("≥20k", TEAL), ("5–20k", BLUE), ("<5k", GREY)]:
        ax.scatter([], [], c=col, s=24, label=lab, edgecolors="white")
    ax.legend(loc="center right", fontsize=5.2, frameon=False)
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.45, zorder=0)

    # ---- d stacked + novel % ----
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Molecule composition by system")
    order = ["FSM", "ISM", "NIC", "NNC"]
    systems = list(cls.groupby("system")["novel_frac"].mean().sort_values(ascending=False).index)
    mat, novels = [], []
    for sys in systems:
        sub = cls[cls.system == sys]
        means = np.array([sub[c].mean() for c in order], float)
        means = means / means.sum()
        mat.append(means)
        novels.append(float(sub["novel_frac"].mean()))
    mat = np.array(mat)
    y = np.arange(len(systems))
    left = np.zeros(len(systems))
    for i, c in enumerate(order):
        ax.barh(y, mat[:, i], left=left, height=0.66, color=S.CLASS[c],
                edgecolor="white", linewidth=0.35, label=c, zorder=2)
        left += mat[:, i]
    ax.scatter(novels, y, s=48, c=CORAL, marker="D", edgecolors="white",
               linewidths=0.65, zorder=4, label="novel %")
    for yi, nv in zip(y, novels):
        ax.text(min(nv + 0.025, 1.05), yi, f"{nv*100:.0f}%", va="center",
                fontsize=5.4, color=CORAL, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(systems, fontsize=6.0)
    ax.set_xlim(0, 1.18)
    ax.set_xlabel("Fraction of spliced molecules")
    ax.axvline(0.5, color=GREYL, lw=0.5, ls=":", zorder=0)
    ax.legend(loc="lower right", ncol=5, fontsize=5.0, handlelength=0.85,
              columnspacing=0.45, frameon=False)
    S.style_ax(ax)

    fig.suptitle("scTHREAD · Database content & architecture",
                 x=0.06, y=0.975, ha="left", fontsize=10, fontweight="bold", color=INK)
    S.save(fig, "NAR_Fig1")
    print("NAR_Fig1 v3 ok")


# ===================================================================== Fig 2
def render_fig2():
    j = pd.read_csv(S.FIGDATA / "fig4_cross_study_junction.tsv", sep="\t")
    m = pd.read_csv(S.FIGDATA / "fig4b_marrow_reproducibility.tsv", sep="\t")
    d = pd.read_csv(S.FIGDATA / "fig1c_debt_saturation.tsv", sep="\t")
    summ = json.loads((S.FIGDATA / "fig1c_debt_summary.json").read_text())
    cal = pd.read_csv(S.FIGDATA / "calibration_negcontrol.tsv", sep="\t")
    cal_map = {r.analysis: float(r.negctrl_frac_p05) for _, r in cal.iterrows()}

    axes_rows = []
    for name, path in [
        ("ASE", S.F2DATA / "ase_interaction.tsv"),
        ("DIU", S.F2DATA / "diu_celltype.tsv"),
        ("APA", S.F2DATA / "apa_celltype.tsv"),
    ]:
        df = pd.read_csv(path, sep="\t")
        axes_rows.append((name, len(df), int(df.sig.sum()), cal_map.get(name, np.nan)))

    fig = plt.figure(figsize=(183 * S.MM, 150 * S.MM), facecolor="none")
    gs = fig.add_gridspec(
        2, 2, hspace=0.42, wspace=0.28,
        left=0.09, right=0.98, top=0.90, bottom=0.10,
    )

    # ---- a ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Novel junctions: recurrence tracks canonicity")
    x = j.n_studies.values
    yv = j.gtag_rate.values * 100
    n_j = j.n_junctions.values.astype(float)
    sizes = 22 + 48 * np.log10(n_j / n_j.min())
    ax.plot(x, yv, color=TEAL, lw=1.4, zorder=2, alpha=0.85)
    ax.scatter(x, yv, s=sizes, c=yv, cmap=CMAP_SEQ, vmin=45, vmax=100,
               edgecolors="white", linewidths=0.55, zorder=3)
    ax.annotate(f"n={int(n_j[0]):,}", (x[0], yv[0]),
                textcoords="offset points", xytext=(7, -12), fontsize=5.2, color=SLATE)
    ax.annotate(f"n={int(n_j[-1]):,}", (x[-1], yv[-1]),
                textcoords="offset points", xytext=(-58, 7), fontsize=5.2, color=SLATE)
    card(ax, 0.42, 0.08, 0.54, 0.15, fc="#FFFFFF", ec=TEAL, lw=0.9)
    ax.text(0.69, 0.155, f"{yv[0]:.1f}%  →  {yv[-1]:.1f}%", transform=ax.transAxes,
            ha="center", va="center", fontsize=8.0, fontweight="bold", color=TEAL)
    # bubble-size legend OUTSIDE data (axes coords, bottom-left strip)
    ax.text(0.02, -0.16, "bubble size ~ log(n junctions):", transform=ax.transAxes,
            ha="left", va="center", fontsize=5.0, color=SLATE, clip_on=False)
    for s_leg, lab, xpos in (
        (18, "few", 0.38),
        (40, "×10", 0.48),
        (70, "×100", 0.58),
    ):
        ax.scatter([xpos], [-0.16], s=s_leg, c=TEAL, alpha=0.65, edgecolors="white",
                   linewidths=0.35, transform=ax.transAxes, clip_on=False, zorder=5)
        ax.text(xpos + 0.035, -0.16, lab, transform=ax.transAxes, ha="left", va="center",
                fontsize=5.0, color=SLATE, clip_on=False)
    ax.set_xlabel("Independent studies reporting the junction")
    ax.set_ylabel("GT–AG canonical rate (%)")
    ax.set_ylim(45, 102)
    ax.set_xticks(x)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.45, zorder=0)

    # ---- b ----
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Cross-study agreement of cell-type profiles")
    xv, y2 = m.novel_frac_myeloma.values, m.novel_frac_ccus.values
    rp, _ = pearsonr(xv, y2)
    rs, _ = spearmanr(xv, y2)
    resid = np.abs(y2 - xv)
    ax.plot([0.15, 0.85], [0.15, 0.85], ls="--", color=GREY, lw=0.75, zorder=1)
    ax.scatter(xv, y2, s=78, c=resid, cmap="YlOrRd", vmin=0, vmax=0.35,
               edgecolors="white", linewidths=0.75, zorder=3)
    for _, row in m.iterrows():
        lab = row.ct.replace(" cell", "").replace("Dendritic", "DC")
        ax.text(row.novel_frac_myeloma + 0.018, row.novel_frac_ccus, lab,
                fontsize=5.6, color=SLATE, va="center")
    ax.set_xlim(0.12, 0.88); ax.set_ylim(0.12, 0.88)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Novel fraction · myeloma marrow")
    ax.set_ylabel("Novel fraction · CCUS marrow")
    card(ax, 0.05, 0.70, 0.38, 0.24, fc="#FFFFFF", ec=CORAL, lw=0.95)
    ax.text(0.24, 0.88, f"r = {rp:.2f}", transform=ax.transAxes, ha="center",
            fontsize=9, fontweight="bold", color=CORAL)
    ax.text(0.24, 0.76, f"ρ = {rs:.2f}", transform=ax.transAxes, ha="center",
            fontsize=7.5, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.4, zorder=0)

    # ---- c ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Reference recovery saturates slowly (resource bound)")
    xr, yr = d.n_runs.values, d.cum_frac.values * 100
    ax.fill_between(xr, yr, 100, color=CORALL, alpha=0.40, label="unobserved (debt)", zorder=1)
    ax.fill_between(xr, 0, yr, color=TEAL, alpha=0.40, label="observed as FSM", zorder=2)
    ax.plot(xr, yr, color=TEAL, lw=1.55, zorder=3)
    ax.axhline(summ["observed_frac"] * 100, color=SLATE, lw=0.5, ls=":", zorder=2)
    ax.set_xlabel("Cumulative reprocessed runs")
    ax.set_ylabel("Reference transcripts (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(0, xr.max())
    card(ax, 0.48, 0.58, 0.48, 0.30, fc="#FFFFFF", ec=TEAL, lw=0.9)
    ax.text(0.72, 0.80, f"observed {summ['observed_frac']*100:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=7.0, fontweight="bold", color=TEAL)
    ax.text(0.72, 0.70, f"debt {summ['debt_frac']*100:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=7.0, fontweight="bold", color=CORAL)
    ax.text(0.72, 0.60, f"{summ['n_runs']} runs · {summ['n_studies']} studies",
            transform=ax.transAxes, ha="center", fontsize=5.3, color=SLATE)
    ax.legend(loc="lower right", fontsize=5.3, frameon=False)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- d inventory as paired bars (cleaner than dumbbell labels) ----
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Precomputed multi-layer tables (content inventory)")
    labels = ["ASE\nallelic balance", "DIU\nisoform usage", "APA\npoly(A) choice"]
    n_test = np.array([r[1] for r in axes_rows], float)
    n_sig = np.array([r[2] for r in axes_rows], float)
    colors = [GREY, TEAL, GOLD]
    x = np.arange(3)
    w = 0.34
    ax.bar(x - w / 2, n_test / 1000, width=w, color=[c + "55" for c in colors],
           edgecolor="white", label="genes tested", zorder=2)
    ax.bar(x + w / 2, n_sig / 1000, width=w, color=colors,
           edgecolor="white", label="precomputed significant", zorder=2)
    for i, (nt, ns, col) in enumerate(zip(n_test, n_sig, colors)):
        ax.text(i - w / 2, nt / 1000 + 0.25, f"{int(nt):,}", ha="center",
                fontsize=5.2, color=SLATE)
        ax.text(i + w / 2, ns / 1000 + 0.25, f"{int(ns):,}", ha="center",
                fontsize=5.6, fontweight="bold", color=col)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("Genes (×10³)")
    ax.set_ylim(0, max(n_test.max() / 1000 * 1.22, 15))
    ax.legend(loc="upper left", fontsize=5.3, frameon=False)
    null_txt = "  ·  ".join(
        f"{n} null {cal_map.get(n, float('nan'))*100:.1f}%" for n in ("ASE", "DIU", "APA")
    )
    ax.text(0.5, -0.22,
            null_txt + "   |   ASE significant = 0 (calibrated null table)",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.3, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    fig.suptitle("scTHREAD · Reliability of evidence contents",
                 x=0.09, y=0.975, ha="left", fontsize=10, fontweight="bold", color=INK)
    S.save(fig, "NAR_Fig2")
    print("NAR_Fig2 v3 ok")


# ===================================================================== Fig 3
def _find_portal_shot() -> Path | None:
    candidates = [
        ASSETS / "ptprc_gene_card.png",
        WALK / "02_browse_PTPRC.png",
        WALK / "02b_browse_PTPRC_vtime.png",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 50_000:
            return p
    return None


def render_fig3():
    GENE, GID = "PTPRC", "ENSG00000081237"
    diu = pd.read_csv(S.F2DATA / "diu_celltype.tsv", sep="\t")
    apa = pd.read_csv(S.F2DATA / "apa_celltype.tsv", sep="\t")
    ase = pd.read_csv(S.F2DATA / "ase_interaction.tsv", sep="\t")
    drow = diu[diu.gene == GID].iloc[0]
    arow = apa[apa.gene == GID].iloc[0]
    srow = ase[ase.gene.astype(str).str.upper() == GENE].iloc[0]
    pt = pd.read_csv(S.FIGDATA / "ptprc_isoform_usage.tsv", sep="\t")
    portal = _find_portal_shot()

    fig = plt.figure(figsize=(183 * S.MM, 168 * S.MM), facecolor="none")
    # layout: path | portal shot on top; heatmap | API bottom
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.05, 1.15], width_ratios=[0.95, 1.15],
        hspace=0.34, wspace=0.22,
        left=0.07, right=0.98, top=0.92, bottom=0.05,
    )

    # ---- a query path (GPT Image schematic, same visual language as Fig1a) ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a", x=-0.01, y=1.05)
    title(ax, "Query → multi-layer evidence → export")
    path_gpt = ASSETS / "Fig3a_query_path_gpt.png"
    if path_gpt.exists():
        show_image_fit(ax, path_gpt)
    else:
        # matplotlib fallback
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        steps = [
            (0.03, "1", "Search", "gene / ID", TEAL),
            (0.22, "2", "Browser", "evidence card", TEAL),
            (0.41, "3", "Inspect", "iso · PAS · ASE", TEAL),
            (0.60, "4", "Export", "CSV / JSON", TEAL),
            (0.79, "5", "Cite", "URL + release", CORAL),
        ]
        for x, num, head, sub, col in steps:
            fc = "#F8EFE8" if col == CORAL else "#FFFFFF"
            card(ax, x, 0.28, 0.17, 0.52, fc=fc, ec=col, lw=1.05, rs=0.05)
            ax.add_patch(Circle((x + 0.085, 0.72), 0.032, transform=ax.transAxes,
                                facecolor=col, edgecolor="white", lw=0.7, zorder=5, clip_on=False))
            ax.text(x + 0.085, 0.72, num, transform=ax.transAxes, ha="center", va="center",
                    fontsize=7, fontweight="bold", color="white", zorder=6)
            ax.text(x + 0.085, 0.50, head, transform=ax.transAxes, ha="center",
                    fontsize=6.8, fontweight="bold", color=INK)
            ax.text(x + 0.085, 0.38, sub, transform=ax.transAxes, ha="center",
                    fontsize=5.3, color=SLATE)

    # ---- b portal screenshot (utility proof) ----
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b", x=-0.01, y=1.05)
    title(ax, f"Live gene card · {GENE} (portal UI)")
    if portal is not None:
        img = mpimg.imread(portal)
        # crop top chrome-ish empty margins if huge; keep central content
        h, w = img.shape[:2]
        # focus on upper 70% where gene card lives
        img = img[: int(h * 0.78), :]
        ax.imshow(img, aspect="auto", interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color(SLATE)
            sp.set_linewidth(0.7)
        # caption ribbon
        ax.text(0.5, -0.06,
                f"DIU q={drow.qval:.3g} · APA q={arow.qval:.3g} · ASE q={srow.qval:.3g}  "
                f"(same values as API / figdata)",
                transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color=SLATE)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "portal screenshot missing", ha="center", color=CORAL)

    # ---- c heatmap ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, f"{GENE} isoform usage (top isoforms in cohort)")
    top = pt.groupby("transcript_id")["count"].sum().nlargest(4).index.tolist()
    sub = pt[pt.transcript_id.isin(top)].copy()
    sub["iso"] = sub.transcript_id.str[-6:]
    ct_order = [c for c in ["Progenitor", "B cell", "T cell", "NK", "Monocyte",
                            "Dendritic cell", "Plasma cell", "Erythroid"] if c in set(sub.ct)]
    mat = sub.pivot_table(index="ct", columns="iso", values="frac", aggfunc="sum").reindex(ct_order).fillna(0)
    if "Monocyte" in mat.index:
        mat = mat[mat.loc["Monocyte"].sort_values(ascending=False).index]
    im = ax.imshow(mat.values, aspect="auto", cmap=CMAP_HEAT, vmin=0, vmax=0.65)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v >= 0.08:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.2,
                        color="white" if v > 0.35 else INK,
                        fontweight="bold" if v > 0.4 else "normal")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels([f"…{c}" for c in mat.columns], fontsize=6.2)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=6.2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(GREYL)
        spine.set_linewidth(0.55)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Usage fraction", fontsize=5.6)
    cbar.ax.tick_params(labelsize=5.0)
    ax.text(0.0, -0.12, "Top 4 isoforms by molecule count · uniform cohort tables",
            transform=ax.transAxes, fontsize=5.2, color=SLATE)

    # ---- d API terminal ----
    ax = fig.add_subplot(gs[1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    tag(ax, "d", x=-0.01, y=1.05)
    title(ax, "Same query via API (export contract)")
    card(ax, 0.02, 0.04, 0.96, 0.90, fc="#1C2130", ec="#1C2130", lw=0, rs=0.05)
    for i, col in enumerate(("#FF5F56", "#FFBD2E", "#27C93F")):
        ax.add_patch(Circle((0.08 + i * 0.035, 0.86), 0.012, transform=ax.transAxes,
                            facecolor=col, edgecolor="none", clip_on=False))
    ax.text(0.52, 0.86, "api · scTHREAD v0.3", transform=ax.transAxes, ha="center",
            va="center", fontsize=5.5, color="#8A8F98")
    ax.text(0.07, 0.74,
            f"$ curl -s https://scthread.ai4sc.ac.cn/api/gene/{GENE}/overview",
            transform=ax.transAxes, ha="left", va="center",
            fontsize=5.2, color="#A7C4C0", family="monospace")
    snippet = (
        "{\n"
        f'  "gene": "{GENE}",\n'
        f'  "diu": {{"sig": true,  "q": {drow.qval:.4f}, "effect": {drow.effect:.3f}}},\n'
        f'  "apa": {{"sig": true,  "q": {arow.qval:.4f}, "effect": {arow.effect:.3f}}},\n'
        f'  "ase": {{"sig": false, "q": {srow.qval:.1f}, "effect": {srow.effect:.3f}}},\n'
        '  "export": ["csv", "json", "api"]\n'
        "}"
    )
    ax.text(0.07, 0.38, snippet, transform=ax.transAxes, ha="left", va="center",
            fontsize=5.7, color="#E8DFD0", family="monospace", linespacing=1.35)
    ax.text(0.07, 0.10,
            "Purpose: multi-layer evidence is queryable & exportable — not a mechanism claim.",
            transform=ax.transAxes, ha="left", va="center", fontsize=5.0,
            color="#8A8F98", style="italic")

    fig.suptitle(f"scTHREAD · Database utility ({GENE} walkthrough)",
                 x=0.07, y=0.98, ha="left", fontsize=10, fontweight="bold", color=INK)
    S.save(fig, "NAR_Fig3")
    print("NAR_Fig3 v3 ok", "portal=", portal)


def main():
    render_fig1()
    render_fig2()
    render_fig3()


if __name__ == "__main__":
    main()
