#!/usr/bin/env python3
"""NAR Fig. 2 — Reliability of stored evidence (quality certification, not discovery)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

OUT = "NAR_Fig2"


def panel_a_junction(ax):
    j = pd.read_csv(S.FIGDATA / "fig4_cross_study_junction.tsv", sep="\t")
    ax.plot(j.n_studies, j.gtag_rate * 100, "-o", color=S.TEAL, lw=1.2,
            ms=4.5, markerfacecolor=S.TEAL, markeredgecolor="white", markeredgewidth=0.5)
    # size annotation for n_junctions at ends
    for _, r in j.iloc[[0, -1]].iterrows():
        ax.annotate(
            f"n={int(r.n_junctions):,}",
            (r.n_studies, r.gtag_rate * 100),
            textcoords="offset points", xytext=(6, -10 if r.n_studies == 1 else 6),
            fontsize=5.5, color=S.SLATE,
        )
    ax.set_xlabel("Independent studies reporting the novel junction")
    ax.set_ylabel("GT–AG canonical rate (%)")
    ax.set_title("Stored novel junctions: recurrence tracks canonicity", pad=4)
    ax.set_ylim(45, 100)
    ax.set_xticks(j.n_studies)
    r1, r13 = j.iloc[0], j.iloc[-1]
    ax.text(
        0.98, 0.08,
        f"{r1.gtag_rate*100:.1f}% (1 study) → {r13.gtag_rate*100:.1f}% (13 studies)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=S.TEAL,
        fontweight="bold",
    )
    S.panel_label(ax, "a")
    S.style_ax(ax)


def panel_b_marrow(ax):
    m = pd.read_csv(S.FIGDATA / "fig4b_marrow_reproducibility.tsv", sep="\t")
    x, y = m.novel_frac_myeloma.values, m.novel_frac_ccus.values
    r_p, _ = pearsonr(x, y)
    r_s, _ = spearmanr(x, y)
    ax.scatter(x, y, s=42, color=S.CORAL, edgecolor="white", linewidth=0.6, zorder=3)
    for _, row in m.iterrows():
        ax.text(row.novel_frac_myeloma + 0.012, row.novel_frac_ccus,
                row.ct.replace(" cell", "").replace("Dendritic", "DC"),
                fontsize=5.3, color=S.SLATE, va="center")
    lo = min(x.min(), y.min()) - 0.05
    hi = max(x.max(), y.max()) + 0.08
    ax.plot([lo, hi], [lo, hi], ls="--", color=S.GREY, lw=0.7, zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Novel-isoform fraction (myeloma marrow)")
    ax.set_ylabel("Novel-isoform fraction (CCUS marrow)")
    ax.set_title("Cross-study agreement of stored cell-type profiles", pad=4)
    ax.text(
        0.04, 0.96, f"Pearson r = {r_p:.2f}\nSpearman ρ = {r_s:.2f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=6, color=S.INK,
        fontweight="bold",
    )
    S.panel_label(ax, "b")
    S.style_ax(ax)


def panel_c_debt(ax):
    d = pd.read_csv(S.FIGDATA / "fig1c_debt_saturation.tsv", sep="\t")
    summ = json.loads((S.FIGDATA / "fig1c_debt_summary.json").read_text())
    ax.plot(d.n_runs, d.cum_frac * 100, color=S.BLUE, lw=1.4)
    ax.fill_between(d.n_runs, d.cum_frac * 100, color=S.BLUE, alpha=0.12)
    ax.axhline(summ["observed_frac"] * 100, color=S.GREY, lw=0.5, ls=":")
    ax.set_xlabel("Cumulative reprocessed runs")
    ax.set_ylabel("Reference transcripts seen as FSM (%)")
    ax.set_title("Reference coverage still unsaturated (resource bound)", pad=4)
    ax.set_ylim(0, 100)
    ax.text(
        0.98, 0.18,
        f"observed {summ['observed_frac']*100:.1f}%\ndebt {summ['debt_frac']*100:.1f}%\n"
        f"({summ['n_runs']} runs / {summ['n_studies']} studies)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.8, color=S.SLATE,
    )
    S.panel_label(ax, "c")
    S.style_ax(ax)


def panel_d_precomputed(ax):
    """Precomputed analysis tables available as database content (inventory, not discovery)."""
    rows = []
    for label, path, kind in [
        ("Allelic balance\n(ASE tables)", S.F2DATA / "ase_interaction.tsv", "ase"),
        ("Isoform usage\n(DIU tables)", S.F2DATA / "diu_celltype.tsv", "diu"),
        ("poly(A) choice\n(APA tables)", S.F2DATA / "apa_celltype.tsv", "apa"),
    ]:
        df = pd.read_csv(path, sep="\t")
        rows.append({
            "label": label,
            "n": len(df),
            "n_sig": int(df.sig.sum()),
            "frac_p05": float((df.pval < 0.05).mean()),
            "kind": kind,
        })
    cal = pd.read_csv(S.FIGDATA / "calibration_negcontrol.tsv", sep="\t")
    cal_map = {r.analysis: r.negctrl_frac_p05 for _, r in cal.iterrows()}

    x = np.arange(len(rows))
    n = np.array([r["n"] for r in rows], dtype=float)
    n_sig = np.array([r["n_sig"] for r in rows], dtype=float)
    colors = [S.GREY, S.TEAL, S.GOLD]

    ax.bar(x - 0.18, n / 1000, width=0.34, color=[c + "55" for c in colors],
           edgecolor="white", label="genes tested")
    ax.bar(x + 0.18, n_sig / 1000, width=0.34, color=colors,
           edgecolor="white", label="precomputed significant")
    ax.set_xticks(x)
    ax.set_xticklabels([r["label"] for r in rows], fontsize=5.8)
    ax.set_ylabel("Genes (×10³)")
    ax.set_title("Precomputed multi-layer tables (content inventory)", pad=4)
    ax.set_ylim(0, max(n.max() / 1000 * 1.28, 14))
    ax.legend(loc="upper right", fontsize=5.4, handlelength=0.9, frameon=False)

    # null calibration + ASE note below plot area via xlabel-style annotation
    notes = []
    for r in rows:
        key = {"ase": "ASE", "diu": "DIU", "apa": "APA"}[r["kind"]]
        nc = cal_map.get(key, np.nan)
        notes.append(f"{key} null {nc*100:.1f}%")
    ax.text(
        0.0, -0.28,
        "Permutation null (frac p<0.05): " + " · ".join(notes)
        + "  |  ASE sig = 0 (calibrated null table, not a discovery panel).",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.2, color=S.SLATE,
        clip_on=False,
    )
    # value labels on significant bars
    for xi, ns, nt in zip(x, n_sig, n):
        ax.text(xi + 0.18, ns / 1000 + 0.25, f"{int(ns):,}", ha="center",
                va="bottom", fontsize=5.2, color=S.INK)
        ax.text(xi - 0.18, nt / 1000 + 0.25, f"{int(nt):,}", ha="center",
                va="bottom", fontsize=5.0, color=S.SLATE)
    S.panel_label(ax, "d")
    S.style_ax(ax)


def main():
    fig = plt.figure(figsize=(183 * S.MM, 140 * S.MM))
    gs = fig.add_gridspec(
        2, 2, hspace=0.45, wspace=0.30,
        left=0.09, right=0.98, top=0.91, bottom=0.10,
    )
    panel_a_junction(fig.add_subplot(gs[0, 0]))
    panel_b_marrow(fig.add_subplot(gs[0, 1]))
    panel_c_debt(fig.add_subplot(gs[1, 0]))
    panel_d_precomputed(fig.add_subplot(gs[1, 1]))

    fig.suptitle(
        "Reliability of scTHREAD evidence contents",
        x=0.08, y=0.985, ha="left", fontsize=9, fontweight="bold", color=S.INK,
    )
    path = S.save(fig, OUT)
    print("wrote", path)


if __name__ == "__main__":
    main()
