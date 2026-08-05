"""Shared style for scTHREAD NAR Database Issue figures (not GB Paper1)."""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("SCTHREAD_PROJECT_ROOT", "/gpfs/home/fuzc/project/scTHREAD"))
FIGDATA = ROOT / "results/paper1/figdata"
F2DATA = ROOT / "results/paper1/f2_grammar/figdata"
OUTDIR = ROOT / "NAR_database/figures"

_ARIAL = Path(os.environ.get("SCTHREAD_FONT", "/gpfs/home/fuzc/lib/Arial.ttf"))
if _ARIAL.exists():
    fm.fontManager.addfont(str(_ARIAL))
    _FAM = fm.FontProperties(fname=str(_ARIAL)).get_name()
else:
    _FAM = "DejaVu Sans"

MM = 1 / 25.4

# Editorial muted palette (colorblind-aware)
INK = "#1C2130"
SLATE = "#3D405B"
TEAL = "#2F6E6B"
CORAL = "#C1503A"
GOLD = "#C08A2E"
BLUE = "#3D6F9B"
GREY = "#8A8F98"
CREAM = "#F4F1EC"
SOFT = "#E8EEF0"

CLASS = {
    "FSM": "#3D6F9B",
    "ISM": "#8FB0C9",
    "NIC": "#D08A4C",
    "NNC": "#C1503A",
}

plt.rcParams.update({
    "font.family": _FAM,
    "font.size": 7,
    "axes.titlesize": 7.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.55,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 150,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "savefig.transparent": True,
    "figure.facecolor": "none",
    "axes.facecolor": "none",
    "axes.edgecolor": INK,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
})


def panel_label(ax, letter: str, x: float = -0.08, y: float = 1.08) -> None:
    """Nature: 8 pt bold upright lowercase a, b, c…"""
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=8, fontweight="bold", style="normal",
        family="sans-serif",
        va="top", ha="right", color=INK, clip_on=False,
    )


def save(fig, stem: str) -> Path:
    """Save PDF/SVG vectors plus a transparent 450-DPI PNG preview."""
    OUTDIR.mkdir(parents=True, exist_ok=True)
    base = OUTDIR / stem
    with matplotlib.rc_context({"savefig.bbox": None, "savefig.pad_inches": 0}):
        for ext in ("pdf", "svg", "png"):
            fig.savefig(
                f"{base}.{ext}",
                dpi=450,
                transparent=True,
                bbox_inches=None,
            )
    plt.close(fig)
    return Path(f"{base}.pdf")


def style_ax(ax) -> None:
    ax.set_axisbelow(True)
    ax.tick_params(length=2.2, width=0.55, pad=1.5)
    for sp in ax.spines.values():
        sp.set_color(INK)
        sp.set_linewidth(0.55)
