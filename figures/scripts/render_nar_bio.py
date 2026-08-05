#!/usr/bin/env python3
"""Biology-first NAR main figures (dense content).

Fig1 Biological landscape of the resource
Fig2 Biological structure of stored evidence
Fig3 Biological gene walkthrough (PTPRC/CD45)

All numbers from frozen figdata + sample_registry. No invented values.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nar_style as S

INK, SLATE = S.INK, S.SLATE
TEAL, CORAL, GOLD, GREY, CREAM = S.TEAL, S.CORAL, S.GOLD, S.GREY, S.CREAM
SOFT, BLUE = S.SOFT, S.BLUE
GREYL = "#D9DCE0"
CORALL = "#E4B4A6"
PURPLE = "#6B5B95"

ASSETS = Path(__file__).resolve().parents[1] / "assets"
REG = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong") + "/docs/sample_registry.tsv")
P0DATA = Path(__file__).resolve().parents[2] / "tables" / "p0_biological_unit_rerun"

# --- manuscript snapshot pin -------------------------------------------------
# sample_registry.tsv keeps growing (469 done runs on 27 Jul 2026 -> 488 today),
# so rendering from it silently drifts the figure off the manuscript denominators.
# Catalog panels therefore read the modality-audited frozen manifest by default;
# set SCTHREAD_CATALOG_SCOPE=live to reproduce the old live behaviour.
SNAPSHOT_MANIFEST = (
    Path(__file__).resolve().parents[2]
    / "tables" / "release_20260803" / "release_manifest.tsv"
)
SNAPSHOT_SHA256 = "d34b4bba637970d3571a983c7db6d21a778523274ed057fc42015c1dce1e97d2"
CATALOG_SCOPE = os.environ.get("SCTHREAD_CATALOG_SCOPE", "snapshot").lower()


def registry_done() -> pd.DataFrame:
    """Return the IsoQuant-complete catalog rows backing the Fig.1 panels.

    Default scope is the frozen release (453 runs, 34 studies) built under the
    single-cell long-read scope rule, which is a property of the data and not of
    the archive it came from. Frozen values win over the live registry;
    live-only columns (isoquant_route etc.) are joined back by run.

    Cell counts are NOT summed from these rows: the authority is per study, in
    study_cells(). Fourteen studies take an author-supplied or STARsolo total
    that has no per-run decomposition, so a per-run sum understates or
    overstates them.
    """
    live = pd.read_csv(REG, sep="\t", dtype=str)
    live = live[live["isoquant_status"].fillna("").str.lower().eq("done")].copy()
    if CATALOG_SCOPE == "live":
        return live

    digest = hashlib.sha256(SNAPSHOT_MANIFEST.read_bytes()).hexdigest()
    if digest != SNAPSHOT_SHA256:
        raise SystemExit(
            f"snapshot manifest checksum changed: {digest} != {SNAPSHOT_SHA256}"
        )
    snap = pd.read_csv(SNAPSHOT_MANIFEST, sep="\t", dtype=str)
    extra = [c for c in live.columns if c not in snap.columns] + ["srr"]
    out = snap.merge(live[extra].drop_duplicates("srr"), on="srr", how="left")
    out["isoquant_status"] = "done"
    return out


def catalog_scope_note() -> str:
    return (
        "frozen release 20260803 (single-cell long-read scope)"
        if CATALOG_SCOPE != "live"
        else "live sample_registry"
    )

CMAP_HEAT = LinearSegmentedColormap.from_list(
    "bio_heat", ["#F7F3EC", "#CFE0D9", "#6FA39A", "#2F6E6B", "#1A3F45"]
)
CMAP_NOVEL = LinearSegmentedColormap.from_list(
    "bio_novel", ["#F4EFE6", "#E8C48A", "#C1503A", "#7A2E1E"]
)
CMAP_COV = LinearSegmentedColormap.from_list(
    "bio_cov", ["#F0F0F0", "#A7C4C0", "#2F6E6B", "#12233A"]
)


def tag(ax, letter, x=-0.02, y=1.08):
    """Nature panel label: 8 pt bold upright."""
    S.panel_label(ax, letter, x=x, y=y)


def title(ax, s, pad=3.5):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=7.4, color=INK, pad=pad)


def card(ax, x, y, w, h, fc=CREAM, ec=GREYL, lw=0.7, rs=0.03):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.008,rounding_size={rs}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=ax.transAxes, clip_on=False,
    ))


def show_img(ax, path: Path):
    img = mpimg.imread(path)
    ax.imshow(img, aspect="equal", interpolation="lanczos")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def _mode_or_blank(s: pd.Series) -> str:
    s = s.dropna().astype(str)
    s = s[s.str.len() > 0]
    return s.mode().iloc[0] if len(s) else ""


# Tissue / system maps for ALL isoquant-done GSE (catalog atlas)
# Prefer biological systems only — no "Benchmark"/"Methods" display category.
TISSUE_MAP = {
    "GSE307660": "Myeloma marrow",
    "GSE276974": "CCUS marrow",
    "GSE288222": "Adult heart",
    "GSE178175": "Frontal cortex",
    "GSE274249": "iPSC→cortical neuron",
    "GSE283629": "iNeuron trajectory",
    "GSE283658": "Fibroblast→HSC",
    "GSE289790": "Prostate cancer",
    "GSE301658": "Pediatric glioma",
    "GSE303762": "Cross-platform BM mix",
    "GSE212945": "scNanoGPS multi-sample",
    "GSE248118": "Human ovary (Isosceles)",
    "GSE224045": "sc long-read (misc)",
    "GSE130708": "Mouse brain (ONT)",
    "GSE292324": "Splicing dynamics (AML-like)",
    "GSE314176": "Mouse brain translation",
    "GSE255520": "Mouse retina",
    "GSE295352": "Mouse brain long-read",
    "GSE295353": "Mouse long-read (ONT)",
    "GSE295932": "Human sc long-read",
    "GSE309071": "Mouse endothelium (SS2)",
    "GSE114157": "Mouse SS2 wells",
    "GSE140890": "Mouse SS2 plate",
    "GSE185554": "Mouse SS2 plate",
    "GSE274527": "Mouse SS2 plate",
    "GSE214231": "CAR-T engineering",
    "GSE250381": "Mouse preimplantation embryo",
    "GSE252344": "Human PacBio",
    "GSE252416": "Genotype–transcript",
    "GSE289428": "TLDR-seq full-length RNA",
    "GSE76026": "Mouse PacBio wells",
    # own / commercial (map for atlas when isoquant completes)
    "OWN_ASE_scONT": "Mouse gastrulation (scONT E6.5–E8.5)",
    "benagen_ovary_Young2": "Human ovary (snap-frozen nuclei)",
    "benagen_human_ovary": "Human ovary (snap-frozen nuclei)",
    "GSA_mouse_testis_11week": "Mouse testis (11 week)",
    "GSE158450": "Mouse hippocampus",
}
SYSTEM_MAP = {
    "GSE307660": "Blood/marrow", "GSE276974": "Blood/marrow", "GSE303762": "Blood/marrow",
    "GSE292324": "Blood/marrow",
    "GSE178175": "Brain", "GSE274249": "Brain", "GSE283629": "Brain",
    "GSE314176": "Brain", "GSE255520": "Brain", "GSE295352": "Brain", "GSE130708": "Brain",
    "GSE288222": "Heart",
    "GSE289790": "Cancer", "GSE301658": "Cancer",
    # Reproductive (human ovary)
    "GSE248118": "Reproductive",
    "benagen_ovary_Young2": "Reproductive",
    # datasets that entered under the single-cell long-read scope rule
    "benagen_human_ovary": "Reproductive",
    "GSA_mouse_testis_11week": "Reproductive",
    "GSE158450": "Brain",
    # Embryo / gastrulation (mouse)
    "GSE250381": "Embryo",
    "OWN_ASE_scONT": "Embryo",
    # former method-only accessions: rehome to biology when known, else Other
    "GSE212945": "Other", "GSE224045": "Other",
    "GSE289428": "Other", "GSE295932": "Other", "GSE295353": "Other",
    "GSE283658": "Differentiation",
    "GSE309071": "Differentiation",
    "GSE114157": "Smart-seq2", "GSE140890": "Smart-seq2", "GSE185554": "Smart-seq2",
    "GSE274527": "Smart-seq2", "GSE76026": "Smart-seq2",
    "GSE214231": "Other", "GSE252344": "Other", "GSE252416": "Other",
}

BIO_GROUP_TO_SYSTEM = {
    "组织-分化轨迹": "Differentiation",
    "疾病-临床": "Disease/clinical",
    "脑-神经-翻译": "Brain",
    # do not surface as "Methods/Benchmark" on figures
    "方法-Benchmark": "Other",
    "剪接-APA-调控": "Splicing/APA",
    "非RNA模态-跨物种": "Other",
    "forward_cross_C57xDBA": "Embryo",
}


def load_catalog_bio() -> pd.DataFrame:
    """IsoQuant-complete catalog runs (snapshot-pinned) → per-GSE atlas."""
    done = registry_done()
    for c in ("isoquant_cells", "isoquant_transcripts", "annotation_cells", "isoquant_barcodes_raw"):
        if c in done.columns:
            done[c] = pd.to_numeric(done[c], errors="coerce")

    agg = done.groupby("gse").agg(
        n_runs=("srr", "nunique"),
        isoquant_cells=("isoquant_cells", "sum"),
        isoquant_transcripts=("isoquant_transcripts", "sum"),
        species=("species", _mode_or_blank),
        platform=("platform", _mode_or_blank),
        biology_group=("biology_group", _mode_or_blank),
        description=("description", lambda s: s.dropna().astype(str).iloc[0] if len(s.dropna()) else ""),
        route=("isoquant_route", _mode_or_blank),
        method=("isoquant_cells_method", _mode_or_blank),
    ).reset_index()

    # English tissue / system
    def resolve_system(r):
        if r.gse in SYSTEM_MAP:
            return SYSTEM_MAP[r.gse]
        bg = r.biology_group or ""
        if bg in BIO_GROUP_TO_SYSTEM:
            # refine disease-clinical by description
            if bg == "疾病-临床":
                d = (r.description or "").lower()
                if "骨髓" in (r.description or "") or "marrow" in d or "ccus" in d or "骨髓瘤" in (r.description or ""):
                    return "Blood/marrow"
                if "癌" in (r.description or "") or "瘤" in (r.description or "") or "cancer" in d or "glioma" in d:
                    return "Cancer"
                if "心" in (r.description or "") or "heart" in d:
                    return "Heart"
            return BIO_GROUP_TO_SYSTEM[bg]
        if str(r.species).lower() == "mouse":
            return "Mouse"
        return "Other"

    agg["tissue"] = agg.gse.map(TISSUE_MAP)
    # fill tissue from description when missing
    miss = agg.tissue.isna() | (agg.tissue.astype(str).str.len() == 0)
    agg.loc[miss, "tissue"] = agg.loc[miss, "description"].where(
        agg.loc[miss, "description"].astype(str).str.len() > 0, agg.loc[miss, "gse"]
    )
    agg["system"] = agg.apply(resolve_system, axis=1)

    # optional layer coverage for labeled subset studies
    cov_path = S.FIGDATA / "ed7_coverage_matrix.tsv"
    if cov_path.exists():
        cov = pd.read_csv(cov_path, sep="\t")
        agg = agg.merge(cov, on="gse", how="left")
    lab_path = S.FIGDATA / "ed1_labeled_frac_bystudy.tsv"
    if lab_path.exists():
        lab = pd.read_csv(lab_path, sep="\t").rename(columns={"labeled_frac": "labeled_frac_all"})
        agg = agg.merge(lab, on="gse", how="left")
    ed1_path = S.FIGDATA / "ed1_cohort_qc.tsv"
    if ed1_path.exists():
        ed1 = pd.read_csv(ed1_path, sep="\t")
        keep = [c for c in ["gse", "novel_pct", "labeled_pct", "admitted_reads_M"] if c in ed1.columns]
        if len(keep) > 1:
            agg = agg.merge(ed1[keep], on="gse", how="left")

    agg["n_cells"] = agg["isoquant_cells"].fillna(0).astype(int)
    return agg.sort_values("n_cells", ascending=False)


def load_study_bio() -> pd.DataFrame:
    """Backward-compatible alias: labeled-cohort merge when needed elsewhere."""
    return load_catalog_bio()


STUDY_CELLS = SNAPSHOT_MANIFEST.parent / "release_study_cells.tsv"
# The cell authority is pinned like the run manifest: it is the only source of
# the headline cell count, so a silent edit here would change every figure.
STUDY_CELLS_SHA256 = "cb67bfffbbc87d8d4840c9e63b266a87d6398393a8be9fdac2b38363449af50d"


def study_cells() -> pd.DataFrame:
    """Authoritative per-study cell counts for the pinned release.

    The winning method differs by study (author total, STARsolo call, pipeline
    caller, one well per cell), so this table - not a sum over run rows - is the
    cell authority. `cells_agree_with_pipeline` records where the two coincide.
    """
    digest = hashlib.sha256(STUDY_CELLS.read_bytes()).hexdigest()
    if CATALOG_SCOPE != "live" and digest != STUDY_CELLS_SHA256:
        raise SystemExit(
            f"study cell table checksum changed: {digest} != {STUDY_CELLS_SHA256}"
        )
    return pd.read_csv(STUDY_CELLS, sep="\t")


def catalog_headline() -> dict:
    """Headline KPIs for the pinned catalog scope (snapshot by default)."""
    done = registry_done()
    study = study_cells()
    per_run = done.assign(
        c=pd.to_numeric(done.get("isoquant_cells"), errors="coerce").fillna(0)
    )
    # One study spans both species (human cortex plus a mouse hippocampus run,
    # confirmed against ENA), and its counts are per-run reproducible, so it is
    # split by run; every other study is assigned whole.
    cells = {"human": 0, "mouse": 0}
    for _, row in study.iterrows():
        if row["species"] in cells:
            cells[row["species"]] += int(row["cells"])
        else:
            sub = per_run[per_run["gse"] == row["gse"]]
            for sp, value in sub.groupby("species")["c"].sum().items():
                if sp in cells:
                    cells[sp] += int(value)
    species = done["species"].fillna("")
    return {
        "n_runs": int(done["srr"].nunique()),
        "n_studies": int(study["gse"].nunique()),
        "n_cells": int(study["cells"].sum()),
        "n_human": int((species == "human").sum()),
        "n_mouse": int((species == "mouse").sum()),
        "cells_human": cells["human"],
        "cells_mouse": cells["mouse"],
        "n_ont": int((done["platform"].fillna("") == "ONT").sum()),
        "n_pacbio": int((done["platform"].fillna("") == "PacBio").sum()),
    }


# ===================================================================== Fig 1
def render_fig1(stem: str = "NAR_Fig1"):
    """Fig.1 content & construction (de-overlap with GA / Fig.2 evidence).

    Focus: catalog composition, construction path, study atlas, species×platform.
    Avoid: portal walkthrough (GA/Fig.3); novelty/GT–AG (Fig.2); heavy dual-scope
    "do not mix denominators" framing (cell-type labels are progress, not a product claim).
    """
    cat = catalog_headline()
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())  # labeled subset
    studies = load_catalog_bio()
    markers = pd.read_csv(S.FIGDATA / "ed5_marker_heatmap.tsv", sep="\t")
    done = registry_done()
    done["isoquant_cells"] = pd.to_numeric(done.get("isoquant_cells"), errors="coerce")

    # Biological systems only — no Methods/Benchmark display alias
    sys_display = {}

    colors_sys = {
        "Blood/marrow": TEAL, "Brain": PURPLE, "Heart": CORAL, "Cancer": GOLD,
        "Differentiation": BLUE, "Reproductive": "#C75B8F", "Embryo": "#D4A017",
        "Smart-seq2": "#8B7355",
        "Disease/clinical": "#B85C38", "Splicing/APA": "#5B7C99",
        "Mouse": "#8B7355", "Other": "#A0A4AB",
    }

    fig = plt.figure(figsize=(183 * S.MM, 230 * S.MM), facecolor="none")
    gs = fig.add_gridspec(
        3, 2, height_ratios=[0.82, 1.40, 1.0],
        hspace=0.42, wspace=0.28,
        left=0.10, right=0.98, top=0.925, bottom=0.048,
    )

    # ---- a: biological systems (FULL catalog, isoquant_cells) ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Catalog composition by biological system")
    sys_order = [
        "Blood/marrow", "Brain", "Heart", "Cancer", "Differentiation",
        "Reproductive", "Embryo",
        "Smart-seq2", "Disease/clinical", "Splicing/APA", "Other", "Mouse",
    ]
    sys_order = [s for s in sys_order if s in set(studies.system)]
    for s in sorted(set(studies.system) - set(sys_order)):
        sys_order.append(s)
    rows = []
    for s in sys_order:
        sub = studies[studies.system == s]
        rows.append((s, int(sub.n_cells.sum()), int(sub.gse.nunique()), int(sub.n_runs.sum())))
    rows = sorted(rows, key=lambda r: r[1])  # small→large for barh
    y = np.arange(len(rows))
    cells = np.array([r[1] for r in rows], dtype=float) / 1000.0
    ax.barh(y, cells, color=[colors_sys.get(r[0], GREY) for r in rows],
            edgecolor="white", height=0.72, zorder=2)
    for yi, (name, ncell, ngse, nrun) in enumerate(rows):
        ax.text(cells[yi] + max(cells) * 0.02, yi,
                f"{ngse} stud. · {nrun} runs · {ncell:,} cells",
                va="center", fontsize=5.2, color=SLATE)
    ax.set_yticks(y)
    ax.set_yticklabels([sys_display.get(r[0], r[0]) for r in rows], fontsize=6.2)
    ax.set_xlabel("Cells (×10³)")
    ax.set_xlim(0, max(cells.max() * 1.55, 1))
    ax.text(0.98, 0.98,
            f"Catalog: {cat['n_studies']} datasets · {cat['n_runs']} runs · "
            f"{cat['n_cells']:,} cells",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
            fontweight="bold", color=TEAL)
    S.style_ax(ax)
    ax.grid(True, axis="x", color=GREYL, lw=0.4, zorder=0)

    # ---- b: construction path + catalog scale (no dual-scope warning tile) ----
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Uniform construction of the catalog")
    ax.axis("off")
    steps = [
        ("1 Public\nsc long-read", TEAL),
        ("2 Uniform\nIsoQuant", BLUE),
        ("3 Cell calling\n+ RNA layers", GOLD),
        ("4 Catalog\nrelease", CORAL),
    ]
    for i, (lab, col) in enumerate(steps):
        x0 = 0.02 + i * 0.245
        card(ax, x0, 0.68, 0.22, 0.26, fc="#FFFFFF", ec=col, lw=1.0, rs=0.04)
        ax.text(x0 + 0.11, 0.81, lab, transform=ax.transAxes, ha="center", va="center",
                fontsize=5.5, fontweight="bold", color=col, linespacing=1.2)
        if i < 3:
            ax.annotate("", xy=(x0 + 0.235, 0.81), xytext=(x0 + 0.215, 0.81),
                        xycoords=ax.transAxes, textcoords=ax.transAxes,
                        arrowprops=dict(arrowstyle="->", color=GREY, lw=0.9))

    card(ax, 0.02, 0.08, 0.96, 0.50, fc="#E8F2F1", ec=TEAL, lw=1.05, rs=0.03)
    ax.text(0.05, 0.48, "Single-cell long-read catalog, one denominator throughout",
            transform=ax.transAxes, fontsize=6.2, fontweight="bold", color=TEAL)
    ax.text(0.05, 0.28,
            f"{cat['n_runs']} runs · {cat['n_studies']} datasets  ·  "
            f"{cat['n_cells']:,} cells\n"
            f"Human {cat['n_human']} · Mouse {cat['n_mouse']} runs  ·  "
            f"ONT {cat['n_ont']} · PacBio {cat['n_pacbio']}\n"
            f"RNA layers: expression · isoform · poly(A) · allele",
            transform=ax.transAxes, fontsize=5.8, color=INK, va="top", linespacing=1.4)
    ax.text(0.05, 0.12,
            "Every cell carries a cell-type assignment: author labels where published, "
            "one uniform procedure otherwise.",
            transform=ax.transAxes, fontsize=5.3, color=SLATE, va="top")

    # ---- c: FULL catalog study atlas ----
    ax = fig.add_subplot(gs[1, :])
    tag(ax, "c", x=-0.01)
    title(ax, "Dataset atlas: biological system, platform and cells")
    st = studies.copy()
    st["sp_s"] = st.species.fillna("?").astype(str).str.lower().map(
        lambda s: "mouse" if "mouse" in s or s == "mm" else ("human" if "human" in s else s[:6])
    )
    st["label"] = st.apply(
        lambda r: f"{r.gse}  {str(r.tissue)[:28]}"
        + (f"  [{r.sp_s}]" if r.sp_s and r.sp_s != "human" else ""),
        axis=1,
    )
    # matrix: platform ONT/PacBio + optional layer presence
    layer_cols = [c for c in ["ASE", "DIU", "APA", "F2jct"] if c in st.columns]
    # build display matrix: col0=log10 cells scaled, then layers if any
    n = len(st)
    n_layer = len(layer_cols)
    # use cells as continuous bar via imshow of layers only if present; always show cells text
    if n_layer:
        mat = st[layer_cols].fillna(0).values.astype(float)
        mat_norm = mat / np.maximum(mat.max(axis=0, keepdims=True), 1)
        im = ax.imshow(mat_norm, aspect="auto", cmap=CMAP_COV, vmin=0, vmax=1,
                       extent=(-0.5, n_layer - 0.5, n - 0.5, -0.5))
        ax.set_xticks(range(n_layer))
        ax.set_xticklabels(["ASE", "DIU", "APA", "Junction"][:n_layer], fontsize=6.2)
        for i in range(n):
            for j in range(n_layer):
                v = int(mat[i, j])
                if v > 0:
                    ax.text(j, i, str(v), ha="center", va="center", fontsize=4.6,
                            color="white" if mat_norm[i, j] > 0.45 else INK, fontweight="bold")
        x_text = n_layer - 0.5 + 0.2
        x_right = n_layer + 5.5
    else:
        ax.set_xlim(-0.5, 6)
        x_text = 0.0
        x_right = 6.5
        ax.set_xticks([])

    ax.set_yticks(range(n))
    ax.set_yticklabels(st.label.values, fontsize=4.8)
    for i, (_, r) in enumerate(st.iterrows()):
        plat = str(r.platform) if pd.notna(r.platform) else "?"
        method = str(r.method)[:18] if pd.notna(r.method) else ""
        labf = r.get("labeled_frac_all", np.nan)
        lab_s = f"  label {labf*100:.0f}%" if pd.notna(labf) else ""
        ax.text(x_text + 0.15, i,
                f"{int(r.n_cells):,} cells · {int(r.n_runs)} runs · {plat}{lab_s}",
                va="center", ha="left", fontsize=4.7, color=SLATE)
        ax.plot(-0.75 if n_layer else -0.35, i, "s",
                color=colors_sys.get(r.system, GREY), ms=4.2, clip_on=False)
    ax.set_xlim((-1.0 if n_layer else -0.6), x_right)
    ax.set_ylim(n - 0.5, -0.5)
    if n_layer:
        ax.set_xlabel(
            "Precomputed layer coverage (run counts) for labeled-subset GSE only; blank = not available",
            fontsize=5.6,
        )
        # summary from ed7 (honest: only labeled GSE, not full catalog)
        cov_sum_path = Path(__file__).resolve().parents[2] / "tables/fig1_layer_coverage_summary.tsv"
        if cov_sum_path.exists():
            cs = pd.read_csv(cov_sum_path, sep="\t")
            bits = []
            for _, r in cs.iterrows():
                bits.append(f"{r.layer} {int(r.n_gse_with_any)}/{int(r.n_gse_total)}")
            ax.text(
                0.0, -0.12,
                "Labeled GSE with any precomputed runs: " + " · ".join(bits),
                transform=ax.transAxes, fontsize=5.0, color=SLATE, clip_on=False,
            )
    handles = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=colors_sys[s],
                          markersize=5.5, label=sys_display.get(s, s))
               for s in sys_order if s in colors_sys and s in set(st.system)]
    ax.legend(handles=handles, loc="lower right", ncol=4, fontsize=4.8, frameon=False,
              bbox_to_anchor=(1.0, -0.14))
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)

    # ---- d: species × platform (FULL catalog; unique vs GA) ----
    ax = fig.add_subplot(gs[2, 0])
    tag(ax, "d"); title(ax, "Species × platform composition")
    done2 = done.copy()
    done2["sp"] = done2["species"].fillna("?").str.lower().map(
        lambda s: "mouse" if "mouse" in s else ("human" if "human" in s else "other")
    )
    done2["plat"] = done2["platform"].fillna("?").astype(str)
    done2.loc[~done2["plat"].isin(["ONT", "PacBio"]), "plat"] = "other"
    piv_sp = (
        done2.groupby(["sp", "plat"])["isoquant_cells"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=["human", "mouse", "other"], columns=["ONT", "PacBio", "other"])
        .fillna(0)
    )
    # drop empty rows/cols
    piv_sp = piv_sp.loc[(piv_sp.sum(axis=1) > 0), [c for c in piv_sp.columns if piv_sp[c].sum() > 0]]
    x = np.arange(len(piv_sp.index))
    bottom = np.zeros(len(piv_sp.index))
    plat_cols = {"ONT": TEAL, "PacBio": PURPLE, "other": GREY}
    for plat in piv_sp.columns:
        vals = piv_sp[plat].values / 1000.0
        ax.bar(x, vals, bottom=bottom, color=plat_cols.get(plat, GREY),
               edgecolor="white", width=0.62, label=plat, zorder=2)
        for xi, v, b in zip(x, vals, bottom):
            if v > max(piv_sp.values.sum(axis=1).max() / 1000 * 0.06, 5):
                ax.text(xi, b + v / 2, f"{v:.0f}k", ha="center", va="center",
                        fontsize=5.4, color="white", fontweight="bold")
        bottom = bottom + vals
    ax.set_xticks(x)
    ax.set_xticklabels([s.capitalize() for s in piv_sp.index], fontsize=6.5)
    ax.set_ylabel("Cells (×10³)")
    ax.legend(loc="upper right", fontsize=5.5, frameon=False, title="platform", title_fontsize=5.2)
    ax.text(0.02, 0.98,
            f"Human {cat['n_human']} · Mouse {cat['n_mouse']} runs\n"
            f"ONT {cat['n_ont']} · PacBio {cat['n_pacbio']}",
            transform=ax.transAxes, ha="left", va="top", fontsize=5.3, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- e: marker validation (LABELED subset only — supports table denominators) ----
    ax = fig.add_subplot(gs[2, 1])
    tag(ax, "e"); title(ax, "Cell-type labels (marker check where annotations exist)")
    piv = markers.pivot(index="marker", columns="cell_type", values="logcpm")
    ct_ord = [c for c in ["Progenitor", "Erythroid", "Monocyte", "Dendritic cell",
                          "B cell", "Plasma cell", "T cell", "NK"] if c in piv.columns]
    mk_ord = [m for m in ["CD34", "GATA1", "HBB", "CD14", "LYZ", "MPO", "IRF8",
                          "MS4A1", "MZB1", "CD3D", "CD8A", "NKG7"] if m in piv.index]
    piv = piv.reindex(index=mk_ord, columns=ct_ord)
    im = ax.imshow(piv.values, aspect="auto", cmap=CMAP_HEAT)
    ax.set_xticks(range(len(ct_ord)))
    ax.set_xticklabels(ct_ord, rotation=40, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(mk_ord)))
    ax.set_yticklabels(mk_ord, fontsize=5.8)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("log CPM", fontsize=5.4)
    cbar.ax.tick_params(labelsize=4.8)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color(GREYL)
        sp.set_linewidth(0.5)
    ax.text(0.0, -0.22,
            "Mean log-CPM of canonical markers by cell type; label plausibility "
            "check, not a statistical test",
            transform=ax.transAxes, fontsize=5.1, color=SLATE)

    fig.suptitle(
        "scTHREAD · Content and construction of the single-cell long-read catalog",
        x=0.09, y=0.978, ha="left", fontsize=10, fontweight="bold", color=INK,
    )
    fig.text(0.09, 0.008,
             f"Catalog (a–d): {catalog_scope_note()} → "
             f"{cat['n_runs']} runs / {cat['n_studies']} datasets / "
             f"{cat['n_cells']:,} cells.  "
             f"Panel e: cell-type annotations where available.  "
             f"Evidence characterization → Fig.2; gene query walkthrough → Fig.3.",
             fontsize=5.1, color=SLATE)
    S.save(fig, stem)
    print(f"{stem} bio ok  catalog={cat['n_runs']} runs / {cat['n_cells']:,} cells")


# ===================================================================== Fig 2
def render_fig2(stem: str = "NAR_Fig2", observed_suffix: str = ""):
    """Evidence characterization / content inventory (NAR Fig.2 role).

    stem: output basename under figures/ (use NAR_Fig2_v2 to avoid clobbering).
    """
    j = pd.read_csv(S.FIGDATA / "fig4_cross_study_junction.tsv", sep="\t")
    m = pd.read_csv(S.FIGDATA / "fig4b_marrow_reproducibility.tsv", sep="\t")
    d = pd.read_csv(S.FIGDATA / "fig1c_debt_saturation.tsv", sep="\t")
    summ = json.loads((S.FIGDATA / "fig1c_debt_summary.json").read_text())
    validation = pd.read_csv(
        P0DATA / "validation_summary.tsv", sep="\t", keep_default_na=False
    )
    null_validation = validation[validation.dataset == "null"]
    null_counts = null_validation.groupby("analysis").size().to_dict()
    if null_counts != {"apa": 3, "ase": 3, "diu": 3}:
        raise ValueError(f"expected three null seeds per axis, got {null_counts}")
    null_ranges = {
        analysis.upper(): (
            float(group.raw_p_lt_0_05_fraction.min()),
            float(group.raw_p_lt_0_05_fraction.max()),
        )
        for analysis, group in null_validation.groupby("analysis")
    }
    novel_ct = pd.read_csv(S.FIGDATA / "fig2a_novel_by_celltype.tsv", sep="\t")
    cls = pd.read_csv(S.FIGDATA / "fig1b_classification_all.tsv", sep="\t")
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())

    axes_rows = []
    for name, path in [
        ("ASE", P0DATA / f"ase_observed{observed_suffix}.tsv"),
        ("DIU", P0DATA / f"diu_observed{observed_suffix}.tsv"),
        ("APA", P0DATA / f"apa_observed{observed_suffix}.tsv"),
    ]:
        df = pd.read_csv(path, sep="\t")
        axes_rows.append((name, len(df), int(df.sig.sum())))

    # slightly taller full-page multi-panel (NAR Database Issue density)
    fig = plt.figure(figsize=(183 * S.MM, 210 * S.MM), facecolor="none")
    gs = fig.add_gridspec(
        3, 2, height_ratios=[1.05, 1.0, 1.05],
        hspace=0.42, wspace=0.30,
        left=0.09, right=0.98, top=0.915, bottom=0.082,
    )

    # ---- a cell-type novelty summary across two diseases ----
    ax = fig.add_subplot(gs[0, 0])
    tag(ax, "a"); title(ax, "Cell-type novel-isoform fractions (two marrow studies)")
    marrow_tab = pd.read_csv(S.FIGDATA / "fig4b_marrow_reproducibility.tsv", sep="\t")
    marrow_tab = marrow_tab.sort_values("novel_frac_myeloma")
    shared = list(marrow_tab.ct.values)
    x = np.arange(len(shared))
    w = 0.38
    my_v = (marrow_tab.novel_frac_myeloma * 100).values
    cc_v = (marrow_tab.novel_frac_ccus * 100).values
    ax.bar(x - w / 2, my_v, width=w, color=TEAL, edgecolor="white", label="Myeloma marrow", zorder=2)
    ax.bar(x + w / 2, cc_v, width=w, color=CORAL, edgecolor="white", label="CCUS marrow", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace(" cell", "") for c in shared], rotation=30, ha="right", fontsize=5.8)
    ax.set_ylabel("Novel-isoform molecules (%)")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper left", fontsize=5.5, frameon=False)
    ax.text(0.98, 0.95,
            "Stored cell-type structure · not a discovery claim",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.3, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- b cross-study scatter ----
    ax = fig.add_subplot(gs[0, 1])
    tag(ax, "b"); title(ax, "Cross-study agreement of cell-type novelty ranks")
    xv, yv = m.novel_frac_myeloma.values, m.novel_frac_ccus.values
    rp, _ = pearsonr(xv, yv)
    rs, _ = spearmanr(xv, yv)
    sizes = 48 + 140 * (xv + yv) / 2
    ax.plot([0.15, 0.85], [0.15, 0.85], ls="--", color=GREY, lw=0.75, zorder=1)
    ax.scatter(xv, yv, s=sizes, c=(xv + yv) / 2, cmap=CMAP_NOVEL, vmin=0.2, vmax=0.8,
               edgecolors="white", linewidths=0.7, zorder=3)
    label_offsets = {
        "Dendritic cell": (-22, 10),
        "B cell": (7, 11),
        "Monocyte": (7, -11),
    }
    for _, row in m.iterrows():
        lab = row.ct.replace(" cell", "").replace("Dendritic", "DC")
        offset = label_offsets.get(row.ct, (7, 0))
        ax.annotate(
            lab,
            (row.novel_frac_myeloma, row.novel_frac_ccus),
            xytext=offset,
            textcoords="offset points",
            fontsize=5.6,
            color=SLATE,
            va="center",
        )
    ax.set_xlim(0.12, 0.88); ax.set_ylim(0.12, 0.88)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Novel fraction · myeloma marrow")
    ax.set_ylabel("Novel fraction · CCUS marrow")
    card(ax, 0.05, 0.70, 0.42, 0.24, fc="#FFFFFF", ec=CORAL, lw=0.95)
    ax.text(0.26, 0.88, f"r = {rp:.2f}", transform=ax.transAxes, ha="center",
            fontsize=9.5, fontweight="bold", color=CORAL)
    ax.text(0.26, 0.76, f"ρ = {rs:.2f}", transform=ax.transAxes, ha="center",
            fontsize=7.4, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, color=GREYL, lw=0.4, zorder=0)

    # ---- d junction canonicity as an evidence attribute ----
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "d"); title(ax, "Junction evidence: GT–AG rate vs cross-study recurrence")
    x = j.n_studies.values
    yrate = j.gtag_rate.values * 100
    n_j = j.n_junctions.values.astype(float)
    sizes = 20 + 50 * np.log10(n_j / n_j.min())
    ax.plot(x, yrate, color=TEAL, lw=1.55, zorder=2)
    ax.scatter(x, yrate, s=sizes, c=yrate, cmap=CMAP_HEAT, vmin=45, vmax=100,
               edgecolors="white", linewidths=0.5, zorder=3)
    ax.annotate(f"{int(n_j[0]):,} private", (x[0], yrate[0]),
                textcoords="offset points", xytext=(8, 10), fontsize=5.2, color=SLATE)
    ax.annotate(f"{int(n_j[-1]):,} core", (x[-1], yrate[-1]),
                textcoords="offset points", xytext=(-50, 8), fontsize=5.2, color=SLATE)
    card(ax, 0.40, 0.10, 0.56, 0.16, fc="#FFFFFF", ec=TEAL, lw=0.9)
    ax.text(0.68, 0.18, f"{yrate[0]:.1f}% → {yrate[-1]:.1f}% GT–AG",
            transform=ax.transAxes, ha="center", va="center", fontsize=7.5,
            fontweight="bold", color=TEAL)
    ax.set_xlabel("Independent studies reporting the junction")
    ax.set_ylabel("GT–AG canonical rate (%)")
    ax.set_ylim(45, 102)
    ax.set_xticks(x)
    ax.text(0.02, -0.18, "x = cross-study recurrence in stored junction tables (not depth alone)",
            transform=ax.transAxes, fontsize=5.2, color=SLATE, clip_on=False)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- e reference recovery / annotation coverage ----
    ax = fig.add_subplot(gs[2, 0])
    tag(ax, "e"); title(ax, "Reference transcript recovery across cumulative runs")
    xr, yr = d.n_runs.values, d.cum_frac.values * 100
    ax.fill_between(xr, yr, 100, color=CORALL, alpha=0.35, zorder=1)
    ax.fill_between(xr, 0, yr, color=TEAL, alpha=0.35, zorder=2)
    ax.plot(xr, yr, color=TEAL, lw=1.55, zorder=3)
    ax.set_xlabel("Cumulative reprocessed runs")
    ax.set_ylabel("Reference transcripts seen as FSM (%)")
    ax.set_ylim(0, 100)
    ax.set_xlim(0, xr.max())
    card(ax, 0.48, 0.55, 0.48, 0.35, fc="#FFFFFF", ec=TEAL, lw=0.85)
    ax.text(0.72, 0.82, f"observed {summ['observed_frac']*100:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=6.8, fontweight="bold", color=TEAL)
    ax.text(0.72, 0.70, f"unseen {summ['debt_frac']*100:.1f}%",
            transform=ax.transAxes, ha="center", fontsize=6.8, fontweight="bold", color=CORAL)
    ax.text(0.72, 0.58, f"{summ['n_runs']} runs · {summ['n_studies']} datasets",
            transform=ax.transAxes, ha="center", fontsize=5.2, color=SLATE)
    sys_nov = cls.groupby("system").novel_frac.mean().sort_values(ascending=False)
    ax.text(0.02, 0.95, "Novel fraction by system:", transform=ax.transAxes,
            fontsize=5.3, color=SLATE, va="top")
    ytxt = 0.88
    for sys, nv in sys_nov.items():
        ax.text(0.02, ytxt, f"  {sys}: {nv*100:.0f}%", transform=ax.transAxes,
                fontsize=5.2, color=INK, va="top", fontfamily="sans-serif")
        ytxt -= 0.07
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- f precomputed ASE/DIU/APA inventory ----
    ax = fig.add_subplot(gs[2, 1])
    tag(ax, "f"); title(ax, "Biological-unit-aware ASE / DIU / APA inventory")
    labels = ["Allelic balance\n(ASE)", "Isoform usage\n(DIU)", "poly(A) choice\n(APA)"]
    n_test = np.array([r[1] for r in axes_rows], float)
    n_sig = np.array([r[2] for r in axes_rows], float)
    colors = [GREY, TEAL, GOLD]
    x = np.arange(3)
    w = 0.34
    ax.bar(x - w / 2, n_test / 1000, width=w, color=[c + "55" for c in colors],
           edgecolor="white", label="genes tested", zorder=2)
    ax.bar(x + w / 2, n_sig / 1000, width=w, color=colors,
           edgecolor="white", label="cell-type differential", zorder=2)
    for i, (nt, ns, col) in enumerate(zip(n_test, n_sig, colors)):
        ax.text(i - w / 2, nt / 1000 + 0.25, f"{int(nt):,}", ha="center", fontsize=5.2, color=SLATE)
        ax.text(i + w / 2, ns / 1000 + 0.25, f"{int(ns):,}", ha="center", fontsize=5.6,
                fontweight="bold", color=col)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("Genes (×10³)")
    ax.set_ylim(0, max(n_test.max() / 1000 * 1.22, 12))
    ax.legend(loc="upper left", fontsize=5.3, frameon=False)
    null_txt = "  ·  ".join(
        f"{name} null {null_ranges[name][0]*100:.1f}–{null_ranges[name][1]*100:.1f}%"
        for name in ("ASE", "DIU", "APA")
    )
    # centred on panel f's axes this ran off the right page edge and the last
    # word was cut mid-character; it belongs on the figure footer, left-aligned
    # with the other footer line, where its length cannot overflow
    null_footer = null_txt + "  |  3 full null seeds/axis; 0 FDR discoveries"
    ax.text(0.98, 0.95,
            "2 marrow studies\n25 runs → 19 independent sources",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.3, color=SLATE)
    S.style_ax(ax)
    ax.grid(True, axis="y", color=GREYL, lw=0.4, zorder=0)

    # ---- c molecule composition by system ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "c"); title(ax, "Stored splice-class composition across systems")
    order = ["FSM", "ISM", "NIC", "NNC"]
    systems = list(cls.groupby("system")["novel_frac"].mean().sort_values(ascending=False).index)
    mat, novels = [], []
    for sys in systems:
        sub = cls[cls.system == sys]
        means = np.array([sub[c].mean() for c in order], float)
        means = means / means.sum()
        mat.append(means)
        novels.append(float(sub.novel_frac.mean()))
    mat = np.array(mat)
    y = np.arange(len(systems))
    left = np.zeros(len(systems))
    for i, c in enumerate(order):
        ax.barh(y, mat[:, i], left=left, height=0.66, color=S.CLASS[c],
                edgecolor="white", linewidth=0.3, label=c, zorder=2)
        left += mat[:, i]
    ax.scatter(novels, y, s=46, c=CORAL, marker="D", edgecolors="white",
               linewidths=0.6, zorder=4, label="novel %")
    for yi, nv in zip(y, novels):
        ax.text(
            nv + 0.025,
            yi,
            f"{nv*100:.0f}%",
            va="center",
            fontsize=5.3,
            color=CORAL,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 0.7},
        )
    ax.set_yticks(y)
    ax.set_yticklabels(systems, fontsize=6.0)
    ax.set_xlim(0, 1.42)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("Fraction of spliced molecules")
    ax.axvline(0.5, color=GREYL, lw=0.5, ls=":", zorder=0)
    ax.legend(loc="center left", bbox_to_anchor=(0.73, 0.5), ncol=1,
              fontsize=5.0, handlelength=0.8, labelspacing=0.25, frameon=False)
    S.style_ax(ax)

    fig.suptitle(
        "scTHREAD · Recurrence, negative controls and stored evidence",
        x=0.09, y=0.975, ha="left", fontsize=10, fontweight="bold", color=INK,
    )
    fig.text(
        0.09, 0.028,
        # the permutation inventory is panel f; panel e is the recovery curve
        "Panels a–d: frozen Paper1 evidence tables. Panel f: corrected source-level "
        "restricted permutations across seven shared marrow lineages; FDR<0.05 and effect≥0.20.",
        fontsize=5.2, color=SLATE, ha="left",
    )
    fig.text(0.09, 0.012, null_footer, fontsize=5.2, color=SLATE, ha="left")
    S.save(fig, stem)
    print(f"{stem} bio ok")


# ===================================================================== Fig 3
def render_fig3(stem: str = "NAR_Fig3"):
    """PTPRC utility walkthrough (NAR Fig.3): query → gene card → usage → export.

    stem: output basename (use NAR_Fig3_v2 to avoid clobbering prior binaries).
    """
    GENE, GID = "PTPRC", "ENSG00000081237"
    diu = pd.read_csv(S.F2DATA / "diu_celltype.tsv", sep="\t")
    apa = pd.read_csv(S.F2DATA / "apa_celltype.tsv", sep="\t")
    ase = pd.read_csv(S.F2DATA / "ase_interaction.tsv", sep="\t")
    drow = diu[diu.gene == GID].iloc[0]
    arow = apa[apa.gene == GID].iloc[0]
    srow = ase[ase.gene.astype(str).str.upper() == GENE].iloc[0]
    pt = pd.read_csv(S.FIGDATA / "ptprc_isoform_usage.tsv", sep="\t")
    portal = ASSETS / "ptprc_gene_card.png"
    path_gpt = ASSETS / "Fig3a_query_path_gpt.png"
    # optional API snapshot numbers for export panel
    api_ov = Path(__file__).resolve().parents[2] / "tables" / "PTPRC_api_overview.json"
    api_note = ""
    if api_ov.exists():
        try:
            api_j = json.loads(api_ov.read_text())
            api_note = f"API snapshot keys: {', '.join(list(api_j)[:6])}"
        except Exception:
            api_note = ""

    fig = plt.figure(figsize=(183 * S.MM, 210 * S.MM), facecolor="none")
    gs = fig.add_gridspec(
        3, 2, height_ratios=[0.88, 1.18, 1.08],
        hspace=0.40, wspace=0.26,
        left=0.08, right=0.98, top=0.915, bottom=0.05,
    )

    # ---- a product walkthrough header ----
    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    tag(ax, "a", x=-0.01, y=1.05)
    title(ax, f"Utility case: query {GENE} (CD45) multi-layer evidence in scTHREAD")
    cards = [
        (0.01, "Gene (textbook)",
         f"{GENE} / CD45 has lineage-linked\nisoform programs (RA / RO / RB).\n"
         "Used here as a portal demo gene.",
         TEAL),
        (0.26, "What the database returns",
         "Full-length isoform usage by cell type\n+ junction / PAS / ASE layers\nin one gene card (not gene-only counts).",
         BLUE),
        (0.51, "Precomputed maps (tables)",
         f"DIU sig (q={drow.qval:.3g}, effect={drow.effect:.2f})\n"
         f"APA sig (q={arow.qval:.3g}, effect={arow.effect:.2f})\n"
         f"ASE not CT-diff (q={srow.qval:.2g})",
         GOLD),
        (0.76, "Access path",
         "Search → Browser multi-layer card\n→ isoform matrix / junctions\n→ CSV / JSON / API export",
         CORAL),
    ]
    for x, head, body, col in cards:
        card(ax, x, 0.08, 0.23, 0.78, fc="#FFFFFF", ec=col, lw=1.05, rs=0.04)
        ax.text(x + 0.115, 0.72, head, transform=ax.transAxes, ha="center",
                fontsize=6.5, fontweight="bold", color=col)
        ax.text(x + 0.115, 0.38, body, transform=ax.transAxes, ha="center", va="center",
                fontsize=5.5, color=INK, linespacing=1.3)

    # ---- b portal ----
    ax = fig.add_subplot(gs[1, 0])
    tag(ax, "b"); title(ax, f"Live multi-layer gene card · {GENE}")
    if portal.exists():
        img = mpimg.imread(portal)
        ax.imshow(img, aspect="auto", interpolation="lanczos")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_color(SLATE); sp.set_linewidth(0.7)
    else:
        ax.axis("off"); ax.text(0.5, 0.5, "portal shot missing", ha="center", color=CORAL)
    ax.text(0.5, -0.06, "https://scthread.ai4sc.ac.cn/browse?query=PTPRC",
            transform=ax.transAxes, ha="center", fontsize=5.4, color=SLATE)

    # ---- c isoform usage heatmap from stored table ----
    ax = fig.add_subplot(gs[1, 1])
    tag(ax, "c"); title(ax, f"{GENE} isoform usage by cell type (stored table)")
    top = pt.groupby("transcript_id")["count"].sum().nlargest(5).index.tolist()
    sub = pt[pt.transcript_id.isin(top)].copy()
    name_map = {
        "ENST00000367364": "…367364 (mono/T/plasma-high)",
        "ENST00000697630": "…697630 (B/NK/prog-high)",
        "ENST00000484135": "…484135",
        "ENST00000697635": "…697635",
        "ENST00000442510": "…442510",
    }
    sub["iso"] = sub.transcript_id.map(lambda t: name_map.get(t, "…" + t[-6:]))
    ct_ord = [c for c in ["Progenitor", "B cell", "T cell", "NK", "Monocyte",
                          "Dendritic cell", "Plasma cell", "Erythroid"] if c in set(sub.ct)]
    mat = sub.pivot_table(index="ct", columns="iso", values="frac", aggfunc="sum").reindex(ct_ord).fillna(0)
    if "Monocyte" in mat.index:
        mat = mat[mat.loc["Monocyte"].sort_values(ascending=False).index]
    im = ax.imshow(mat.values, aspect="auto", cmap=CMAP_HEAT, vmin=0, vmax=0.65)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v >= 0.08:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.0,
                        color="white" if v > 0.35 else INK,
                        fontweight="bold" if v > 0.4 else "normal")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=25, ha="right", fontsize=5.0)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=6.0)
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Usage fraction", fontsize=5.4)
    cbar.ax.tick_params(labelsize=4.8)
    ax.text(0.0, -0.28,
            "Source: figdata/ptprc_isoform_usage.tsv · utility demo of multi-layer export,\n"
            "not a mechanism claim. Monocyte/T/plasma vs B/NK/progenitor programs as stored.",
            transform=ax.transAxes, fontsize=5.3, color=SLATE)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(GREYL); sp.set_linewidth(0.5)

    # ---- d query path GPT ----
    ax = fig.add_subplot(gs[2, 0])
    tag(ax, "d"); title(ax, "Reproducible query path (portal, no local reanalysis)")
    if path_gpt.exists():
        show_img(ax, path_gpt)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "query-path schematic missing", ha="center", color=CORAL)

    # ---- e API + biology export contract ----
    ax = fig.add_subplot(gs[2, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    tag(ax, "e"); title(ax, "Export the same multi-layer biology via API")
    card(ax, 0.02, 0.05, 0.96, 0.88, fc="#1C2130", ec="#1C2130", lw=0, rs=0.05)
    for i, col in enumerate(("#FF5F56", "#FFBD2E", "#27C93F")):
        ax.add_patch(Circle((0.08 + i * 0.035, 0.84), 0.012, transform=ax.transAxes,
                            facecolor=col, edgecolor="none", clip_on=False))
    ax.text(0.52, 0.84, "api · scTHREAD · gene biology export",
            transform=ax.transAxes, ha="center", fontsize=5.4, color="#8A8F98")
    ax.text(0.07, 0.72, f"$ curl -s .../api/gene/{GENE}/overview",
            transform=ax.transAxes, ha="left", fontsize=5.2, color="#A7C4C0", family="monospace")
    snippet = (
        "{\n"
        f'  "gene": "{GENE} (CD45)",\n'
        f'  "lineage_isoform_DIU": {{"sig": true, "q": {drow.qval:.4f}}},\n'
        f'  "APA": {{"sig": true, "q": {arow.qval:.4f}}},\n'
        f'  "ASE_CT": {{"sig": false, "q": {srow.qval:.1f}}},\n'
        '  "layers": ["isoform","junction","PAS","allele"],\n'
        '  "export": ["csv","json","api"]\n'
        "}"
    )
    ax.text(0.07, 0.36, snippet, transform=ax.transAxes, ha="left", va="center",
            fontsize=5.5, color="#E8DFD0", family="monospace", linespacing=1.3)
    ax.text(0.07, 0.10,
            "Purpose: export multi-layer biology for a gene — not a new mechanism claim."
            + (f"  {api_note}" if api_note else ""),
            transform=ax.transAxes, ha="left", fontsize=4.8, color="#8A8F98", style="italic")

    fig.suptitle(
        f"scTHREAD · Database utility: {GENE}/CD45 multi-layer query & export",
        x=0.08, y=0.975, ha="left", fontsize=10, fontweight="bold", color=INK,
    )
    fig.text(
        0.08, 0.012,
        "Live URL: https://scthread.ai4sc.ac.cn/browse?query=PTPRC · numbers from f2_grammar tables "
        "+ figdata/ptprc_isoform_usage.tsv · utility demo only",
        fontsize=5.2, color=SLATE, ha="left",
    )
    S.save(fig, stem)
    print(f"{stem} bio ok")


def main(argv: list[str] | None = None):
    import argparse
    p = argparse.ArgumentParser(description="Render scTHREAD NAR biology figures")
    p.add_argument("--fig", choices=["1", "2", "3", "all"], default="all")
    p.add_argument("--stem1", default="NAR_Fig1", help="Fig1 output stem (e.g. NAR_Fig1_v2)")
    p.add_argument("--stem2", default="NAR_Fig2", help="Fig2 output stem (e.g. NAR_Fig2_v2)")
    p.add_argument("--stem3", default="NAR_Fig3", help="Fig3 output stem (e.g. NAR_Fig3_v2)")
    p.add_argument(
        "--observed-suffix",
        default="",
        help="Optional suffix for observed tables, for example _9999",
    )
    args = p.parse_args(argv)
    if args.fig in ("1", "all"):
        render_fig1(stem=args.stem1)
    if args.fig in ("2", "all"):
        render_fig2(stem=args.stem2, observed_suffix=args.observed_suffix)
    if args.fig in ("3", "all"):
        render_fig3(stem=args.stem3)


if __name__ == "__main__":
    main()
