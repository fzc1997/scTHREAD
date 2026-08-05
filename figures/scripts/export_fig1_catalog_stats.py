#!/usr/bin/env python3
"""Export Fig.1 catalog statistics from sample_registry + figdata (no invented numbers).

Writes TSVs under NAR_database/tables/fig1_*.tsv and a dated KPI freeze.
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_nar_bio as R  # noqa: E402
import nar_style as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "tables"
REG = R.REG
OUT_PREFIX = TABLES / "fig1"


def _mode_or_blank(s: pd.Series) -> str:
    s = s.dropna().astype(str)
    s = s[s.str.len() > 0]
    return s.mode().iloc[0] if len(s) else ""


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    today = date.today().isoformat().replace("-", "")

    # --- headline ---
    cat = R.catalog_headline()
    ch = json.loads((S.FIGDATA / "cohort_headline.json").read_text())
    headline = pd.DataFrame(
        [
            {"metric": "catalog_isoquant_done_runs", "value": cat["n_runs"], "source": R.catalog_scope_note()},
            {"metric": "catalog_gse", "value": cat["n_studies"], "source": "unique gse among done"},
            {"metric": "catalog_isoquant_cells", "value": cat["n_cells"], "source": "sum isoquant_cells"},
            {"metric": "catalog_human_runs", "value": cat["n_human"], "source": "species==human"},
            {"metric": "catalog_mouse_runs", "value": cat["n_mouse"], "source": "species==mouse"},
            {"metric": "catalog_ont_runs", "value": cat["n_ont"], "source": "platform==ONT"},
            {"metric": "catalog_pacbio_runs", "value": cat["n_pacbio"], "source": "platform==PacBio"},
            {"metric": "labeled_studies", "value": ch["n_studies"], "source": "cohort_headline.json"},
            {"metric": "labeled_runs", "value": ch["n_runs"], "source": "cohort_headline.json"},
            {"metric": "labeled_annotated_cells", "value": ch["n_annotated_cells"], "source": "cohort_headline.json"},
            {"metric": "isoforms_catalog_bound", "value": ">200k", "source": "qualitative bound (portal/docs); not recounted this freeze"},
        ]
    )
    headline.to_csv(f"{OUT_PREFIX}_catalog_headline.tsv", sep="\t", index=False)

    # --- dual scope ---
    dual = pd.DataFrame(
        [
            {
                "scope": "catalog",
                "definition": f"IsoQuant-complete runs, {R.catalog_scope_note()}",
                "n_studies": cat["n_studies"],
                "n_runs": cat["n_runs"],
                "n_cells": cat["n_cells"],
                "cell_definition": "isoquant_cells",
            },
            {
                "scope": "labeled_analysis_subset",
                "definition": "donor x cell-type labels for DIU/APA/ASE tables",
                "n_studies": ch["n_studies"],
                "n_runs": ch["n_runs"],
                "n_cells": ch["n_annotated_cells"],
                "cell_definition": "annotated cells (cohort_headline)",
            },
        ]
    )
    dual.to_csv(f"{OUT_PREFIX}_dual_scope.tsv", sep="\t", index=False)

    # --- system composition ---
    studies = R.load_catalog_bio()
    sys_disp = {}  # biological names as-is; no Methods/Benchmark display alias
    sys_rows = []
    for system, sub in studies.groupby("system"):
        # never export Benchmark/Methods as a figure category
        if system in ("Benchmark", "Methods"):
            continue
        sys_rows.append(
            {
                "system": system,
                "system_display": sys_disp.get(system, system),
                "n_gse": int(sub.gse.nunique()),
                "n_runs": int(sub.n_runs.sum()),
                "isoquant_cells": int(sub.n_cells.sum()),
            }
        )
    sys_df = pd.DataFrame(sys_rows).sort_values("isoquant_cells", ascending=False)
    sys_df.to_csv(f"{OUT_PREFIX}_system_composition.tsv", sep="\t", index=False)

    # --- species × platform ---
    done = R.registry_done()
    done["isoquant_cells"] = pd.to_numeric(done.get("isoquant_cells"), errors="coerce").fillna(0)
    done["species_bin"] = done["species"].fillna("?").str.lower().map(
        lambda s: "mouse" if "mouse" in s else ("human" if "human" in s else "other")
    )
    done["platform_bin"] = done["platform"].fillna("?").astype(str)
    done.loc[~done["platform_bin"].isin(["ONT", "PacBio"]), "platform_bin"] = "other"

    sp_plat = (
        done.groupby(["species_bin", "platform_bin"], dropna=False)
        .agg(n_runs=("srr", "nunique"), isoquant_cells=("isoquant_cells", "sum"))
        .reset_index()
    )
    sp_plat["isoquant_cells"] = sp_plat["isoquant_cells"].astype(int)
    sp_plat.to_csv(f"{OUT_PREFIX}_species_platform.tsv", sep="\t", index=False)

    # --- study atlas ---
    atlas = studies.copy()
    atlas["system_display"] = atlas["system"].map(lambda s: sys_disp.get(s, s))
    keep = [
        c
        for c in [
            "gse",
            "tissue",
            "system",
            "system_display",
            "species",
            "platform",
            "n_runs",
            "n_cells",
            "biology_group",
            "description",
            "ASE",
            "DIU",
            "APA",
            "F2jct",
            "labeled_frac_all",
        ]
        if c in atlas.columns
    ]
    atlas[keep].sort_values("n_cells", ascending=False).to_csv(
        f"{OUT_PREFIX}_study_atlas.tsv", sep="\t", index=False
    )

    # --- rolling-view atlas (datasets served by the portal but outside the snapshot) ---
    # The manuscript snapshot is frozen and RNA-only; datasets such as the mouse
    # gastrulation cohort (OWN_ASE_scONT) are served by the rolling web view and
    # must be reported from a separate table so they never enter snapshot totals.
    live = pd.read_csv(REG, sep="\t", dtype=str)
    live = live[live["isoquant_status"].fillna("").str.lower().eq("done")].copy()
    live["isoquant_cells"] = pd.to_numeric(live["isoquant_cells"], errors="coerce").fillna(0)
    rolling = live[~live["srr"].isin(set(R.registry_done()["srr"]))]
    rolling_atlas = (
        rolling.groupby("gse", as_index=False)
        .agg(n_runs=("srr", "nunique"),
             n_cells=("isoquant_cells", "sum"),
             species=("species", _mode_or_blank),
             platform=("platform", _mode_or_blank))
        .sort_values("n_cells", ascending=False)
    )
    rolling_atlas["n_cells"] = rolling_atlas["n_cells"].astype(int)
    rolling_atlas["scope"] = "rolling_web_view_not_in_manuscript_snapshot"
    rolling_atlas.to_csv(f"{OUT_PREFIX}_rolling_view_atlas.tsv", sep="\t", index=False)

    # --- layer coverage (labeled GSE only) ---
    cov_path = S.FIGDATA / "ed7_coverage_matrix.tsv"
    if cov_path.exists():
        cov = pd.read_csv(cov_path, sep="\t")
        cov.to_csv(f"{OUT_PREFIX}_layer_coverage_labeled.tsv", sep="\t", index=False)
        # summary
        layer_cols = [c for c in cov.columns if c != "gse"]
        summ = []
        for c in layer_cols:
            v = pd.to_numeric(cov[c], errors="coerce").fillna(0)
            summ.append(
                {
                    "layer": c,
                    "n_gse_with_any": int((v > 0).sum()),
                    "n_gse_total": int(len(cov)),
                    "sum_run_counts": int(v.sum()),
                    "source": str(cov_path),
                }
            )
        pd.DataFrame(summ).to_csv(f"{OUT_PREFIX}_layer_coverage_summary.tsv", sep="\t", index=False)

    # --- optional build funnel from registry status columns ---
    # Upstream acquisition status is a property of the whole registry, not of the
    # frozen snapshot, so this block is explicitly labelled as live.
    reg = pd.read_csv(REG, sep="\t", dtype=str)
    funnel_rows = []
    if "status" in reg.columns:
        st = reg["status"].fillna("").str.lower().value_counts()
        for k, v in st.items():
            funnel_rows.append({"stage": f"registry_status:{k or 'blank'}", "n": int(v),
                                "note": f"all live registry rows ({today})"})
    funnel_rows.append({"stage": "isoquant_done_runs", "n": cat["n_runs"],
                        "note": R.catalog_scope_note()})
    funnel_rows.append({"stage": "labeled_runs", "n": ch["n_runs"], "note": "cohort_headline analysis subset"})
    pd.DataFrame(funnel_rows).to_csv(f"{OUT_PREFIX}_build_funnel.tsv", sep="\t", index=False)

    # --- dated freeze ---
    freeze = TABLES / f"Table_S_frozen_kpis_{today}.tsv"
    freeze_df = pd.DataFrame(
        [
            {"kpi": "catalog_isoquant_done_runs", "value": cat["n_runs"], "note": R.catalog_scope_note()},
            {"kpi": "catalog_gse", "value": cat["n_studies"], "note": "unique GSE among done"},
            {"kpi": "catalog_isoquant_cells", "value": cat["n_cells"], "note": "sum isoquant_cells"},
            {"kpi": "catalog_human_runs", "value": cat["n_human"], "note": ""},
            {"kpi": "catalog_mouse_runs", "value": cat["n_mouse"], "note": ""},
            {"kpi": "catalog_ont_runs", "value": cat["n_ont"], "note": ""},
            {"kpi": "catalog_pacbio_runs", "value": cat["n_pacbio"], "note": ""},
            {"kpi": "analysis_studies", "value": ch["n_studies"], "note": "cohort_headline.json"},
            {"kpi": "analysis_runs", "value": ch["n_runs"], "note": "cohort_headline.json"},
            {"kpi": "analysis_annotated_cells", "value": ch["n_annotated_cells"], "note": "cohort_headline.json"},
            {"kpi": "isoforms_catalog", "value": ">200k", "note": "qualitative bound; not re-counted this freeze"},
        ]
    )
    freeze_df.to_csv(freeze, sep="\t", index=False)

    # supersede pointer
    for old in TABLES.glob("Table_S_frozen_kpis_*.tsv"):
        if old.name == freeze.name:
            continue
        note = old.with_name(old.stem + "_SUPERSEDED.txt")
        if not note.exists() and old.name.startswith("Table_S_frozen_kpis_2026"):
            note.write_text(f"SUPERSEDED by {freeze.name}\n")

    manifest = {
        "date": today,
        "catalog": cat,
        "labeled": ch,
        "outputs": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("fig1_*.tsv"))
        + [str(freeze.relative_to(ROOT))],
    }
    (TABLES / "fig1_export_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
