#!/usr/bin/env python3
"""Derive the tissue inventory that backs the graphical-abstract icon strip.

Joins three sources that are already the authority elsewhere in the project:

* ``metadata/annotation/authoritative_cellcount_clean_20260802.tsv`` -- the 34
  in-scope studies and their authoritative study-level cell counts (453 runs,
  923,389 cells).  Cell counts are study-level and must not be summed per run.
* ``docs/sample_registry.tsv`` -- run-level GEO ``source_name`` and
  ``characteristics`` used to read off the biological source.
* ``figures/scripts/prepare_ga_v9_sources.py`` -- ``SYSTEM_MAP``, the frozen
  assignment of each study to one of six mutually exclusive release classes.

Output: ``tables/ga_tissue_icon_map.tsv``, one row per study, with the icon that
represents it in the strip.  Nothing here is inferred beyond the registry text;
studies whose source is a cultured line are labelled as such rather than being
attributed to the organ the line was derived from.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import os
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCLONG = Path(os.environ.get("SCTHREAD_SCLONG_ROOT", "/gpfs/home/fuzc/Seq_Database/scLong"))
AUTHORITATIVE = SCLONG / "metadata/annotation/authoritative_cellcount_clean_20260802.tsv"
REGISTRY = SCLONG / "docs/sample_registry.tsv"
V9_SOURCES = ROOT / "figures/scripts/prepare_ga_v9_sources.py"

# study -> (biological source as written in the registry, icons it backs).
# A study may back more than one icon (e.g. a cohort that mixes primary tumour
# nuclei with cancer cell lines), so the per-icon cell totals printed at the end
# double-count those studies and do not sum to 923,389.  The per-study rows do.
#
# Recorded explicitly rather than parsed, because the registry free text is
# ambiguous in exactly the places that matter: GSE212945 lists "Skin" for runs
# that are melanoma cell lines (A375/MEL1/MEL2), GSE283658 lists "neonatal
# foreskin" for the BJ fibroblast line, and GSE303762 lists lung-cancer cell
# lines rather than lung tissue.  None of those are primary skin or lung.
STUDY_TISSUE = {
    "GSE178175": ("Frontal cortex; hippocampus", ["Brain"]),
    "GSE314176": ("Hippocampus", ["Brain"]),
    "GSE158450": ("Prefrontal cortex; hippocampus", ["Brain"]),
    "GSE130708": ("E18 brain", ["Brain"]),
    "GSE76026": ("Oligodendrocyte lineage", ["Brain"]),
    "GSE283629": ("iPSC-derived neurons and neural progenitors", ["Brain"]),
    "GSE274249": ("iPSC-derived neurons and neural progenitors", ["Brain"]),
    "GSE301658": ("Brain tumour", ["Brain"]),
    "GSE255520": ("Retina", ["Retina"]),
    "GSE114157": ("Cochlea, outer hair cells", ["Cochlea"]),
    "GSE276974": ("Bone marrow, CCUS", ["Marrow"]),
    "GSE307660": ("Bone marrow, myeloma", ["Marrow"]),
    "GSE292324": ("Bone marrow, AML", ["Marrow"]),
    "GSE214231": ("Donor TCR-negative T cells, CAR library", ["Blood"]),
    "GSE295352": ("Spleen", ["Spleen"]),
    "GSE252416": ("Lymph node", ["Lymph"]),
    "GSE252344": ("Lymph node", ["Lymph"]),
    "GSE288222": ("Heart, left ventricle", ["Heart"]),
    "GSE224045": ("Lung adenocarcinoma", ["Lung"]),
    "GSE303762": ("Lung-cancer cell lines", ["Lung", "Cell line"]),
    "GSE212945": ("Renal-cell-carcinoma nuclei; melanoma and lung lines", ["Kidney", "Cell line"]),
    "GSE289790": ("Prostate", ["Prostate"]),
    "GSE248118": ("Ovarian-cancer cell lines", ["Cell line"]),
    "GSE289428": ("HeLa-S3", ["Cell line"]),
    "GSE295932": ("Jurkat T-cell line, CRISPR", ["Cell line"]),
    "GSE295353": ("Pancreatic islet", ["Islet"]),
    "GSA_mouse_testis_11week": ("Testis", ["Testis"]),
    "benagen_human_ovary": ("Ovary", ["Ovary"]),
    "OWN_ASE_scONT": ("Gastrulating embryo, E6.5-E8.5", ["Embryo"]),
    "GSE250381": ("Preimplantation embryo", ["Embryo"]),
    "GSE309071": ("Embryo yolk sac and AGM region", ["Embryo"]),
    "GSE185554": ("Embryo AGM region", ["Embryo"]),
    "GSE274527": ("Embryo yolk sac and AGM region", ["Embryo"]),
    "GSE283658": ("Bone marrow; BJ foreskin fibroblast; reprogrammed HSC", ["Marrow", "Cell line"]),
}


def read_tsv(path: Path) -> list[dict]:
    raw = path.read_text(encoding="utf-8", errors="replace").replace("\r", "")
    return list(csv.DictReader(raw.split("\n"), delimiter="\t"))


def load_system_map() -> dict[str, str]:
    spec = importlib.util.spec_from_file_location("ga_v9_sources", V9_SOURCES)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return dict(module.SYSTEM_MAP)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=ROOT / "tables/ga_tissue_icon_map.tsv")
    parser.add_argument("--show-sources", action="store_true")
    args = parser.parse_args()

    studies = [r for r in read_tsv(AUTHORITATIVE) if r.get("gse")]
    registry = [r for r in read_tsv(REGISTRY) if r.get("srr")]
    system_map = load_system_map()

    by_study = defaultdict(list)
    for row in registry:
        by_study[row["gse"]].append(row)

    missing = [s["gse"] for s in studies if s["gse"] not in STUDY_TISSUE]
    if missing:
        raise SystemExit(f"studies without a recorded tissue: {missing}")

    rows = []
    for study in sorted(studies, key=lambda s: -int(s["cell_count"])):
        gse = study["gse"]
        runs = by_study.get(gse, [])
        sources = Counter(
            (r.get("geo_source_name") or "").strip() for r in runs if (r.get("geo_source_name") or "").strip()
        )
        species = sorted({(r.get("species") or "").strip() for r in runs if (r.get("species") or "").strip()})
        tissue, icons = STUDY_TISSUE[gse]
        rows.append(
            {
                "study": gse,
                "cells": int(study["cell_count"]),
                "runs": len(runs),
                "species": "/".join(species) or "n/a",
                "system": system_map.get(gse, "n/a"),
                "tissue": tissue,
                "icon_labels": ";".join(icons),
                "registry_source_name": "; ".join(k for k, _ in sources.most_common(3)) or "n/a",
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    total = sum(r["cells"] for r in rows)
    print(f"{len(rows)} studies, {total:,} cells -> {args.out}")
    if total != 923_389:
        raise SystemExit(f"cell total {total} does not match the frozen 923,389")

    per_icon = defaultdict(lambda: [0, 0])
    for row in rows:
        for icon in row["icon_labels"].split(";"):
            per_icon[icon][0] += row["cells"]
            per_icon[icon][1] += 1
    print("\nicon\tcells\tstudies   (multi-tissue studies counted under each icon)")
    for icon, (cells, n) in sorted(per_icon.items(), key=lambda kv: -kv[1][0]):
        print(f"{icon}\t{cells:,}\t{n}")

    if args.show_sources:
        print("\nregistry source_name per study:")
        for row in rows:
            print(f"  {row['study']:28s} {row['registry_source_name']}")


if __name__ == "__main__":
    main()
