#!/usr/bin/env python3
"""Verify that every figure, panel and supplementary callout in the draft resolves.

A dangling "Supplementary Table S14" or a panel letter that no longer exists after a
figure is redrawn is invisible to prose review but is exactly what a reviewer trips
over. This checks the manuscript's callouts against the artefacts actually present.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
WORKBOOK = ROOT / "manuscript/supplementary/scTHREAD_Supplementary_Tables.xlsx"
FIGDIR = ROOT / "figures"
SF_LEGENDS = ROOT / "docs/NAR_SF_LEGENDS.md"


def _expand(spec: str) -> set[str]:
    """'a', 'b,c' and 'a-f' all become the set of panel letters they name."""
    letters: set[str] = set()
    for part in spec.split(","):
        if re.search(r"[–-]", part):
            a, b = re.split(r"[–-]", part)
            letters.update(chr(c) for c in range(ord(a), ord(b) + 1))
        else:
            letters.add(part)
    return letters


def _declared_panels(path: Path, heading: str) -> dict[int, set[str]]:
    """Panel letters each legend declares, keyed by figure number.

    Legends mark panels as **(a)** or **(b,c)**; a figure with no such marker
    (a screenshot gallery, say) maps to an empty set, which correctly rejects
    any panel-level callout against it.
    """
    if not path.is_file():
        return {}
    text = path.read_text(encoding="utf-8")
    sections = list(re.finditer(heading, text, re.M))
    out: dict[int, set[str]] = {}
    for i, m in enumerate(sections):
        end = sections[i + 1].start() if i + 1 < len(sections) else len(text)
        letters: set[str] = set()
        for pm in re.finditer(r"\*\*\(([a-z](?:,[a-z])*(?:[–-][a-z])?)\)\*\*",
                              text[m.end():end]):
            letters |= _expand(pm.group(1))
        out[int(m.group(1))] = letters
    return out


def main() -> int:
    text = DRAFT.read_text(encoding="utf-8")
    body = text.split("## References", 1)[0]
    problems: list[str] = []

    # --- supplementary tables ---------------------------------------------
    cited_tables = set()
    for m in re.finditer(r"Supplementary Tables? (S\d+)(?:[–-](S?\d+))?", body):
        lo = int(m.group(1)[1:])
        hi = int(re.sub(r"\D", "", m.group(2))) if m.group(2) else lo
        cited_tables.update(f"S{n}" for n in range(lo, hi + 1))
    try:
        import pandas as pd
        sheets = pd.ExcelFile(WORKBOOK).sheet_names
    except Exception as exc:  # noqa: BLE001
        problems.append(f"cannot read supplementary workbook: {exc}")
        sheets = []
    present_tables = {s.split("_", 1)[0] for s in sheets if re.fullmatch(r"S\d+_.*", s)}
    for t in sorted(cited_tables, key=lambda s: int(s[1:])):
        if present_tables and t not in present_tables:
            problems.append(f"cited but missing from workbook: Supplementary Table {t}")
    for t in sorted(present_tables - cited_tables, key=lambda s: int(s[1:])):
        problems.append(f"in workbook but never cited: Supplementary Table {t}")

    # --- supplementary figures --------------------------------------------
    cited_sf = set()
    for m in re.finditer(r"Supplementary Figures? (S\d+)(?:[–-](S?\d+))?", body):
        lo = int(m.group(1)[1:])
        hi = int(re.sub(r"\D", "", m.group(2))) if m.group(2) else lo
        cited_sf.update(range(lo, hi + 1))
    for n in sorted(cited_sf):
        if not (FIGDIR / f"NAR_SF{n}.pdf").is_file():
            problems.append(f"cited but no file: Supplementary Figure S{n} "
                            f"(expected figures/NAR_SF{n}.pdf)")

    # A figure sitting in figures/ that nothing points at is either a leftover
    # from a previous numbering or a callout that was dropped in an edit. Both
    # have happened; the tables above are already checked this way.
    present_sf = {int(p.stem.removeprefix("NAR_SF"))
                  for p in FIGDIR.glob("NAR_SF*.pdf")
                  if re.fullmatch(r"NAR_SF\d+", p.stem)}
    for n in sorted(present_sf - cited_sf):
        problems.append(f"in figures/ but never cited: Supplementary Figure S{n}")

    # NAR numbers supplementary figures in order of first mention. Ranges such
    # as "Supplementary Figures S1-S4" are inventory statements, not narrative
    # first mentions, so they must not set the order.
    first_seen: dict[int, int] = {}
    for m in re.finditer(r"Supplementary Figure (S\d+)(?![–-]S?\d)", body):
        first_seen.setdefault(int(m.group(1)[1:]), m.start())
    order = [n for n, _ in sorted(first_seen.items(), key=lambda kv: kv[1])]
    if order != sorted(order):
        problems.append(
            f"supplementary figures are not numbered in order of first mention: "
            f"first mentions run {order}, expected {sorted(order)}")

    # --- supplementary figure panels --------------------------------------
    sf_panels = _declared_panels(SF_LEGENDS, r"^## Supplementary Figure (\d+)")
    for m in re.finditer(r"Supplementary Figure S(\d+)([a-z](?:,[a-z])*(?:[–-][a-z])?)\b",
                         body):
        n, spec = int(m.group(1)), m.group(2)
        if n not in sf_panels:
            problems.append(f"Supplementary Figure S{n} has no legend in "
                            f"{SF_LEGENDS.name}, so panel '{spec}' cannot be checked")
            continue
        missing = _expand(spec) - sf_panels[n]
        if missing:
            problems.append(f"Supplementary Figure S{n} callout references panel(s) "
                            f"{sorted(missing)} not defined in its legend")

    # --- main figure panels -----------------------------------------------
    legends = {}
    for m in re.finditer(r"^### Figure (\d)\.(.*?)(?=^### |\Z)", text, re.M | re.S):
        legends[m.group(1)] = m.group(2)
    # Legends write ranges with either an en dash or a plain hyphen; the callout
    # side already accepts both, so the declaration side must too.
    declared = {fig: set(re.findall(r"\(([a-z](?:,[a-z])*(?:[–-][a-z])?)\)", body_txt))
                for fig, body_txt in legends.items()}
    panel_letters = {}
    for fig, groups in declared.items():
        letters = set()
        for g in groups:
            for part in g.split(","):
                if re.search(r"[–-]", part):
                    a, b = re.split(r"[–-]", part)
                    letters.update(chr(c) for c in range(ord(a), ord(b) + 1))
                else:
                    letters.add(part)
        panel_letters[fig] = letters
    for m in re.finditer(r"Figure (\d)([a-z](?:,[a-z])*(?:[–-][a-z])?)", body):
        fig, spec = m.group(1), m.group(2)
        wanted = set()
        for part in re.split(r",", spec):
            if re.search(r"[–-]", part):
                a, b = re.split(r"[–-]", part)
                wanted.update(chr(c) for c in range(ord(a), ord(b) + 1))
            else:
                wanted.add(part)
        missing = wanted - panel_letters.get(fig, set())
        if missing:
            problems.append(f"Figure {fig} callout references panel(s) "
                            f"{sorted(missing)} not defined in its legend")

    print(f"cited supplementary tables : {len(cited_tables)}")
    print(f"cited supplementary figures: {len(cited_sf)}")
    print(f"main figures with legends  : {sorted(legends)}")
    for fig in sorted(panel_letters):
        print(f"  Figure {fig} panels in legend: {''.join(sorted(panel_letters[fig]))}")
    if problems:
        print("\nPROBLEMS")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nALL CROSS-REFERENCES RESOLVE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
