#!/usr/bin/env python3
"""Check that every panel is cited, and that citations follow panel order.

check_crossrefs.py answers "does every callout point at a panel that exists".
This answers the two questions that miss: does every panel that exists get
called out at all, and does the text walk the panels in the order the figure
draws them. A panel the text never mentions is a panel the reader has no reason
to look at; a text that jumps c -> a -> b makes the reader hunt.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"


def expand(spec: str) -> list[str]:
    out: list[str] = []
    for part in spec.split(","):
        part = part.strip()
        if re.search(r"[-–—]", part):
            a, b = re.split(r"[-–—]", part)[:2]
            out.extend(chr(c) for c in range(ord(a.strip()), ord(b.strip()) + 1))
        elif part:
            out.append(part)
    return out


def legend_panels(text: str, heading: str) -> list[str]:
    """Panels a figure's legend defines, in the order the legend introduces them."""
    m = re.search(rf"^### {re.escape(heading)}(.*?)(?=^### |^## |\Z)", text, re.M | re.S)
    if not m:
        return []
    seen: list[str] = []
    # "poly(A)" is not a panel marker
    for spec in re.findall(r"(?<!poly)\(([A-Za-z](?:\s*,\s*[A-Za-z])*(?:\s*[-–—]\s*[A-Za-z])?)\)", m.group(1)):
        for letter in expand(spec):
            if letter not in seen:
                seen.append(letter)
    return seen


def main() -> int:
    text = DRAFT.read_text(encoding="utf-8")
    body = text.split("\n## References", 1)[0]
    # strip legends: a legend mentioning its own panels is not a text callout
    body_only = re.sub(r"^### (Figure|Supplementary Figure)\b.*?(?=^### |^## |\Z)",
                       "", body, flags=re.M | re.S)

    problems: list[str] = []
    print(f"{'figure':<12} {'panels in legend':<22} {'cited in text':<22} status")
    print("-" * 78)

    # ---------------- main figures ----------------
    fig_headings = dict(re.findall(r"^### (Figure (\d)\..*)$", text, re.M))
    for heading, number in sorted(((h, n) for h, n in fig_headings.items()), key=lambda x: x[1]):
        defined = legend_panels(text, heading)
        cited: list[str] = []
        for m in re.finditer(rf"Figure\s+{number}([A-Za-z](?:\s*,\s*[A-Za-z])*(?:\s*[-–—]\s*[A-Za-z])?)", body_only):
            for letter in expand(m.group(1)):
                if letter not in cited:
                    cited.append(letter)
        missing = [p for p in defined if p not in cited]
        extra = [p for p in cited if p not in defined]
        ordered = cited == sorted(cited, key=lambda c: defined.index(c) if c in defined else 99)
        status = []
        if missing:
            status.append(f"NEVER CITED: {','.join(missing)}")
            problems.append(f"Figure {number}: panels never cited in the text: {','.join(missing)}")
        if extra:
            status.append(f"NOT IN LEGEND: {','.join(extra)}")
            problems.append(f"Figure {number}: cited panels absent from the legend: {','.join(extra)}")
        if cited and not ordered:
            status.append(f"note: cited {''.join(cited)}, drawn {''.join(defined)}")
        print(f"Figure {number:<5} {''.join(defined) or '-':<22} {''.join(cited) or '-':<22} "
              f"{'; '.join(status) or 'ok'}")

    # ---------------- supplementary figures ----------------
    print()
    sf_cited: dict[int, list[str]] = {}
    for m in re.finditer(r"Supplementary\s+Figures?\s+S(\d+)([a-z](?:\s*,\s*[a-z])*(?:\s*[-–—]\s*[a-z])?)?", body_only):
        n = int(m.group(1))
        letters = expand(m.group(2)) if m.group(2) else []
        cur = sf_cited.setdefault(n, [])
        for letter in letters:
            if letter not in cur:
                cur.append(letter)
    legend_file = ROOT / "docs/NAR_SF_LEGENDS.md"
    sf_legends = legend_file.read_text(encoding="utf-8") if legend_file.is_file() else ""
    # legend headings read "## Supplementary Figure 1 | ..."; the text cites
    # "Supplementary Figure S1". Anchoring on the heading keeps the two apart,
    # and keeps the superseded-pack note in the preamble out of the tally.
    heads = [(int(m.group(1)), m.start())
             for m in re.finditer(r"^## Supplementary Figure (\d+)\b", sf_legends, re.M)]
    bounds = {n: (start, heads[i + 1][1] if i + 1 < len(heads) else len(sf_legends))
              for i, (n, start) in enumerate(heads)}
    if sf_cited and not bounds:
        problems.append(
            f"{legend_file.name} defines no '## Supplementary Figure N' headings, "
            "so no supplementary panel can be checked"
        )
    for n in sorted(set(sf_cited) | set(bounds)):
        pdf = ROOT / f"figures/NAR_SF{n}.pdf"
        defined: list[str] = []
        if n in bounds:
            start, end = bounds[n]
            for spec in re.findall(r"\*\*\(([a-z])\)\*\*|\(([a-z](?:\s*,\s*[a-z])*(?:\s*[-–—]\s*[a-z])?)\)",
                                   sf_legends[start:end]):
                for letter in expand(spec[0] or spec[1]):
                    if letter not in defined:
                        defined.append(letter)
        cited = sf_cited.get(n, [])
        status = []
        if not pdf.is_file():
            status.append("NO FILE")
            problems.append(f"Supplementary Figure S{n} has no figures/NAR_SF{n}.pdf")
        if n not in sf_cited:
            status.append("FIGURE NEVER CITED")
            problems.append(f"Supplementary Figure S{n} is never cited in the text")
        elif defined:
            missing = [p for p in defined if p not in cited]
            if missing:
                status.append(f"NEVER CITED: {','.join(missing)}")
                problems.append(
                    f"Supplementary Figure S{n}: panels never cited in the text: {','.join(missing)}"
                )
            if cited != sorted(cited, key=lambda c: defined.index(c) if c in defined else 99):
                status.append(f"note: cited {''.join(cited)}, drawn {''.join(defined)}")
        print(f"SF{n:<10} {''.join(defined) or '-':<22} {''.join(cited) or '(whole figure)':<22} "
              f"{'; '.join(status) or 'ok'}")

    # ---------------- figure order of first mention ----------------
    print()
    order = []
    for m in re.finditer(r"Figure\s+(\d)", body_only):
        if m.group(1) not in order:
            order.append(m.group(1))
    if order != sorted(order):
        problems.append(f"main figures are first mentioned in the order {','.join(order)}")
        print(f"main-figure first-mention order: {','.join(order)}  <-- NOT ASCENDING")
    else:
        print(f"main-figure first-mention order: {','.join(order)}  ok")
    sf_order = []
    for m in re.finditer(r"Supplementary\s+Figures?\s+S(\d+)", body_only):
        if m.group(1) not in sf_order:
            sf_order.append(m.group(1))
    ints = [int(x) for x in sf_order]
    if ints != sorted(ints):
        problems.append(f"supplementary figures are first mentioned in the order {','.join(sf_order)}")
        print(f"supplementary first-mention order: {','.join(sf_order)}  <-- NOT ASCENDING")
    else:
        print(f"supplementary first-mention order: {','.join(sf_order)}  ok")

    if problems:
        print("\nPROBLEMS")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("\nEVERY PANEL IS CITED, IN FIGURE ORDER")
    return 0


if __name__ == "__main__":
    sys.exit(main())
