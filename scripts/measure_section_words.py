#!/usr/bin/env python3
"""Compare our per-section word counts with the published NAR Database Issue corpus.

Corpus counts come from `pdftotext` in reading order with running heads, download
banners and the reference list stripped. Two-column typesetting means a little
figure-caption text bleeds into the nearest section, so corpus numbers are upper
bounds (roughly 5-10% high); our manuscript is measured exactly from Markdown.
"""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "docs/NAR_Database_Issue_2026"
DRAFT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"

# canonical section -> heading spellings seen in the corpus
SECTIONS = {
    "Introduction": [r"introduction"],
    "Materials and methods": [r"materials and methods", r"methods", r"material and methods"],
    "Results": [r"results", r"results and discussion"],
    "Database content": [r"database content", r"database content and usage"],
    "Discussion": [
        r"discussion", r"discussion and conclusions?", r"discussion and perspectives",
        r"conclusions?", r"conclusions? and future extensions",
        r"conclusions? and future directions", r"conclusion and future perspective",
        r"summary and future directions", r"summary and future perspectives",
        r"summary and perspectives",
    ],
}
END = [r"acknowledgements?", r"supplementary data", r"data availability",
       r"funding", r"conflict of interest", r"author contributions", r"references"]

NOISE = re.compile(
    r"^\s*(downloaded from|nucleic acids research|©|©|https?://\S+$|"
    r"d\d{2,4}(\s*[-–]\s*d?\d{2,4})?\s*$|\d+\s*$|"
    r"received:|revised:|accepted:|advance access|database issue)", re.I)


def words(text: str) -> int:
    return len(re.findall(r"\b[\w′'’-]+\b", text))


def clean(lines: list[str]) -> list[str]:
    return [l for l in lines if not NOISE.match(l.strip())]


def heading_index(lines: list[str], patterns: list[str], start: int = 0) -> int | None:
    for i in range(start, len(lines)):
        s = lines[i].strip().rstrip(":").strip()
        if 0 < len(s) <= 60 and any(re.fullmatch(p, s, re.I) for p in patterns):
            return i
    return None


def unwrap_columns(pdf: Path) -> str:
    """Re-serialise a two-column page as left column then right column.

    Plain reading-order extraction interleaves the columns badly enough that
    section headings end up adjacent to one another, which silently drops whole
    sections. `-layout` keeps horizontal position, so the page can be split at
    the gutter and each column read in turn.
    """
    raw = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                         capture_output=True, text=True).stdout
    pages = []
    for page in raw.split("\f"):
        rows = page.split("\n")
        if not rows:
            continue
        width = max((len(r) for r in rows), default=0)
        if width < 60:                       # single-column page
            pages.append(page)
            continue
        # the gutter is the column with the most blank characters near the middle
        lo, hi = int(width * 0.35), int(width * 0.65)
        blanks = [sum(1 for r in rows if len(r) <= c or r[c] == " ")
                  for c in range(lo, hi)]
        gutter = lo + max(range(len(blanks)), key=lambda i: blanks[i]) if blanks else width // 2
        left = "\n".join(r[:gutter].rstrip() for r in rows)
        right = "\n".join(r[gutter:].rstrip() for r in rows)
        pages.append(left + "\n" + right)
    return "\n".join(pages)


def _measure(raw: str, pdf: Path) -> dict:
    lines = clean(raw.splitlines())

    marks: list[tuple[int, str]] = []
    for name, pats in SECTIONS.items():
        idx = heading_index(lines, pats)
        if idx is not None:
            marks.append((idx, name))
    tail = None
    for pat in END:
        idx = heading_index(lines, [pat], start=max((m[0] for m in marks), default=0) + 1)
        if idx is not None:
            tail = idx if tail is None else min(tail, idx)
    marks.sort()

    out: dict[str, int] = {}
    for pos, (idx, name) in enumerate(marks):
        stop = marks[pos + 1][0] if pos + 1 < len(marks) else (tail or len(lines))
        out[name] = words("\n".join(lines[idx + 1:stop]))

    abs_start = heading_index(lines, [r"abstract"])
    if abs_start is not None:
        abs_stop = heading_index(lines, [r"graphical abstract", r"introduction"], abs_start + 1)
        out["Abstract"] = words("\n".join(lines[abs_start + 1:abs_stop or abs_start + 30]))
    flat = subprocess.run(["pdftotext", str(pdf), "-"],
                          capture_output=True, text=True).stdout
    out["References"] = len(re.findall(r"^\s*\d{1,3}\.\s+\S", flat, re.M))
    out["Body"] = sum(v for k, v in out.items() if k not in {"Abstract", "References"})
    return out


def measure_pdf(pdf: Path) -> dict:
    """Measure twice and keep the reading that recovered more of the paper.

    Neither extraction is reliable on every PDF: reading order sometimes emits
    two headings back to back and drops the section between them, while the
    gutter split fails on pages with spanning figures. Taking the variant that
    resolves more sections avoids silently reporting a truncated section.
    """
    flat = subprocess.run(["pdftotext", str(pdf), "-"],
                          capture_output=True, text=True).stdout
    candidates = [_measure(flat, pdf), _measure(unwrap_columns(pdf), pdf)]
    body_sections = ("Introduction", "Materials and methods", "Results",
                     "Database content", "Discussion")
    return max(candidates,
               key=lambda m: (sum(1 for s in body_sections if m.get(s)), m.get("Body", 0)))


def measure_draft() -> dict:
    text = DRAFT.read_text(encoding="utf-8")
    body = text.split("## References", 1)[0]
    blocks = re.split(r"^## ", body, flags=re.M)[1:]
    out: dict[str, int] = {}
    for block in blocks:
        head, _, rest = block.partition("\n")
        name = head.strip()
        # subsection text belongs to its parent section
        out[name] = words(rest)
    out["References"] = len(re.findall(r"^\d{1,3}\.\s+\S",
                                       text.split("## References", 1)[1], re.M))
    keep = ["Introduction", "Materials and methods", "Results", "Discussion"]
    out["Body"] = sum(out.get(k, 0) for k in keep)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "tables/section_word_counts.tsv")
    args = parser.parse_args()

    cols = ["Abstract", "Introduction", "Materials and methods", "Results",
            "Database content", "Discussion", "Body", "References"]
    rows = []
    for pdf in sorted(CORPUS.glob("*.pdf")):
        m = measure_pdf(pdf)
        rows.append({"paper": pdf.stem, **{c: m.get(c, "") for c in cols}})
    draft = measure_draft()
    rows.append({"paper": "scTHREAD (this manuscript)",
                 **{c: draft.get(c, "") for c in cols}})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        fh.write("\t".join(["paper", *cols]) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in ["paper", *cols]) + "\n")

    width = max(len(r["paper"]) for r in rows)
    header = f"{'paper':<{width}}  " + "  ".join(f"{c[:11]:>11}" for c in cols)
    print(header); print("-" * len(header))
    for r in rows:
        print(f"{r['paper']:<{width}}  " + "  ".join(f"{str(r.get(c,'') or '-'):>11}" for c in cols))
    print(f"\nwritten: {args.out}")


if __name__ == "__main__":
    main()
