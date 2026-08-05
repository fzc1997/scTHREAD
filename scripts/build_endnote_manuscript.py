#!/usr/bin/env python3
"""Build an EndNote-ready manuscript: temporary citations plus a matching library.

The submitted manuscript carries plain numeric citations - [1], [2,3], [4-6] -
and its own numbered reference list. That is correct for submission but cannot
be reformatted: changing the reference style, or adding one reference, means
renumbering every citation by hand.

This writes a second pair of files that EndNote can drive:

  manuscript/endnote/scTHREAD_NAR_Database_Issue_endnote.md
      the same text with every [n] replaced by {Surname, Year #n}, and the
      numbered reference list removed, since EndNote regenerates it
  manuscript/endnote/scTHREAD_NAR_references.xml
      an EndNote XML library whose rec-number IS the manuscript's citation
      number, so #n resolves to the intended work

The record number is the matching key. Importing the .ris or .enw instead would
let EndNote assign its own numbers, and every #n in the text would then point at
whatever happened to land on that number.

Workflow for the author: import the .xml into EndNote, open the built docx in
Word with that library open, then Update Citations and Bibliography.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
OUTDIR = ROOT / "manuscript/endnote"

ENTRY = re.compile(
    r"^(?P<n>\d+)\.\s+(?P<authors>.+?)\.\s+(?P<title>.+?[.?!])\s+"
    r"\*(?P<journal>[^*]+?)\.?\*\s*"
    r"(?P<year>\d{4})?;?(?P<volume>[^(;:]*)?(?:\((?P<issue>[^)]*)\))?"
    r"(?::(?P<pages>[^.]*))?\.?\s*(?:doi:(?P<doi>10\.\S+?))?\.?$"
)


def parse_references(text: str) -> dict[int, dict]:
    block = text.split("\n## References", 1)[1]
    block = re.split(r"^## ", block, maxsplit=1, flags=re.M)[0]
    refs: dict[int, dict] = {}
    for line in block.splitlines():
        if not re.match(r"^\d+\. ", line.strip()):
            continue
        m = ENTRY.match(line.strip())
        if not m:
            raise SystemExit(f"unparsed reference: {line[:70]}")
        d = {k: (v or "").strip() for k, v in m.groupdict().items()}
        d["n"] = int(d["n"])
        d["title"] = d["title"].rstrip(".")
        d["doi"] = d["doi"].rstrip(".") if d["doi"] else ""
        refs[d["n"]] = d
    return refs



def ascii_punct(value: str) -> str:
    """EndNote matches the citation text against the library character for
    character. Crossref returns curly apostrophes and en dashes in names -
    "Al’Khafaji" - and a citation carrying one will not match a record
    carrying the other. Normalise both sides to ASCII."""
    return (value.replace("\u2019", "'").replace("\u2018", "'")
                 .replace("\u201c", '"').replace("\u201d", '"')
                 .replace("\u2013", "-").replace("\u2014", "-")
                 .replace("\u00a0", " "))


def surname(authors: str) -> str:
    """First author's family name, as EndNote matches on it."""
    first = authors.split(",")[0].strip()
    first = re.sub(r"\bet al\b\.?", "", first).strip()
    parts = first.split()
    # entries read "Barrett T", "Arzalluz-Luque A", "O'Keefe S": family first
    return ascii_punct(parts[0] if parts else first)


def expand(spec: str) -> list[int]:
    out: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if re.search(r"[-–]", part):
            lo, hi = re.split(r"[-–]", part)[:2]
            out.extend(range(int(lo), int(hi) + 1))
        elif part:
            out.append(int(part))
    return out


def temporary_citation(spec: str, refs: dict[int, dict]) -> str:
    items = []
    for n in expand(spec):
        r = refs[n]
        items.append(f"{surname(r['authors'])}, {r['year']} #{n}")
    return "{" + "; ".join(items) + "}"


def endnote_xml(refs: dict[int, dict]) -> str:
    rows = []
    for n in sorted(refs):
        r = refs[n]
        names = [a.strip() for a in r["authors"].split(",") if a.strip()]
        names = [a for a in names if not re.fullmatch(r"et al\.?", a)]
        authors = "".join(
            f"<author>{escape(ascii_punct(a))}</author>" for a in names
        )
        rows.append(
            "<record>"
            f"<rec-number>{n}</rec-number>"
            '<ref-type name="Journal Article">17</ref-type>'
            f"<contributors><authors>{authors}</authors></contributors>"
            f"<titles><title>{escape(r['title'])}</title>"
            f"<secondary-title>{escape(r['journal'])}</secondary-title></titles>"
            f"<periodical><full-title>{escape(r['journal'])}</full-title></periodical>"
            + (f"<pages>{escape(r['pages'])}</pages>" if r["pages"] else "")
            + (f"<volume>{escape(r['volume'])}</volume>" if r["volume"] else "")
            + (f"<number>{escape(r['issue'])}</number>" if r["issue"] else "")
            + f"<dates><year>{escape(r['year'])}</year></dates>"
            + (f"<electronic-resource-num>{escape(r['doi'])}</electronic-resource-num>"
               if r["doi"] else "")
            + "</record>"
        )
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<xml><records>\n'
            + "\n".join(rows) + "\n</records></xml>\n")


def main() -> int:
    text = DRAFT.read_text(encoding="utf-8")
    refs = parse_references(text)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    head, sep, tail = text.partition("\n## References")
    # drop the numbered list itself; keep whatever sections follow it
    rest = re.split(r"^## ", tail, maxsplit=1, flags=re.M)
    after = ("\n## " + rest[1]) if len(rest) > 1 else ""

    swapped = 0

    def sub(m: re.Match) -> str:
        nonlocal swapped
        swapped += 1
        return temporary_citation(m.group(1), refs)

    body = re.sub(r"\[([\d]+(?:\s*[,–-]\s*\d+)*)\]", sub, head + after)

    md = OUTDIR / "scTHREAD_NAR_Database_Issue_endnote.md"
    md.write_text(body, encoding="utf-8")
    (OUTDIR / "scTHREAD_NAR_references.xml").write_text(endnote_xml(refs), encoding="utf-8")

    docx = OUTDIR / "scTHREAD_NAR_Database_Issue_endnote.docx"
    subprocess.run(
        ["pandoc", str(md), "--from=gfm+superscript+fenced_divs+raw_attribute",
         "--to=docx", f"--reference-doc={ROOT}/manuscript/style/reference.docx",
         f"--output={docx}"], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/fix_docx_table_widths.py"),
                    str(docx)], check=True)

    leftover = re.findall(r"\[\d+(?:\s*[,–-]\s*\d+)*\]", body)
    print(f"references parsed: {len(refs)}")
    print(f"citation groups converted: {swapped}")
    print(f"numeric citations remaining: {len(leftover)} {leftover[:3]}")
    for p in (md, OUTDIR / "scTHREAD_NAR_references.xml", docx):
        print(f"  {p.relative_to(ROOT)}  {p.stat().st_size:,} bytes")
    return 1 if leftover else 0


if __name__ == "__main__":
    sys.exit(main())
