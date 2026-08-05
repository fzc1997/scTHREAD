#!/usr/bin/env python3
"""Export the manuscript reference list to EndNote (.enw) and RIS (.ris).

The list in the markdown abbreviates author lists with "et al.", which is right
for the paper and wrong for a reference library, so every entry is enriched from
Crossref by its DOI. Enrichment is best-effort: an entry whose DOI cannot be
resolved falls back to what the markdown states, and the run reports which ones
those were rather than silently shipping a thin record.

    python3 scripts/export_references.py [--offline] [--outdir manuscript/references]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
MAILTO = "fuzc@blsa.com.cn"  # Crossref asks for a contact; it buys the polite pool

ENTRY = re.compile(
    # a title may end in "?" or "!", not only ".", and one of these does:
    # "Single-cell RNAseq for the study of isoforms-how is that possible?"
    r"^(?P<n>\d+)\.\s+(?P<authors>.+?)\.\s+(?P<title>.+?[.?!])\s+"
    r"\*(?P<journal>[^*]+?)\.?\*\s*"
    r"(?P<year>\d{4})?;?(?P<volume>[^(;:]*)?(?:\((?P<issue>[^)]*)\))?"
    r"(?::(?P<pages>[^.]*))?\.?\s*"
    r"(?:doi:(?P<doi>10\.\S+?))?\.?$",
    re.S,
)


def parse(text: str) -> list[dict]:
    block = text.split("\n## References", 1)[1]
    block = re.split(r"^## ", block, maxsplit=1, flags=re.M)[0]
    out = []
    for line in block.splitlines():
        line = line.strip()
        if not re.match(r"^\d+\. ", line):
            continue
        m = ENTRY.match(line)
        if not m:
            out.append({"n": int(line.split(".", 1)[0]), "raw": line, "parse_failed": True})
            continue
        d = {k: (v.strip() if isinstance(v, str) else v) for k, v in m.groupdict().items()}
        d["n"] = int(d["n"])
        # the sentence-ending period is punctuation of the reference list, not
        # part of the title; "?" and "!" are part of the title and stay
        if d.get("title", "").endswith("."):
            d["title"] = d["title"][:-1]
        d["raw"] = line
        if d.get("doi"):
            d["doi"] = d["doi"].rstrip(".")
        out.append(d)
    return out


def crossref(doi: str) -> dict | None:
    url = f"https://api.crossref.org/works/{doi}?mailto={MAILTO}"
    try:
        with urllib.request.urlopen(url, timeout=25) as fh:
            return json.load(fh)["message"]
    except (urllib.error.URLError, TimeoutError, KeyError, json.JSONDecodeError, OSError):
        return None


def enrich(ref: dict, msg: dict) -> dict:
    ref = dict(ref)
    authors = []
    for a in msg.get("author", []) or []:
        family, given = a.get("family"), a.get("given")
        if family and given:
            authors.append(f"{family}, {given}")
        elif family:
            authors.append(family)
        elif a.get("name"):
            authors.append(a["name"])
    if authors:
        ref["author_list"] = authors
    for key, path in (("title", "title"), ("journal", "container-title")):
        value = msg.get(path)
        if isinstance(value, list) and value:
            ref[key] = value[0]
    for key in ("volume", "issue", "page"):
        if msg.get(key):
            ref["pages" if key == "page" else key] = msg[key]
    issued = (msg.get("issued") or {}).get("date-parts") or [[None]]
    # NAR Database-issue entries are cited by their issue year; Crossref reports
    # the advance-access date, so never let it overwrite a year already parsed
    if issued[0][0] and not ref.get("year"):
        ref["year"] = str(issued[0][0])
    if msg.get("abstract"):
        ref["abstract"] = re.sub(r"<[^>]+>", "", msg["abstract"]).strip()
    ref["enriched"] = True
    return ref


def authors_of(ref: dict) -> list[str]:
    if ref.get("author_list"):
        return ref["author_list"]
    raw = (ref.get("authors") or "").replace(", et al", "").replace(" et al", "")
    return [a.strip() for a in raw.split(",") if a.strip()]


def write_enw(refs: list[dict], path: Path) -> None:
    """EndNote's own tagged format."""
    lines = []
    for ref in refs:
        lines.append("%0 Journal Article")
        for author in authors_of(ref):
            lines.append(f"%A {author}")
        if ref.get("title"):
            lines.append(f"%T {ref['title']}")
        if ref.get("journal"):
            lines.append(f"%J {ref['journal']}")
        if ref.get("year"):
            lines.append(f"%D {ref['year']}")
        if ref.get("volume"):
            lines.append(f"%V {ref['volume']}")
        if ref.get("issue"):
            lines.append(f"%N {ref['issue']}")
        if ref.get("pages"):
            lines.append(f"%P {ref['pages']}")
        if ref.get("doi"):
            lines.append(f"%R {ref['doi']}")
            lines.append(f"%U https://doi.org/{ref['doi']}")
        if ref.get("abstract"):
            lines.append(f"%X {ref['abstract']}")
        lines.append(f"%M {ref['n']}")  # keeps the manuscript's citation number
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def write_ris(refs: list[dict], path: Path) -> None:
    """RIS, which every manager including EndNote imports."""
    lines = []
    for ref in refs:
        lines.append("TY  - JOUR")
        for author in authors_of(ref):
            lines.append(f"AU  - {author}")
        if ref.get("title"):
            lines.append(f"TI  - {ref['title']}")
        if ref.get("journal"):
            lines.append(f"JO  - {ref['journal']}")
        if ref.get("year"):
            lines.append(f"PY  - {ref['year']}")
        if ref.get("volume"):
            lines.append(f"VL  - {ref['volume']}")
        if ref.get("issue"):
            lines.append(f"IS  - {ref['issue']}")
        if ref.get("pages"):
            pages = str(ref["pages"]).replace("–", "-")
            start, _, end = pages.partition("-")
            lines.append(f"SP  - {start.strip()}")
            if end.strip():
                lines.append(f"EP  - {end.strip()}")
        if ref.get("doi"):
            lines.append(f"DO  - {ref['doi']}")
            lines.append(f"UR  - https://doi.org/{ref['doi']}")
        if ref.get("abstract"):
            lines.append(f"AB  - {ref['abstract']}")
        lines.append(f"ID  - ref{ref['n']}")
        lines.append("ER  - ")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=ROOT / "manuscript/references")
    ap.add_argument("--offline", action="store_true", help="skip Crossref enrichment")
    args = ap.parse_args()

    refs = parse(MANUSCRIPT.read_text(encoding="utf-8"))
    print(f"parsed {len(refs)} references")

    failed_parse = [r["n"] for r in refs if r.get("parse_failed")]
    no_doi = [r["n"] for r in refs if not r.get("doi")]
    not_enriched: list[int] = []

    if not args.offline:
        for ref in refs:
            if not ref.get("doi"):
                continue
            msg = crossref(ref["doi"])
            if msg is None:
                not_enriched.append(ref["n"])
            else:
                refs[refs.index(ref)] = enrich(ref, msg)
            time.sleep(0.12)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_enw(refs, args.outdir / "scTHREAD_NAR_references.enw")
    write_ris(refs, args.outdir / "scTHREAD_NAR_references.ris")

    enriched = sum(1 for r in refs if r.get("enriched"))
    print(f"enriched from Crossref: {enriched}/{len(refs)}")
    if failed_parse:
        print(f"  UNPARSED entries (exported from the raw line): {failed_parse}")
    if no_doi:
        print(f"  no DOI in the manuscript: {no_doi}")
    if not_enriched:
        print(f"  DOI present but Crossref did not answer: {not_enriched}")
    for name in ("scTHREAD_NAR_references.enw", "scTHREAD_NAR_references.ris"):
        p = args.outdir / name
        print(f"  {p.relative_to(ROOT)}  {p.stat().st_size:,} bytes")
    return 1 if failed_parse else 0


if __name__ == "__main__":
    sys.exit(main())
