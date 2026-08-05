#!/usr/bin/env python3
"""Resolve manuscript references against the Crossref REST API."""

from __future__ import annotations

import argparse
import json
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANUSCRIPT = ROOT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
DEFAULT_OUTPUT_DIR = ROOT / "tables/citation_audit"
API = "https://api.crossref.org"
USER_AGENT = "scTHREAD-reference-audit/1.0 (mailto:fuzhican@163.com)"


def normalize_title(value: str) -> str:
    value = unicodedata.normalize("NFKD", value)
    value = re.sub(r"<[^>]+>", " ", value)
    value = re.sub(r"[^a-z0-9]+", " ", value.lower())
    return re.sub(r"\s+", " ", value).strip()


def similarity(left: str, right: str) -> float:
    return SequenceMatcher(None, normalize_title(left), normalize_title(right)).ratio()


def request_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def extract_references(text: str) -> list[dict[str, object]]:
    block = re.split(r"\n## ", text.split("## References\n", 1)[1], maxsplit=1)[0]
    rows = []
    pattern = re.compile(
        r"^(\d+)\.\s+(.+?)\.\s+(.+?)[.?!]\s+\*([^*]+)\*\.?\s+(.+)$"
    )
    for line in block.splitlines():
        if not re.match(r"^\d+\.", line):
            continue
        match = pattern.match(line)
        if not match:
            raise RuntimeError(f"Cannot parse reference: {line}")
        number, authors, title, journal, tail = match.groups()
        doi_match = re.search(r"(10\.\d{4,9}/\S+?)(?:[.]?$)", tail, re.I)
        doi = doi_match.group(1).rstrip(".") if doi_match else ""
        rows.append(
            {
                "reference_number": int(number),
                "manuscript_authors": authors,
                "manuscript_title": title,
                "manuscript_journal": journal.rstrip("."),
                "manuscript_tail": tail,
                "manuscript_doi": doi.lower(),
                "original_line": line,
            }
        )
    if len(rows) < 20:
        raise RuntimeError(f"Expected at least 20 references, found {len(rows)}")
    return rows


def resolve(row: dict[str, object]) -> tuple[dict, str, float]:
    doi = str(row["manuscript_doi"])
    if doi:
        encoded = urllib.parse.quote(doi, safe="")
        payload = request_json(f"{API}/works/{encoded}")
        item = payload["message"]
        return item, f"{API}/works/{encoded}", similarity(
            str(row["manuscript_title"]), item.get("title", [""])[0]
        )

    query = urllib.parse.urlencode(
        {
            "query.title": str(row["manuscript_title"]),
            "query.author": str(row["manuscript_authors"]).split(",", 1)[0],
            "rows": 5,
        }
    )
    url = f"{API}/works?{query}"
    items = request_json(url)["message"]["items"]
    if not items:
        query = urllib.parse.urlencode(
            {
                "query.bibliographic": str(row["manuscript_title"]),
                "rows": 5,
            }
        )
        url = f"{API}/works?{query}"
        items = request_json(url)["message"]["items"]
    if not items:
        raise RuntimeError(f"No Crossref results for reference {row['reference_number']}")
    scored = [
        (similarity(str(row["manuscript_title"]), item.get("title", [""])[0]), item)
        for item in items
    ]
    score, item = max(scored, key=lambda pair: pair[0])
    return item, url, score


def first(value: object) -> str:
    if isinstance(value, list):
        result = str(value[0]) if value else ""
    else:
        result = str(value or "")
    return re.sub(r"\s+", " ", result).strip()


def year(item: dict) -> int | str:
    for key in ("published-print", "published", "published-online", "issued"):
        parts = item.get(key, {}).get("date-parts", [])
        if parts and parts[0]:
            return parts[0][0]
    return ""


def author_text(item: dict) -> str:
    authors = []
    for author in item.get("author", []):
        family = author.get("family", "").strip()
        given = author.get("given", "").strip()
        if not family:
            continue
        initials = []
        for token in re.split(r"\s+", given):
            components = [part for part in token.split("-") if part]
            if not components:
                continue
            initials.append("-".join(part[0].upper() for part in components))
        authors.append(f"{family} {''.join(initials)}".strip())
    if len(authors) > 4:
        return ", ".join(authors[:3]) + ", et al"
    return ", ".join(authors)


def formatted_reference(number: int, item: dict) -> str:
    authors = author_text(item)
    title = first(item.get("title")).rstrip(".")
    journal = first(item.get("short-container-title")) or first(
        item.get("container-title")
    )
    if not journal and item.get("institution"):
        journal = first(item["institution"][0].get("name"))
    citation = f"{number}. {authors}. {title}. *{journal}.* {year(item)}"
    volume = str(item.get("volume", "")).strip()
    issue = str(item.get("issue", "")).strip()
    page = str(item.get("page", "")).strip().replace("-", "–")
    article = str(item.get("article-number", "")).strip()
    locator = page or article
    if volume:
        citation += f";{volume}"
        if issue:
            citation += f"({issue})"
        if locator:
            citation += f":{locator}"
    elif locator:
        citation += f":{locator}"
    citation += "."
    doi = str(item.get("DOI", "")).lower()
    if doi:
        citation += f" doi:{doi}."
    return citation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=DEFAULT_MANUSCRIPT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    references = extract_references(args.manuscript.read_text())
    if args.limit:
        references = references[: args.limit]
    output_rows = []
    raw_records = []
    suggested = []
    for row in references:
        item, url, score = resolve(row)
        status = "verified" if score >= 0.95 else "suspicious" if score >= 0.80 else "manual_needed"
        citation = formatted_reference(int(row["reference_number"]), item)
        output_rows.append(
            {
                **row,
                "status": status,
                "title_similarity": round(score, 6),
                "crossref_doi": str(item.get("DOI", "")).lower(),
                "crossref_title": first(item.get("title")),
                "crossref_authors": author_text(item),
                "crossref_journal": first(item.get("container-title")),
                "crossref_short_journal": first(item.get("short-container-title")),
                "crossref_year": year(item),
                "crossref_volume": item.get("volume", ""),
                "crossref_issue": item.get("issue", ""),
                "crossref_page": item.get("page", ""),
                "crossref_article_number": item.get("article-number", ""),
                "crossref_type": item.get("type", ""),
                "crossref_query_url": url,
                "suggested_reference": citation,
            }
        )
        raw_records.append(
            {
                "reference_number": row["reference_number"],
                "query_url": url,
                "message": item,
            }
        )
        suggested.append(citation)
        time.sleep(0.1)

    frame = pd.DataFrame(output_rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = args.output_dir / "crossref_reference_audit.tsv"
    raw_path = args.output_dir / "crossref_reference_metadata.json"
    suggested_path = args.output_dir / "crossref_suggested_references.md"
    frame.to_csv(audit_path, sep="\t", index=False)
    raw_path.write_text(
        json.dumps(
            {
                "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                "api": API,
                "records": raw_records,
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    )
    suggested_path.write_text(
        "# Crossref-resolved references\n\n" + "\n".join(suggested) + "\n"
    )
    print(frame["status"].value_counts().to_dict())
    print(audit_path)
    print(raw_path)
    print(suggested_path)


if __name__ == "__main__":
    main()
