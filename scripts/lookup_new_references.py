#!/usr/bin/env python3
"""Resolve candidate reference titles through Crossref and emit formatted entries.

Companion to audit_references_crossref.py: that script verifies references that are
already in the manuscript, this one verifies candidates *before* they are added, so
no reference is ever written from memory. Titles that do not reach the similarity
threshold are reported as UNRESOLVED and must not be cited.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_references_crossref import (  # noqa: E402
    API,
    formatted_reference,
    normalize_title,
    request_json,
    similarity,
)

MIN_SIMILARITY = 0.90


def resolve_title(title: str) -> tuple[dict | None, float]:
    query = urllib.parse.urlencode(
        {"query.bibliographic": title, "rows": 5, "select": ",".join([
            "DOI", "title", "author", "container-title", "short-container-title",
            "volume", "issue", "page", "article-number", "issued",
            "published-print", "published-online", "published", "type",
        ])}
    )
    payload = request_json(f"{API}/works?{query}")
    best, best_score = None, 0.0
    for item in payload.get("message", {}).get("items", []):
        candidate = item.get("title") or [""]
        score = similarity(normalize_title(title), normalize_title(candidate[0]))
        if score > best_score:
            best, best_score = item, score
    return best, best_score


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--titles", type=Path, required=True,
                        help="one candidate title per line; blank lines and # ignored")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    titles = [
        line.strip()
        for line in args.titles.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]

    resolved, unresolved = [], []
    for index, title in enumerate(titles, start=1):
        try:
            item, score = resolve_title(title)
        except Exception as exc:  # network / API failure must be visible
            unresolved.append({"title": title, "reason": f"lookup failed: {exc}"})
            print(f"[{index:2d}] ERROR   {title[:70]} :: {exc}")
            continue
        if item is None or score < MIN_SIMILARITY:
            unresolved.append({"title": title, "reason": f"best similarity {score:.2f}"})
            print(f"[{index:2d}] UNRESOLVED ({score:.2f}) {title[:70]}")
            continue
        entry = formatted_reference(0, item)[3:]  # drop the placeholder "0. "
        resolved.append({
            "query_title": title,
            "similarity": round(score, 4),
            "doi": str(item.get("DOI", "")).lower(),
            "formatted": entry,
        })
        print(f"[{index:2d}] OK ({score:.2f}) {entry[:110]}")

    args.out.write_text(
        json.dumps({"resolved": resolved, "unresolved": unresolved}, indent=2),
        encoding="utf-8",
    )
    print(f"\nresolved {len(resolved)} / {len(titles)}; unresolved {len(unresolved)}")
    print(f"written: {args.out}")


if __name__ == "__main__":
    main()
