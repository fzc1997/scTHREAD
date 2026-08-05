#!/usr/bin/env python3
"""Renumber manuscript citations into order of first appearance.

The draft is the single source of truth, so new references are appended to the
list with provisional numbers and this script rewrites both the in-text markers
and the reference list into NAR's citation order. Run it after adding or moving
any citation, then re-run tests/test_nar_manuscript_narrative.py.

Scan order = main text (before "## References"), then the material printed after
the list (tables, figure legends), which mirrors the reading order of the article.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

DEFAULT = Path(__file__).resolve().parents[1] / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
CITE = re.compile(r"\[(\d+(?:\s*[–,-]\s*\d+)*)\]")


def expand(token: str) -> list[int]:
    numbers: list[int] = []
    for part in token.split(","):
        part = part.strip()
        if re.search(r"[–-]", part):
            low, high = re.split(r"\s*[–-]\s*", part)
            numbers.extend(range(int(low), int(high) + 1))
        else:
            numbers.append(int(part))
    return numbers


def collapse(numbers: list[int]) -> str:
    numbers = sorted(set(numbers))
    parts, index = [], 0
    while index < len(numbers):
        end = index
        while end + 1 < len(numbers) and numbers[end + 1] == numbers[end] + 1:
            end += 1
        if end - index >= 2:
            parts.append(f"{numbers[index]}–{numbers[end]}")
        else:
            parts.extend(str(numbers[i]) for i in range(index, end + 1))
        index = end + 1
    return "[" + ",".join(parts) + "]"


def renumber_once(text: str) -> tuple[str, int]:
    """One citation-order pass; returns the rewritten text and how many moved."""
    head, refs_block = text.split("## References\n", 1)
    entries = re.findall(r"^(\d+)\.\s+(.*)$", refs_block, re.MULTILINE)
    if not entries:
        raise SystemExit("no numbered references found")
    listed = {int(number): body for number, body in entries}

    # The reference list is the run of numbered lines that opens the block; the
    # tail is everything after it (Table 1 and the figure legends live there).
    # Locating the end by searching for the highest number breaks as soon as that
    # string also occurs in the tail, and silently drops those sections.
    lines = refs_block.split("\n")
    last_entry = -1
    for i, line in enumerate(lines):
        if re.match(r"^\d+\.\s", line):
            last_entry = i
        elif last_entry != -1 and line.strip() and not line.startswith(" "):
            break
    if last_entry == -1:
        raise SystemExit("could not delimit the reference list")
    tail = "\n" + "\n".join(lines[last_entry + 1:])

    order: list[int] = []
    for source in (head, tail):
        for match in CITE.finditer(source):
            for number in expand(match.group(1)):
                if number not in order:
                    order.append(number)

    missing = sorted(set(listed) - set(order))
    unknown = sorted(set(order) - set(listed))
    if unknown:
        raise SystemExit(f"cited but not in the reference list: {unknown}")
    if missing:
        raise SystemExit(f"listed but never cited: {missing}")

    mapping = {old: new for new, old in enumerate(order, start=1)}
    rewrite = lambda source: CITE.sub(  # noqa: E731
        lambda m: collapse([mapping[n] for n in expand(m.group(1))]), source
    )
    new_head, new_tail = rewrite(head), rewrite(tail)
    new_list = "\n".join(
        f"{new}. {listed[old]}" for old, new in sorted(mapping.items(), key=lambda kv: kv[1])
    )
    # new_tail already starts at the blank line that follows the last entry
    rebuilt = f"{new_head}## References\n\n{new_list}{new_tail}"
    return rebuilt, sum(1 for old, new in mapping.items() if old != new)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript", type=Path, default=DEFAULT)
    parser.add_argument("--check", action="store_true",
                        help="fail if renumbering would change the file")
    args = parser.parse_args()

    original = args.manuscript.read_text(encoding="utf-8")
    # Rendering sorts the numbers inside a citation group, which can itself change
    # the order of first appearance, so iterate to a fixed point.
    text, total_moved = original, 0
    for _ in range(10):
        rebuilt, moved = renumber_once(text)
        total_moved += moved
        if rebuilt == text:
            break
        text = rebuilt
    else:
        raise SystemExit("citation order did not converge after 10 passes")

    n_refs = len(re.findall(r"^\d+\.\s", text.split("## References\n", 1)[1], re.MULTILINE))
    if text == original:
        print(f"references already in citation order ({n_refs} entries)")
        return
    if args.check:
        raise SystemExit("references are NOT in citation order; run without --check")
    args.manuscript.write_text(text, encoding="utf-8")
    print(f"renumbered {n_refs} references to citation order")


if __name__ == "__main__":
    main()
