#!/usr/bin/env python3
"""Apply paragraph-level tracked edits to a docx that carries EndNote fields.

docx-track-changes rewrites changed paragraphs word by word, which zeroes the
instrText of any EndNote field inside them - measured elsewhere as 72 fields to
0. This file has 77 EN.CITE fields, so that tool is off limits.

This does the two operations that are safe on a field-heavy document:

  replace  mark a whole paragraph deleted and insert a replacement after it
  delete   mark a whole paragraph deleted

A deletion wraps each existing run in <w:del> and renames w:t to w:delText
without touching run structure, so a field inside a deleted paragraph keeps its
begin / separate / end runs intact: Word shows the citation struck through and
drops it cleanly on accept. Nothing inside a kept paragraph is rewritten.

    python3 scripts/surgical_tracked_edit.py IN.docx OUT.docx edits.json \
        --author "Zhi-Can Fu"

edits.json is a list of {"match": "<unique substring>", "action": "replace",
"text": "<new paragraph text>"} or {"match": ..., "action": "delete"}.
"""
from __future__ import annotations

import argparse
import copy
import json
import shutil
import sys
import zipfile
from pathlib import Path

from lxml import etree

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def para_text(p) -> str:
    return "".join(t.text or "" for t in p.iter(W + "t"))


def mark_deleted(p, author: str, date: str, ident: list[int]) -> None:
    """Wrap every run of a paragraph in <w:del>, preserving run structure."""
    for run in list(p.findall(W + "r")):
        idx = list(p).index(run)
        dele = etree.Element(W + "del")
        dele.set(W + "id", str(ident[0])); ident[0] += 1
        dele.set(W + "author", author)
        dele.set(W + "date", date)
        p.remove(run)
        for t in run.iter(W + "t"):
            t.tag = W + "delText"
            t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        dele.append(run)
        p.insert(idx, dele)
    # and mark the paragraph mark itself deleted, so the paragraph disappears
    pPr = p.find(W + "pPr")
    if pPr is None:
        pPr = etree.SubElement(p, W + "pPr")
        p.remove(pPr); p.insert(0, pPr)
    rPr = pPr.find(W + "rPr")
    if rPr is None:
        rPr = etree.SubElement(pPr, W + "rPr")
    dele = etree.SubElement(rPr, W + "del")
    dele.set(W + "id", str(ident[0])); ident[0] += 1
    dele.set(W + "author", author)
    dele.set(W + "date", date)


def inserted_paragraph(template, text: str, author: str, date: str,
                       ident: list[int]):
    """A new paragraph in the template's style, wholly marked as inserted."""
    p = etree.Element(W + "p")
    tmpl_pPr = template.find(W + "pPr")
    if tmpl_pPr is not None:
        pPr = copy.deepcopy(tmpl_pPr)
        for tag in (W + "rPr",):
            node = pPr.find(tag)
            if node is not None:
                for child in list(node):
                    if child.tag in (W + "del", W + "ins"):
                        node.remove(child)
        p.append(pPr)
    ins = etree.SubElement(p, W + "ins")
    ins.set(W + "id", str(ident[0])); ident[0] += 1
    ins.set(W + "author", author)
    ins.set(W + "date", date)
    run = etree.SubElement(ins, W + "r")
    # take the character formatting of the paragraph's first plain run
    first = template.find(W + "r")
    if first is not None:
        rPr = first.find(W + "rPr")
        if rPr is not None:
            run.append(copy.deepcopy(rPr))
    t = etree.SubElement(run, W + "t")
    t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    t.text = text
    return p


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("edits", type=Path)
    ap.add_argument("--author", required=True)
    ap.add_argument("--date", default="2026-08-05T00:00:00Z")
    args = ap.parse_args()

    edits = json.loads(args.edits.read_text(encoding="utf-8"))
    with zipfile.ZipFile(args.source) as z:
        names = z.namelist()
        blobs = {n: z.read(n) for n in names}
    root = etree.fromstring(blobs["word/document.xml"])
    body = root.find(W + "body")
    ident = [9000]

    before = {k: blobs["word/document.xml"].decode().count(k)
              for k in ("w:fldChar", "w:instrText", "EN.CITE")}

    applied = []
    for edit in edits:
        target = None
        for p in body.findall(W + "p"):
            if edit["match"] in para_text(p):
                target = p
                break
        if target is None:
            raise SystemExit(f"no paragraph matches: {edit['match'][:60]!r}")
        idx = list(body).index(target)
        if edit["action"] == "replace":
            new = inserted_paragraph(target, edit["text"], args.author,
                                     args.date, ident)
            body.insert(idx + 1, new)
        mark_deleted(target, args.author, args.date, ident)
        applied.append((edit["action"], edit["match"][:44]))

    blobs["word/document.xml"] = etree.tostring(
        root, xml_declaration=True, encoding="UTF-8", standalone=True)
    after = {k: blobs["word/document.xml"].decode().count(k)
             for k in ("w:fldChar", "w:instrText", "EN.CITE")}

    tmp = args.output.with_suffix(".tmp.docx")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
        for n in names:
            out.writestr(n, blobs[n])
    shutil.move(tmp, args.output)

    for action, m in applied:
        print(f"  {action:<8} {m}")
    ok = all(before[k] == after[k] for k in before)
    for k in before:
        print(f"  {k}: {before[k]} -> {after[k]} {'ok' if before[k]==after[k] else 'CHANGED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
