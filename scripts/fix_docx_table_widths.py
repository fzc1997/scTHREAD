#!/usr/bin/env python3
"""Give Table 1 a usable column layout in the built docx.

pandoc emits nine equal columns for a nine-column pipe table and ignores the
relative dash widths in the separator row, so every column came out 0.61 inch.
At that width Word breaks mid-word: the header read "scTHRE / AD",
"isoformAt / las [9]", "Program / matic or bulk access", and the first feature
label wrapped to nine lines of one word each.

The fix is a layout one, so it belongs after pandoc rather than in the prose:
widen the feature column, spread the rest over the full text width, and trim
the default cell margins. Run from the build after the docx is written.

    python3 scripts/fix_docx_table_widths.py manuscript/<file>.docx
"""
from __future__ import annotations

import re
import shutil
import sys
import zipfile
from pathlib import Path

# text width of a US Letter page at the template's margins, in twentieths of a
# point: 12240 - 1440 - 1440
TEXT_WIDTH = 9360
FEATURE_COL = 1900          # the widest resource header is "isoformAtlas [9]",
                            # which needs the rest of the width to avoid a
                            # mid-word break
CELL_MARGIN = 40            # 0.028 inch each side, against the 115 default


def fix(document: str) -> tuple[str, str]:
    match = re.search(r"<w:tbl>.*?</w:tbl>", document, re.S)
    if not match:
        return document, "no table found"
    table = match.group(0)
    n = len(re.findall(r"<w:gridCol\b", table))
    if n < 2:
        return document, f"table has {n} columns, nothing to lay out"

    rest = (TEXT_WIDTH - FEATURE_COL) // (n - 1)
    widths = [FEATURE_COL] + [rest] * (n - 1)

    # the grid, which sets the column boundaries
    i = iter(widths)
    fixed = re.sub(r'<w:gridCol w:w="\d+"\s*/>',
                   lambda _: f'<w:gridCol w:w="{next(i)}"/>', table)

    # and every cell's own width, which Word honours over the grid
    per_row = iter(())

    def cell_width(m: re.Match) -> str:
        nonlocal per_row
        try:
            w = next(per_row)
        except StopIteration:
            per_row = iter(widths)
            w = next(per_row)
        return f'<w:tcW w:w="{w}" w:type="dxa"/>'

    rows = []
    for row in re.split(r"(?<=</w:tr>)", fixed):
        per_row = iter(widths)
        rows.append(re.sub(r'<w:tcW w:w="\d+" w:type="\w+"\s*/>', cell_width, row))
    fixed = "".join(rows)

    # fixed layout, full width, and tight margins
    fixed = re.sub(r'<w:tblW[^/]*/>',
                   f'<w:tblW w:w="{TEXT_WIDTH}" w:type="dxa"/>', fixed)
    if "<w:tblLayout" not in fixed:
        fixed = fixed.replace("</w:tblPr>",
                              '<w:tblLayout w:type="fixed"/></w:tblPr>', 1)
    if "<w:tblCellMar" not in fixed:
        fixed = fixed.replace(
            "</w:tblPr>",
            f'<w:tblCellMar><w:left w:w="{CELL_MARGIN}" w:type="dxa"/>'
            f'<w:right w:w="{CELL_MARGIN}" w:type="dxa"/></w:tblCellMar></w:tblPr>', 1)

    return document.replace(table, fixed, 1), (
        f"{n} columns: feature {FEATURE_COL}, others {rest} twips"
    )


def main() -> int:
    path = Path(sys.argv[1])
    tmp = path.with_suffix(".tmp.docx")
    with zipfile.ZipFile(path) as src:
        names = src.namelist()
        blobs = {n: src.read(n) for n in names}
    document = blobs["word/document.xml"].decode("utf-8")
    document, note = fix(document)
    blobs["word/document.xml"] = document.encode("utf-8")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as out:
        for n in names:
            out.writestr(n, blobs[n])
    shutil.move(tmp, path)
    print(f"table layout: {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
