#!/bin/bash
#SBATCH --job-name=nar-manuscript-docx
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/nar_manuscript_docx_%j.out
#SBATCH --error=logs/nar_manuscript_docx_%j.err

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python_bin=python3
"$python_bin" \
  tests/test_nar_manuscript_narrative.py
# Prove the cross-reference checker still rejects the defects it is meant to
# catch before trusting it to pass the draft.
"$python_bin" tests/test_check_crossrefs.py
"$python_bin" scripts/check_crossrefs.py
"$python_bin" scripts/check_panel_coverage.py
# figures were never scanned for superseded denominators; SF3 shipped a
# 469-run screenshot past every check in this gate
"$python_bin" scripts/check_figure_numbers.py
# Callout checks prove references resolve, not that the prose is true. This
# recomputes the anchor denominators from the release manifest and the served
# portal tables and rejects claims that contradict them.
"$python_bin" tests/test_check_anchor_numbers.py
"$python_bin" scripts/check_anchor_numbers.py
# the workbook carries the same denominators as the text; left out of this
# gate it silently kept shipping the superseded 469-run/30-study candidate
"$python_bin" tests/test_nar_supplementary_tables.py \
  manuscript/supplementary/scTHREAD_Supplementary_Tables.xlsx
pandoc manuscript/scTHREAD_NAR_Database_Issue_draft.md \
  --from=gfm+superscript+fenced_divs+raw_attribute \
  --to=docx \
  --reference-doc=manuscript/style/reference.docx \
  --output=manuscript/scTHREAD_NAR_Database_Issue_draft.docx
# pandoc gives a nine-column table nine equal 0.61-inch columns and breaks the
# headers mid-word; this is layout, so it is applied after the prose is written
"$python_bin" scripts/fix_docx_table_widths.py \
  manuscript/scTHREAD_NAR_Database_Issue_draft.docx

docx=manuscript/scTHREAD_NAR_Database_Issue_draft.docx
check_dir=$(mktemp -d)
trap 'rm -rf "$check_dir"' EXIT
unzip -t "$docx"
pandoc "$docx" --to=plain --output="$check_dir/manuscript.txt"
"$python_bin" - "$check_dir/manuscript.txt" <<'PY'
import re
import sys
from pathlib import Path

text = re.sub(r"\s+", " ", Path(sys.argv[1]).read_text())
required = (
    "Supplementary Tables S1–S15",
    "12,471,531 positive run-by-transcript detections",
    "Table 1. Feature comparison between scTHREAD",
    "13 of 13 donors",
    "the three-locus pan-myeloid model was not supported",
    "Every dataset is carried to the same isoform-level end product",
    "cells carry a cell-type assignment wherever cell type is defined",
    "Gene Expression Omnibus, the European Nucleotide Archive",
    "117 records and 736,853 cells",
    "its 68,417 cells are part of the 923,389 reported in the main text",
)
for phrase in required:
    if phrase not in text:
        raise AssertionError(f"DOCX text missing: {phrase}")
for phrase in (
    "Supplementary Tables S1–S10 are provided",
    "469 IsoQuant-complete",
    "frozen candidate manifest",
    "author-provided LAN deployment",
    "Cell-type-specific Malat1 DTU",
    "D291–D303",
    "941,733",
    "Database content and results",
    "membership-only",
    "completeness tier",
    "one-well records",
    "850,938",
    "434 run records",
    "434-row manifest",
):
    if phrase in text:
        raise AssertionError(f"Forbidden legacy wording remains: {phrase}")
print("DOCX NORMALIZED-TEXT QA PASS")
PY
sha256sum "$docx"
