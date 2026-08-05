#!/usr/bin/env bash
# Assemble a clean, publishable scTHREAD repository for the Zenodo deposit.
#
# The working tree lives inside a git repo rooted at the user's home directory
# that also holds several unrelated projects, so it cannot be pushed. This
# copies only what a reader needs to reproduce the paper, by whitelist, and
# skips the large or copyrighted directories listed in
# docs/HANDOFF_zenodo_deposit.md.
set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DST="${1:?usage: build_release_repo.sh <target-dir>}"
[ -e "$DST" ] && { echo "refusing to write into an existing path: $DST" >&2; exit 1; }
mkdir -p "$DST"

copy() {  # copy() <relative path> [pattern]
  local rel="$1" pat="${2:-}"
  [ -e "$SRC/$rel" ] || { echo "  skip (absent): $rel"; return; }
  mkdir -p "$DST/$(dirname "$rel")"
  if [ -n "$pat" ]; then
    mkdir -p "$DST/$rel"
    find "$SRC/$rel" -maxdepth 1 -name "$pat" -type f -exec cp {} "$DST/$rel/" \;
  else
    cp -r "$SRC/$rel" "$DST/$(dirname "$rel")/"
  fi
  echo "  + $rel${pat:+ ($pat)}"
}

echo "assembling into $DST"
copy scripts                     '*.py'
copy scripts                     '*.sh'
copy figures/scripts             '*.py'
copy tests                       '*.py'
copy slurm                       '*.slurm'
copy manuscript/build_manuscript_docx.slurm
copy manuscript/style
# tables/ is 37 MB in full; take only what the shipped scripts actually read.
copy tables/release_20260803
copy tables/fig4_source          '*.tsv'
copy tables                      '*.tsv'
copy tables                      '*.json'
copy tables/release_content      '*.tsv'
copy tables/modality_audit       '*.tsv'
copy tables/p0_biological_unit_rerun '*.tsv'
# the superseded candidate manifest is deliberately NOT shipped: a deposit
# holding two manifests with different counts invites the wrong one to be read.
# audit_manifest_cellcounts.py is the only reader and is a historical audit.
# figure data lives beside the repo, under $SCTHREAD_PROJECT_ROOT
PROJ="${SCTHREAD_PROJECT_ROOT:-/gpfs/home/fuzc/project/scTHREAD}"
for sub in results/paper1/figdata results/paper1/f2_grammar/figdata; do
  if [ -d "$PROJ/$sub" ]; then
    mkdir -p "$DST/$sub"
    find "$PROJ/$sub" -maxdepth 1 -type f \( -name '*.tsv' -o -name '*.json' \) \
      -exec cp {} "$DST/$sub/" \;
    echo "  + $sub ($(find "$DST/$sub" -type f | wc -l) files)"
  else
    echo "  skip (absent): $sub"
  fi
done
copy manuscript/supplementary/scTHREAD_Supplementary_Tables.xlsx
copy docs/Figure_Script_mapping.tsv
# DATABASE_SCALE.md and CD45_CLAIM_CORRECTION.md are internal Chinese working
# notes and DATABASE_SCALE.md still quotes the superseded 434/30/850,938; the
# deposit carries an English scope statement instead. NAR_FIGURE_LEGENDS.md is
# superseded working notes - the submitted legends are in the manuscript.
copy docs/SCOPE_AND_DENOMINATORS.md
copy docs/NAR_SF_LEGENDS.md

# the submission figures, by their delivered names
mkdir -p "$DST/figures"
for f in NAR_Fig2_v4.pdf NAR_Fig3_v8.pdf NAR_Fig4.pdf \
         scTHREAD_graphical_abstract_gpt2_hybrid_v9.pdf \
         NAR_GraphicalAbstract_v2_20260805.pdf; do
  [ -f "$SRC/figures/$f" ] && cp "$SRC/figures/$f" "$DST/figures/" && echo "  + figures/$f"
done
for i in $(seq 1 10); do
  [ -f "$SRC/figures/NAR_SF$i.pdf" ] && cp "$SRC/figures/NAR_SF$i.pdf" "$DST/figures/"
done
echo "  + figures/NAR_SF1..10.pdf"
[ -f "$SRC/figures/fig1_r15/NAR_Fig1_scTHREAD_editable_icons_v4_textfit_r15.svg" ] && \
  cp "$SRC/figures/fig1_r15/NAR_Fig1_scTHREAD_editable_icons_v4_textfit_r15.svg" \
     "$SRC/figures/fig1_r15/NAR_Fig1_r15.pdf" "$DST/figures/" && echo "  + figures/NAR_Fig1_r15.{svg,pdf}"

cat > "$DST/.gitignore" <<'EOF'
__pycache__/
*.pyc
.ipynb_checkpoints/
logs/
*.log
figures/_staging/
figures/_old_sf/
EOF

cat > "$DST/README.md" <<'EOF'
# scTHREAD

Code and frozen release tables for the scTHREAD NAR Database Issue paper.

- Database: https://scthread.ai4sc.ac.cn/
- Frozen release: `tables/release_20260803/` — 453 run records, 34 datasets, 923,389 cells
- Figure provenance: `docs/Figure_Script_mapping.tsv` (which script builds which figure)
- Scope and denominators: `docs/SCOPE_AND_DENOMINATORS.md`

Raw sequencing data are not redistributed here; accessions are listed in the
manuscript and in `tables/release_20260803/release_manifest.tsv`.

## Licence

TO BE SET BEFORE DEPOSIT — see docs/HANDOFF_zenodo_deposit.md section 5.
EOF

sync
echo
echo "size: $(du -sh "$DST" | cut -f1)   files: $(find "$DST" -type f | wc -l)"
echo
echo "NEXT: review the copied files by hand (credentials, internal paths, third-party material),"
echo "      set the licence in README.md, then git init / add / commit / tag."
