#!/bin/bash
#SBATCH --job-name=nar-sf-refresh
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nar_sf_refresh_%j.out
#SBATCH --error=logs/nar_sf_refresh_%j.err
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 figures/scripts/render_nar_sf_compact.py
echo "--- composition figure printed numbers ---"
# The 2026-08-04 compaction made the composition figure SF1 (it was SF3).
txt=$(pdftotext figures/NAR_SF1.pdf - | tr '\n' ' ')
for want in "453" "923,389" "34"; do
  echo "$txt" | grep -qF "$want" || { echo "MISSING from the composition figure: $want"; exit 1; }
  echo "  present: $want"
done
for stale in "434" "850,938" "845,746" "845,781" "469" "452,213" "1,045,231" "isoquant_cells" "annotation_cells" "barcodes_raw" "Labeled analysis subset"; do
  if echo "$txt" | grep -qF "$stale"; then echo "STALE or internal wording on the composition figure: $stale"; exit 1; fi
done
echo "the composition figure carries the release numbers, no stale values and no internal column names"
