#!/bin/bash
#SBATCH --job-name=nar-fig1-snapshot
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/nar_fig1_snapshot_%j.out
#SBATCH --error=logs/nar_fig1_snapshot_%j.err

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python_bin=python3

"$python_bin" figures/scripts/export_fig1_catalog_stats.py
"$python_bin" figures/scripts/render_nar_bio.py --fig 1 --stem1 NAR_Fig1_v2

echo "--- printed numbers on the rebuilt Fig.1 ---"
txt=$(pdftotext figures/NAR_Fig1_v2.pdf - | tr '\n' ' ')
for want in "453" "923,389" "34"; do
  echo "$txt" | grep -qF "$want" || { echo "MISSING from Fig.1: $want"; exit 1; }
  echo "  present: $want"
done
for stale in "434" "850,938" "845,746" "845,781" "941,733" "469" "478"; do
  if echo "$txt" | grep -qF "$stale"; then echo "STALE number still on Fig.1: $stale"; exit 1; fi
done
echo "Fig.1 carries the 20260803 release numbers and no stale ones"
