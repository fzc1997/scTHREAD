#!/bin/bash
#SBATCH --job-name=nar-sf10-render
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=logs/nar_sf10_render_%j.out
#SBATCH --error=logs/nar_sf10_render_%j.err
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 figures/scripts/render_nar_sf10_v3.py
echo "--- printed numbers ---"
pdftotext figures/NAR_SF10.pdf - | tr '\n' ' ' | grep -o -E "55,729|36 cell types|25,621|8,367|68,417" | sort | uniq -c
