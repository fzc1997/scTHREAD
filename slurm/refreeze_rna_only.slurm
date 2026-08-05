#!/bin/bash
#SBATCH --job-name=nar-refreeze-rna
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nar_refreeze_rna_%j.out
#SBATCH --error=logs/nar_refreeze_rna_%j.err

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
py=python3

echo "### completeness tiers (S11)"
"$py" scripts/build_release_completeness_tiers.py
echo "### file contract (S12)"
"$py" scripts/build_release_feature_contract.py
echo "### molecular inventory (S13)"
"$py" scripts/build_release_molecular_inventory.py
echo "### supplementary workbook (S1-S13)"
"$py" scripts/build_nar_supplementary_tables.py
echo "DONE"
