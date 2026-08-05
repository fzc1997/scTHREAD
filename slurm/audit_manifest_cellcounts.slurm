#!/bin/bash
#SBATCH --job-name=nar-cellcount-audit
#SBATCH --partition=BICP_Test
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nar_cellcount_audit_%j.out
#SBATCH --error=logs/nar_cellcount_audit_%j.err
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 scripts/audit_manifest_cellcounts.py
