#!/bin/bash
#SBATCH --job-name=nar-tx-structure
#SBATCH --partition=BICP_Comm_1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/nar_tx_structure_%j.out
#SBATCH --error=logs/nar_tx_structure_%j.err
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 scripts/measure_transcript_structure.py --max-runs 8
