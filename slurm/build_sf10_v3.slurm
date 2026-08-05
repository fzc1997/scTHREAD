#!/bin/bash
#SBATCH --job-name=nar-sf10-v3
#SBATCH --partition=BICP_Comm_1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/nar_sf10_v3_%j.out
#SBATCH --error=logs/nar_sf10_v3_%j.err
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export TMPDIR=/dev/shm
python3 scripts/build_sf10_v3_sources.py
