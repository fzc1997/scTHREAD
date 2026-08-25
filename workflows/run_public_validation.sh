#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
python3 -m compileall -q reproducibility figures/scripts
python3 reproducibility/p0_biological_unit_rerun/validate_results.py
python3 - <<'PY'
import csv
from pathlib import Path
p=Path('tables/processing_manifest_public.tsv')
rows=list(csv.DictReader(p.open(),delimiter='\t'))
assert len(rows)==453
assert all(not any(k.endswith('_path') for k in row) for row in rows)
print({'processing_manifest_rows':len(rows),'status':'PASS'})
PY
