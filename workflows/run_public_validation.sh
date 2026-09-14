#!/usr/bin/env bash
set -euo pipefail
PY="${PYTHON:-python3}"
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
${PY} -m compileall -q database evidence_layers figures
${PY} evidence_layers/biological_unit_validation/validate_results.py
${PY} - <<'PY'
import csv
from pathlib import Path
p = Path('database/run_manifest/processing_manifest_public.tsv')
rows = list(csv.DictReader(p.open(), delimiter='\t'))
assert len(rows) == 453
assert all(not any(k.endswith('_path') for k in row) for row in rows)
print({'processing_manifest_rows': len(rows), 'status': 'PASS'})
PY
${PY} - <<'PY'
import csv
from pathlib import Path
p = Path('database/frozen_release/release_manifest.tsv')
rows = list(csv.DictReader(p.open(), delimiter='\t'))
print({'frozen_release_rows': len(rows)})
assert len(rows) == 453
PY
echo "public package validation PASS"
