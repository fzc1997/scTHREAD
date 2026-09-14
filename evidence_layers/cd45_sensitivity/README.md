# CD45 read-vs-UMI sensitivity

This directory contains the code and summary outputs for the CD45 read-versus-UMI
sensitivity analysis. The recount used the indexed tagged BAMs from the two
source studies and their run/barcode/cell-type annotation tables; raw BAMs are
not redistributed.

## Contract

- The direct exon-3 splice-choice classifier is the same classifier used for
  the read-level CD45 measurement.
- Read counts are retained as the primary direct-junction evidence.
- UMI class counts use `(corrected CB, UB, splice-choice class)`.
- A CB/UB key supporting more than one splice-choice class is a conflict and is
  excluded from UMI class counts rather than resolved by dominance.
- The same `>=100 informative read` run × cell-type eligibility threshold is
  applied before the read/UMI summary is compared.

## Reproduction entry points

From a package root with the raw tagged BAMs and annotation TSVs supplied via
the manuscript accessions:

```bash
python evidence_layers/cd45_sensitivity/measure_cd45_umi_sensitivity.py \
  --bam-root /path/to/GSE276974/runs \
  --annotation /path/to/GSE276974_cell_annotation.tsv \
  --study GSE276974 \
  --output GSE276974.tsv \
  --manifest GSE276974.manifest.json

python evidence_layers/cd45_sensitivity/summarize_cd45_umi_sensitivity.py \
  --input GSE276974.tsv GSE307660.tsv \
  --output-tsv cd45_umi_sensitivity.tsv \
  --output-json cd45_umi_sensitivity.json
```

The frozen summary contains 31 successfully annotated runs, 312,919
informative reads and 114,610 informative UMI class counts. The within-study
CD45RO ordering is preserved and the maximum absolute change is 2.14 percentage
points. This is a sensitivity check of the qualitative splice-choice pattern;
it is not a claim that every transcript-abundance layer is UMI-insensitive.
