# scTHREAD reproducibility package

This directory contains the analysis entry points needed to reproduce the
evidence layers described in the manuscript and the biological-unit-aware
cell-context analyses.

## Entry points

- `reducers/diu_reduce.py`: transcript-by-cell-barcode matrix to gene/transcript/cell-type counts;
- `reducers/apa_reduce.py`: supported 3′ endpoints to PolyASite assignments;
- `reducers/ase_celltype_reduce.py`: cell-type allele-aware reduction;
- `reducers/f2_junction_celltype.py`: junction × cell-type molecule/read reduction;
- `reducers/diu_test.py`, `apa_test.py`, `ase_interaction_test.py`: layer-specific tests and multiple-testing correction;
- `biological_unit_validation/run_celltype_permutation.py`: source-restricted DIU/APA/ASE permutation analysis;
- `biological_unit_validation/study_stratified_replication.py`: study-stratified follow-up;
- `biological_unit_validation/validate_results.py`: independent output and manifest validation.
- `processing_manifest_schema.tsv`: field-level definitions for the run manifest and its counting-unit boundaries.

## Fixed manuscript scope

All release-facing examples and documentation use the manuscript snapshot of
**453 run records, 34 datasets and 923,389 cells**. Run-level pipeline counts,
library-level records and rolling portal additions are separate quantities and
must not replace this scope.

## Examples

`examples/` contains small, real output examples for isoform, three-axis
DIU/APA/ASE summary, junction and allele-aware evidence. They are output
fixtures, not simulated data.

## External inputs

Raw reads and reference bundles are not redistributed. Set the environment
variables in `ENVIRONMENT.md`, use the source accessions in the frozen release
manifest, and record the versions/checksums of external resources before a
rerun. The scripts are designed to fail closed when a required resource is
missing.
