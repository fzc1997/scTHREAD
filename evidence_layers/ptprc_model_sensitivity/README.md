# PTPRC transcript/model truncation sensitivity (2026-09-13)

This package is a scoped, model-level truncation stress test.
It compares three nested filters on real IsoQuant transcript counts and run-local
GTF models from GSE276974 and GSE307660.

## Filters

- F0: all positive PTPRC models represented in the run GTF.
- F1: F0 models whose annotated terminal coordinate matches a human PolyASite 2.0 cluster on the same chromosome/strand within 25 bp on each side.
- F2: F1 reference models with transcript support level 1/2 or `basic` tag, excluding `mRNA_end_NF`/`cds_end_NF` models.
- F3: PAS-supported run-local novel models with count >=3 (descriptive only); none were scorable in this target subset.

## Results

- 34 records were in scope; 33 had positive PTPRC model counts and were scorable.
- F1 retained 37,880 / 151,856 model counts (median run retention 24.1%).
- F2 retained 14,164 / 151,856 model counts (median run retention 8.7%).
- F0 included 113,522 counts (74.8%) from models >10 kb short of the canonical PTPRC 3′ end and without PAS support under this model-end diagnostic.
- RA-like `ENST00000442510` usage increased in 32/33 runs under both F1 and F2; median paired usage change was 0.170 (F1) and 0.577 (F2). Source-level bootstrap 95% CIs are in `ptprc_model_bootstrap_summary.json`.
- Usage ranks among retained common models were highly concordant (median Spearman rho 1.0), but this does not mean abundance was stable because the denominator changed and low-support models were removed.

## Interpretation boundary

The result supports a material truncation/endpoint-completeness risk for full-length PTPRC model usage. It is a sensitivity/stress test, not proof that every removed model is artifactual: PolyASite support is endpoint plausibility, not full-length proof, and the F2 subset retains only a small fraction of counts. PTPRC exon-3 junction/read-level usage remains the more defensible truncation-robust measurement.

## Provenance

- Model analysis script: `run_ptprc_model_truncation_sensitivity_20260913.py` (in this directory)
- Summary script: `summarize_ptprc_model_truncation_sensitivity_20260913.py` (internal provenance; the shipped tables are its frozen output)
- PolyASite 2.0: human cluster atlas (external reference, not redistributed)
- Exact outputs: `ptprc_model_rows.tsv`, `ptprc_model_usage.tsv`, `ptprc_model_usage_by_run.tsv`, `ptprc_model_comparison.tsv`, `ptprc_endpoint_burden_by_filter.tsv`, `ptprc_model_sensitivity_summary.tsv`, `ptprc_model_bootstrap_summary.json`.
- The 1 missing run (`SRR30639550`) is recorded in `ptprc_model_missing.tsv`; it is not imputed.

This follows the field-standard logic of SQANTI3 structural/end QC and Bambu full-length-versus-partial read evidence, but the present implementation is deliberately a scoped PTPRC stress test rather than a claim of full-length truth.
- `ptprc_model_validation_report.md`, `ptprc_model_validation.tsv`, `ptprc_direction_validation.tsv`: independent validation of nested filters, paired shifts, BH correction and direction CIs.


## Reproduction in this repository

`run_ptprc_model_truncation_sensitivity_20260913.py` is the provenance
runner: it needs the run classification table, the run-local IsoQuant GTF
models, the human PolyASite 2.0 cluster atlas and the GRCh38 GENCODE GTF,
none of which are redistributed here (`SCTHREAD_PROJECT_ROOT` and the
documented external inputs locate them on the original machine layout).
`validate_ptprc_model_truncation_sensitivity_20260913.py` re-checks the
frozen outputs shipped beside it; its report is `ptprc_model_validation_report.md`.
