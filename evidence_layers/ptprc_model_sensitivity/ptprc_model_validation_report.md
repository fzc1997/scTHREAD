# PTPRC model truncation sensitivity validation

Status: **PASS** for design and statistical validation of the scoped sensitivity output.

| Filter | Metric | n | Wilcoxon statistic | P | BH q | Median delta | Positive fraction |
|---|---:|---:|---:|---:|---:|---:|---:|
| F1_pas_supported_models | RA_like_usage_delta | 33 | 0 | 7.95e-07 | 1.59e-06 | 0.170 | 0.970 |
| F1_pas_supported_models | RO_like_usage_delta | 33 | 0 | 0.000982 | 0.000982 | 0.000 | 0.424 |
| F2_high_confidence_reference_models | RA_like_usage_delta | 33 | 0 | 7.95e-07 | 1.59e-06 | 0.577 | 0.970 |
| F2_high_confidence_reference_models | RO_like_usage_delta | 33 | 0 | 0.000982 | 0.000982 | 0.000 | 0.424 |

Assumption checks: paired run-level units are fixed by the same run IDs; no cross-study pooling is introduced. The primary effect is median paired delta with the pre-registered 1,000-resample bootstrap CI. BH correction spans the four pre-specified paired tests.
The validated result is a scoped endpoint/completeness sensitivity signal. It is not a claim that all excluded models are technical artifacts, and it does not replace the direct PTPRC exon-3 junction measurement.
