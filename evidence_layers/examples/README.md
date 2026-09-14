# Real output examples — and where each result appears

Small, real output fixtures for every analysis family the manuscript
describes. They are excerpts of frozen outputs, not simulated data, and each
maps to the script family that produced it and to the manuscript surface
where the result is reported.

| Example family | File | Produced by | Appears in |
|---|---|---|---|
| Transcript usage (per gene, per cell type) | `PTPRC_isoform_output.json` | `reducers/diu_reduce.py` + `reducers/diu_test.py` | Figure 3 panel a; `source_data/f2_grammar_figdata/diu_celltype.tsv` |
| Poly(A)-site usage (differential test rows) | `apa_site_output_example.tsv` | `reducers/apa_reduce.py` + `reducers/apa_test.py` | Figure 2 panel f; `source_data/biological_unit_validation/apa_observed_9999.tsv` |
| Allele-aware counts (gene-level imbalance) | `ase_summary.json` | `reducers/ase_celltype_reduce.py` + `reducers/ase_interaction_test.py` | Figure 2 panel f; `source_data/tables/ase_allelic_imbalance_feature_20260802.tsv` |
| Junction evidence (novel-junction replication) | `novel_junction_output.tsv` | `reducers/f2_junction_celltype.py` | SF2; `source_data/tables/junction_first100_live.tsv` |
| Exact-junction validation (sign-flip test) | `source_data/fig4_application/Figure2_source_discovery_exact_test.json` | junction discovery chain | Fig4 application figure; SF5 |
| CD45 read-versus-UMI sensitivity | `../cd45_sensitivity/cd45_umi_sensitivity_20260826.tsv` | `../cd45_sensitivity/measure_cd45_umi_sensitivity.py` | Figure 3 panels c–e; SF7 |
| Run-local transcript models (NIC/NNC catalogue rows) | `run_local_transcript_models_example.tsv` | IsoQuant per-run models (upstream of the catalogue) | Novel-structure counts in `database/release_content/`; SF2 |
| PTPRC model-level truncation sensitivity | `../ptprc_model_sensitivity/ptprc_model_sensitivity_summary.tsv` | `../ptprc_model_sensitivity/run_ptprc_model_truncation_sensitivity_20260913.py` | Manuscript Discussion; response to the platform-heterogeneity comment |

The `apa_site_output_example.tsv` and `run_local_transcript_models_example.tsv`
rows are verbatim excerpts of the referenced frozen tables (top-ranked APA
rows; the first 25 catalogue rows of one run).
