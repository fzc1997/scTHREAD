# Frozen source data

Every table and screenshot the figure renderers and validators read. Nothing
in this directory is recomputed; renderers fail closed when a file is missing,
and the frozen release manifest is checksum-pinned by the Figure 2 renderer.

| Path | Consumed by |
|---|---|
| `catalog_figdata/` | Figure 2 panels, SF1–SF2 (the frozen per-study figdata set behind the catalog, QC and splice-class panels) |
| `f2_grammar_figdata/` | Figure 2 panel f, Figure 3 panel a (DIU / APA / ASE cell-type results from the reducers in `evidence_layers/reducers/`) |
| `biological_unit_validation/` | Figure 2/3 statistics panels, SF5 QA (observed and 9999-permutation null tables plus the study-stratified replication outputs of the permutation framework) |
| `tables/` | Misc frozen tables: study atlas and biological-system map (SF1), ASE allelic-imbalance screen (Figure 2 panel f), CD45 exon-inclusion and PTPRC portal snapshots (Figure 3), the frozen sample-registry snapshot joined by the Figure 2 catalog loader |
| `portal_screenshots/` | SF3 walkthrough tiles and the Figure 3 portal views (captured 2026-08-05 / 2026-08-07 against the public site) |
| `sf4_malat1/` | SF4 cross-species worked example (MALAT1 isoform-pair composition, embedding/cell-type tables, UMAP coordinates, and the cohort annotation tables of the published mouse gastrulation study) |
| `sf5_discovery/` | SF5 nine-locus discovery screen (effects table and exact-test JSON) |
| `fig4_application/` | The junction-usage application figure: donor effects, cohort medians, cell-type profiles, state sensitivity, breadth boundary and the discovery-screen sources |
| `platform_sensitivity/` | Platform-sensitivity summaries for the retained MACF1 study (reference only) |
