# scTHREAD — database construction and figure reproduction package

scTHREAD is a single-cell long-read transcript database for human and mouse
(Oxford Nanopore and PacBio). The live database is freely accessible at
https://scthread.ai4sc.ac.cn/ (no login required). This repository contains the
code and frozen data needed to understand and reproduce two things:

1. **How the database was built** — from raw sequencing reads (FASTQ, retrieved
   from the public archives by accession) through four processing routes to
   the frozen release tables and the per-cell-type evidence layers.
   → [`database/`](database/) and [`evidence_layers/`](evidence_layers/)

2. **How every figure was made** — one folder per manuscript figure with its
   renderer and the exact frozen input tables.
   → [`figures/`](figures/) (inputs under [`source_data/`](source_data/))

## Repository map

| Directory | Contents |
|---|---|
| `database/` | Part 1a — construction: processing-route entry scripts, the per-run manifest, and the frozen release tables (453 runs, 34 datasets, 923,389 cells) |
| `evidence_layers/` | Part 1b — reduction of per-run IsoQuant output into the DIU / APA / allele-specific / junction evidence layers, the biological-unit permutation framework, worked examples, and the CD45 read-versus-UMI sensitivity analysis |
| `figures/` | Part 2 — every manuscript figure: main Figures 1–3, Supplementary Figures S1–S7, the graphical abstract; each with its renderer and QA script |
| `source_data/` | The frozen tables and screenshots the figure renderers and validators read; nothing here is recomputed |
| `workflows/` | One-command validation of the whole package |
| `ENVIRONMENT.md` | Environment variables used by the scripts |

## Quick start

```bash
# validate the package (compiles everything, validates outputs and manifest)
./workflows/run_public_validation.sh

# re-render any figure, e.g. main Figure 2 and Figure 3
python figures/main/fig2/render_nar_bio.py --fig 2 --stem2 /tmp/NAR_Fig2 --observed-suffix _9999
python figures/main/fig3/render_nar_fig3_v8.py --stem /tmp/NAR_Fig3 --observed-suffix _9999

# supplementary figures S1–S4 (one command) and S5
python figures/supplementary/sf1_sf4_sf7/render_nar_sf_compact.py --outdir /tmp/sf
python figures/supplementary/sf5_discovery/render_nar_sf5_discovery_screen.py \
  --effects source_data/sf5_discovery/Figure2_source_discovery_panel.tsv \
  --test-json source_data/sf5_discovery/Figure2_source_discovery_exact_test.json \
  --output-dir /tmp/sf5
```

Figure-specific commands and caveats: [`figures/README.md`](figures/README.md).

## Examples and where each result appears

One real output example for every analysis family, with the producing script
and the manuscript surface it feeds, is indexed in
[`evidence_layers/examples/README.md`](evidence_layers/examples/README.md).

## Associated manuscript

This package accompanies the manuscript *"scTHREAD: a cell-resolved
long-read transcriptome database with linked isoform, poly(A), junction and
allelic evidence"* (Fu Z-C, et al.; under review).

## Frozen release scope

The release reported in the manuscript is **453 run records from 34 datasets
(923,389 cells)**: 399 Oxford Nanopore and 54 PacBio records; 155
barcode-resolved and 298 sample/file-resolved. Renderers fail closed if a
frozen input is missing or a checksum changes; they never re-derive catalog
numbers from rolling sources.

## License

CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
