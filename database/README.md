# Part 1a — Building the database from raw reads

This directory documents how the frozen release was produced from the public
sequencing archives. Raw reads are not redistributed; every run is identified
by accession in the frozen release manifest.

## From FASTQ to the frozen release

Input runs are first split by how cell identity is recovered — the four
processing routes (113 + 42 + 286 + 12 = 453):

| Route | Runs | Software | Entry point |
|---|---|---|---|
| ONT droplet, barcode-resolved | 113 | wf-single-cell 3.3.4 (EPI2ME; minimap2 2.24-r1122 inside), barcode quality ≥ 15 | third-party workflow; per-run assignment in the manifest |
| PacBio, barcode-resolved | 42 (38 droplet-based + 4 pooled barcoded) | custom adapter/UMI extraction with whitelist cell-barcode correction → minimap2 2.31-r1302 `splice:hq` → CB/UB BAM tags | `processing_routes/` |
| ONT, sample/file-resolved (plate Smart-seq2) | 286 | minimap2 2.31-r1302 `splice`; one FASTQ = one cell | `processing_routes/run_smartseq2_isoquant.sh` |
| PacBio, sample/file-resolved (Fluidigm C1) | 12 | minimap2 2.31-r1302 `splice:hq`; 12 files = 6 objects × 2 chemistries | `processing_routes/run_fluidigm_c1_mouse.sh` |

All four routes quantify transcripts and isoforms with **IsoQuant 3.13.1** (alignment with minimap2 2.31-r1302 and SAMtools 1.24 where applied directly; the droplet route runs wf-single-cell 3.3.4 on Nextflow 26.04.4)
(`--read_group tag:CB` for barcode-resolved routes, `--read_group file_name
--fl_data` for sample/file-resolved routes). Shared output: gene counts,
transcript counts, read assignments and transcript models.

## Contents

| Path | Contents |
|---|---|
| `processing_routes/` | The five entry scripts for the three in-house routes (the PacBio extractors, the SAM CB/UB tagger, the Smart-seq2 and Fluidigm C1 runners), with per-route documentation and the exact differences from the in-house originals |
| `run_manifest/` | `processing_manifest_public.tsv` — the redacted per-run layer contract (route, kit, whitelist, output location for all 453 runs) — and the field-level schema |
| `frozen_release/` | The frozen release tables: `release_manifest.tsv` (every source accession), `release_study_cells.tsv` (per-study cell authority), `release_summary.json` |
| `release_content/` | Molecular inventory summaries by stratum for the frozen release |

The evidence layers computed on top of the per-run IsoQuant outputs live in
[`../evidence_layers/`](../evidence_layers/).
