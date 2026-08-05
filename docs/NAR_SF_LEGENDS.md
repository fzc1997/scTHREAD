# scTHREAD — Supplementary Figure legends

**Outputs**: `NAR_database/figures/NAR_SF{1..4}.{pdf,png,svg}`
**Script**: `figures/scripts/render_nar_sf_compact.py`
**Panel provenance**: `docs/NAR_SF_COMPACT_PANEL_MAP.md`
**Style**: Nature-compatible (Arial, 8 pt upright panel labels, transparent PDF).
**Tone**: supporting / QC — **not** primary discovery claims.

**Denominators (use consistently):**

- **Frozen manuscript release**: 453 run records · 34 single-cell long-read
  datasets · **923,389** cells. Cells are authoritative per study, never summed
  over run rows.
- **Layer-audit scope**: the 15 studies carrying a differential-usage or
  allele-aware layer. This is a property of what those studies can support, not
  a processing tier — every dataset reaches the same isoform-level end product.
- **Published mouse cohort**: NGDC accession CRA044500 is **inside** the frozen
  release — six ONT runs of the 453, contributing 68,417 of the 923,389 cells.
  Its harmonized portal embedding holds 55,729 of those cells; the two numbers
  are different views of the same cohort and must not be swapped.
  `OWN_ASE_scONT` is only its internal scTHREAD dataset identifier.

> Superseded on 2026-08-04: the previous ten-figure pack (`NAR_SF1`–`NAR_SF10`,
> renderer `render_nar_sf_all.py`) is archived under
> `figures/superseded_sf_10pack_20260804/`. Five of its panels duplicated main
> Figure 2 and two of its figures are now carried by Supplementary Tables S5 and
> S7; see the panel map for the panel-by-panel disposition.

---

## Supplementary Figure 1 | Composition of the release and evidence-layer coverage

**(a)** Run records by species and platform: 86 human and 313 mouse Oxford
Nanopore records, 31 human and 23 mouse PacBio records.
**(b)** Cells by species. Human data contribute 736,853 cells across 117
records, mouse data 186,536 cells across 336 records; records and cells do not
track each other, because a plate-based run contributes one cell whereas a
droplet run contributes thousands.
**(c)** Cells by biological system, with the number of contributing datasets
against each bar.
**(d)** Barcode candidates retained as cells, by platform. Nanopore records
enter cell calling already resolved, whereas PacBio libraries carry many one- to
two-molecule ambient barcodes, so a much smaller share of candidates becomes
cells. Rendered from the live run manifest, so this ratio moves as records are
reprocessed; it is a QC field and is not a release denominator.
**(e)** Runs carrying each evidence layer, for the 15 studies with a
differential-usage or allele-aware layer. Numbers are runs per layer; each study's release cell count is printed to the
right. Two studies have no allele-aware run and
one has no differential-isoform-usage run, so a layer that is unavailable for a
record is reported as unavailable rather than as zero.
**(f)** Cells by cell-type annotation status, over all 34 datasets. 757,531
cells in 29 datasets carry a cell-type assignment; 165,858 cells in five do not.
Four of those five are a cell-line benchmark (GSE303762, 71,800 cells), a
perturbation screen (GSE295932, 21,623), a cancer cell line (GSE289790, 4,015)
and a method benchmark (GSE289428, 3), for which cell type is not defined. The
fifth is the mouse gastrulation cohort (`OWN_ASE_scONT`, 68,417), whose cell
types are carried by the separately maintained portal embedding shown in
Supplementary Figure 4a rather than by its release records. Annotation status is
recorded per run and is all-or-nothing within a dataset. This replaces a
per-study labelled fraction drawn from
`results/paper1/figdata/ed1_labeled_frac_bystudy.tsv`, which the release
contradicts in both directions.

Footer: the release line `453 run records · 34 datasets · 923,389 cells`.

---

## Supplementary Figure 2 | Splice-class composition and technical controls for novel structures

**(a)** Full-splice-match, incomplete-splice-match, novel-in-catalog and
novel-not-in-catalog composition of spliced molecules, by tissue.
**(b)** Novel fraction (novel-in-catalog plus novel-not-in-catalog) for each of
the 13 studies with full classification, ranging from 41.1% to 73.3%.
**(c)** Novel-isoform fraction against molecules per cell, within cell type and
at run level. Novelty tracks depth only weakly, so the novel fraction is not a
restatement of how deeply a record was sequenced.
**(d)** Depth-matched rarefaction of recurrent novel paths by cell type.
Recurrent structures survive subsampling, so recurrence is not produced by the
deepest records alone.
**(e)** Nested versus conflicting assignments among discordant molecules, for
the four studies with an assignment audit. Conflict rates run from 19.8% to
46.6%. **Method-supporting only** — a technical boundary of isoform assignment,
not a biological discovery claim.

Panels c and d support the depth statements in Results; the cross-study
canonicity axis they are usually read against is main Figure 2d.

---

## Supplementary Figure 3 | Registration-free query path

Six numbered screenshots of the portal, using **PTPRC** as the worked example:
home and catalog, gene search, the multi-layer gene card, per-cell-type isoform
usage, the download and export page, and the published API documentation. The
gene card reports differential isoform usage (effect 0.338, *q* = 0.000879),
differential poly(A)-site usage (effect 0.204, *q* = 0.000601) and a
non-significant allelic interaction (effect 0.194, *q* = 0.949), matching the
values reported in the main text and in Figure 3.

The junction counts on the screenshot do not: the tile shows 320 junctions in
the filtered view against 9,019 database-wide, whereas the main text and Figure
3 report 317 and 8,994. Both are correct for what they describe. The text and
Figure 3 are pinned to the dated API snapshot
`tables/PTPRC_api_overview_corrected_9999_20260727.json`, which is what makes
them reproducible; the screenshot shows the rolling live catalog on its capture
date, which has since grown. Effect sizes and *q*-values are unaffected, because
they come from the frozen analysis rather than from the catalog.

The figure demonstrates that multi-layer evidence is queryable and exportable
without registration and without local reanalysis. Tiles are laid out at a
common per-row height with widths set by each screenshot's own aspect ratio, so
the cropped detail views stay legible.

**Capture provenance**: all six shots were taken on 2026-08-05 against the
public site at device scale factor 2. The four full-page shots (1, 2, 5, 6) come
from `figures/scripts/capture_portal_walkthrough.py` at 1440x900 CSS px; the two
gene-card crops (3, 4) from `figures/scripts/capture_ptprc_live.py` at
1500x1200. The earlier pack was captured
against an internal preview build whose Download page still advertised the
superseded candidate release; both the home page and the Download page now lead
with the frozen release — 453 run records, 34 datasets, 923,389 cells — and
label the rolling live catalog separately. Shots 3 and 4 are cropped sub-views
of the PTPRC gene card and carry no release denominator.

Live site: https://scthread.ai4sc.ac.cn

---

## Supplementary Figure 4 | Cross-species worked example using a published mouse gastrulation cohort

The previously published CRA044500 cohort contributes six Oxford Nanopore runs
spanning embryonic days E6.5, E7.5 and E8.5. The portal serves them as one cell
map of 55,729 cells and 36 cell types, built from the forward and reciprocal
crosses with runs kept separately labelled rather than pooled as replicates; the
published ANCHOR differential-usage analysis of the same cohort uses a separate
25,621-cell, 24-cell-type forward-cross-only embedding, reported on that
denominator in Supplementary Table S8. The two denominators are never mixed.
The cohort is part of the frozen release: its six runs are among the 453 and
its 68,417 registered cells are among the 923,389. The 55,729 figure is the
subset that entered the harmonized embedding these panels are drawn on.

**(a)** Cell-type composition of the portal embedding, twelve largest of 36
cell types.
**(b,c)** The two displayed *Malat1* transcripts, `ENSMUST00000245150`
(12,433 cells, 16,979 molecules) and `ENSMUST00000172812` (9,810 cells,
12,920 molecules), on the shared portal coordinates. Cells without signal are
drawn in grey; colour is log(1 + molecules) capped at the 99th percentile. The
maps localize isoform expression and are not the differential-usage test.
**(d)** Within-gene usage of the same two transcripts by cell type, showing the
six most strongly biased cell types in each direction (maximum absolute
difference 0.213, anterior primitive streak versus neural crest).

**Evidence boundaries.** Usage proportions are descriptive and come from the
tables the portal serves (`tables/sf10_v3/`, built by
`scripts/build_sf10_v3_sources.py`), so the figure cannot drift away from the
live site. 30,931 of the 55,729 cells carry *Malat1* signal and usage fractions
use the full gene denominator. The shorter transcript is nested within the
longer one, so assignment between them depends on read length and 3' coverage.
The available cell-bootstrap pseudo-replicates are not independent embryos, so
no biological-replicate *P* value or false discovery rate is reported for this
cohort.

---

## Notes

1. Full supplementary figure pack for submission: **SF1–SF4**. Every figure is
   cited, and first mention runs S1 → S2 → S3 → S4.
2. Do not cite 4M / ~3M cells as the release size; use **923,389** for the
   frozen release, and label any rolling-registry number as such.
3. Main figures: Fig. 1 content and access · Fig. 2 evidence characterization ·
   Fig. 3 PTPRC utility · Fig. 4 monocyte-versus-T-cell junction usage.
4. The ASE / DIU / APA inventory and its label-permutation nulls are main
   Figure 2f and its footer, backed by Supplementary Table S5. The MS4A1 second
   gene case is Supplementary Table S7. Neither has a supplementary figure any
   more, by design.
