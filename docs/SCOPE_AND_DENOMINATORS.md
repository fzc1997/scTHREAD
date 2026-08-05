# scTHREAD — scope and denominators

This file ships with the code deposit. It states what the frozen release
contains and which number to use for which question. It replaces the internal
working note `DATABASE_SCALE.md`, which is not part of the deposit.

## The release

`tables/release_20260803/` is the frozen snapshot the manuscript reports.

| | value |
|---|---|
| run records | **453** |
| datasets | **34** |
| cells | **923,389** |
| release identifier | `scTHREAD_NAR_release_20260803` |

These three numbers are the only release denominator. They also serve as the
frozen download at
`https://scthread.ai4sc.ac.cn/api/download/table/release_manifest?format=tsv`,
and the checksum of that download is recorded in the `Provenance` sheet of the
supplementary workbook.

## Scope rule

A dataset qualifies on a property of the **data** — single-cell resolution and
long-read sequencing — not on the archive it came from. Records reach scTHREAD
from the Gene Expression Omnibus, the European Nucleotide Archive, the Genome
Sequence Archive and direct submission, and the rule is applied identically to
all of them.

Every dataset is carried to the same isoform-level end product. The four
assay-aware processing routes differ only in how an assay encodes cell identity;
they are not quality or completeness tiers.

## Cells are per study, never summed over runs

`release_study_cells.tsv` is the authority, one row per dataset. Summing the
per-run `isoquant_cells` column of `release_manifest.tsv` gives **1,045,231**,
which is wrong: fourteen datasets carry an author-supplied or STARsolo total
that has no per-run decomposition, so a per-run sum both understates and
overstates individual studies.

Anything that reports a cell count must join `release_study_cells.tsv`.

## Two views of the same data

The portal's Search and Browse pages serve a **rolling registry** that grows as
records are added; the Download page serves the **frozen release** above. Both
are correct for their own question and must never be mixed in one sentence. A
figure or table that shows a rolling number states so on its face.

The published mouse gastrulation cohort NGDC `CRA044500` (internal identifier
`OWN_ASE_scONT`) is inside the frozen release: six of the 453 runs, 68,417 of
the 923,389 cells. Its harmonized portal embedding holds 55,729 of those cells,
and the previously published allele-aware analysis used a 25,621-cell
forward-cross-only embedding. The three figures describe the same cohort at
different stages and are not interchangeable.

## Availability of evidence layers

A record carries the gene, isoform, junction, poly(A) and allele layers its
reads support. Fifteen datasets carry a differential-usage or allele-aware
layer. Where a layer is unavailable for a record it is reported as unavailable,
never as a measured zero, and never as a processing tier.

Cell-type assignments follow the same principle: cells carry one wherever cell
type is defined for the dataset. Cell-line, perturbation and method-benchmark
datasets carry none, because cell type does not apply to them.

## Superseded numbers

If you meet any of these, the artefact producing them predates the frozen
release and must not be used:

| superseded | replaced by |
|---|---|
| 469 runs / 31 or 30 studies / 845,781 cells | 453 / 34 / 923,389 |
| 434 runs / 850,938 cells | 453 / 923,389 |
| `release_candidate_20260801`, `candidate_20260727` | `release_20260803` |
| gene universes 10,494 (isoform) and 13,214 (poly(A) site) | 8,092 and 10,531, from the biological-unit rerun |
| ten supplementary figures `NAR_SF1`–`NAR_SF10` | four, `NAR_SF1`–`NAR_SF4` |
