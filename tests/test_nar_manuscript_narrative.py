#!/usr/bin/env python3
"""Deterministic narrative and consistency checks for the NAR manuscript."""

from __future__ import annotations

import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT / "manuscript/scTHREAD_NAR_Database_Issue_draft.md"
SUPPLEMENTARY_FIGURES = 4  # the compacted pack; see docs/NAR_SF_COMPACT_PANEL_MAP.md


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def word_count(text: str) -> int:
    return len(re.findall(r"\b[\w′'-]+\b", text))


def cited_numbers(text: str) -> set[int]:
    numbers: set[int] = set()
    for block in re.findall(r"\[([0-9,–-]+)\]", text):
        for item in block.split(","):
            if "–" in item or "-" in item:
                lo, hi = re.split(r"[–-]", item)
                numbers.update(range(int(lo), int(hi) + 1))
            else:
                numbers.add(int(item))
    return numbers


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    normalized = re.sub(r"\s+", " ", text)
    title = text.splitlines()[0]
    require(title.startswith("# scTHREAD:"), "Title must begin with scTHREAD")
    require("database" in title.lower(), "Title must identify the database")

    abstract = text.split("## Abstract\n", 1)[1].split("\n## Introduction", 1)[0]
    require(150 <= word_count(abstract) <= 200, "Abstract must contain 150–200 words")
    require("https://scthread.ai4sc.ac.cn/" in abstract, "Abstract URL missing")

    ordered_headings = [
        "## Introduction",
        "## Materials and methods",
        "## Results",
        # Corpus convention is short noun-phrase headings, and exactly one
        # interface subsection - the walkthrough must not be told three times.
        "### Database content",
        "### Most stored reads carry structures absent from the reference",
        "### Evidence fields for unannotated structures",
        "### Differential-usage annotations are precomputed against a null",
        "### Web interface and data access",
        "### Case study: retrieving isoform-level evidence for one gene",
        "## Discussion",
        "## Data availability",
        "## Ethics statement",
    ]
    positions = [text.index(heading) for heading in ordered_headings]
    require(positions == sorted(positions), "NAR narrative section order changed")

    for phrase in (
        "453 run records from 34",
        # how 923,389 is derived: per dataset, with a stated precedence.
        # The tiers must add to the headline, and the run sum (1,045,231)
        # must never be presented as an alternative total.
        "Cell counts are authoritative per dataset rather than",
        "11 datasets, 104,161 cells",
        "4 datasets, 45,172",
        "17 datasets, 773,802",
        "923,389",
        # content characterization that must not silently disappear
        "44.8-45.8% of spliced reads",
        "8% of spliced reads in erythroid cells",
        "80% in malignant plasma cells",
        # the eligible universes are the biological-unit ones (8,092 / 10,531);
        # the run-as-donor tables that gave 10,494 / 13,214 are superseded
        "median gene carries seven distinct isoforms",
        "15 supported poly(A) sites",
        # what the catalog adds beyond the existing annotation
        "1.70 billion classified spliced",
        "28.9% match an annotated transcript end to end",
        "860,127 (36.3%) in two or more",
        "710,535",
        "median of three exons",
        "incomplete molecules produce the",
        # the three stated uses of the resource
        "designing an assay that distinguishes two isoforms of one gene",
        "checking whether a published isoform switch appears in other data",
        "benchmarking isoform-level methods",
        # PTPRC/MS4A1 are known-biology checks, not discoveries
        "expected pattern is fixed by prior work",
        "described in Supplementary Note 1",
        "117 records and 736,853 cells",
        "336 records and 186,536 cells",
        # repositories named in the paper, accessions in the supplement
        "Gene Expression Omnibus, the European Nucleotide Archive and the Genome Sequence Archive",
        # Methods must name the tools that produced the records, with versions
        "wf-single-cell v3.3.4",
        "IsoQuant v3.13.1",
        # the PacBio route recovers barcodes from the reads themselves; the
        # skera/lima/isoseq chain was never run and must not come back
        "CTACACGACGCTCTTCCGATCT",
        "minimap2 v2.31-r1302",
        "SAMtools v1.24",
        "wf-single-cell v3.3.4",
        # ANCHOR produced the allele-aware layer and must be credited as the
        # method, not only as the source of the reused mouse dataset
        "produced with ANCHOR",
        "calls variants de novo",
        "minimap2",
        # one uniform processing route, no data tiers
        "Every dataset is carried to the same isoform-level end product",
        "cells carry a cell-type assignment wherever cell type is defined",
        # one cell denominator per scope: the catalog, the labelled subset,
        # and the mouse portal map. Near-duplicate totals must not come back.
        "2,341,560",
        "no biological-replicate P-value or FDR is reported",
        "FastAPI",
        "DuckDB",
        "gene-level annotation",
        "CRA044500",
        "10.64898/2026.06.08.730656",
        "No new animal experiments were performed for this study",
        "12,471,531",
        "8,698,826",
        "210,341 distinct Ensembl",
        "Supplementary Tables S1–S15",
        "D291–D302",
        "10.1093/nar/gkaf1144",
        # snapshot-scoped composition promoted into Results (must match Fig.1)
        "2,369,193 unannotated junctions",
        "97.0% for the 16,983 seen in all 13",
        # main-text feature comparison replaces the supplementary-only claim
        # SF10 was rebuilt on the current portal embedding
        "55,729 cells and 36 cell types",
        "30,931 cells carry signal",
        # 0.213 is the gap between the two transcripts inside neural crest,
        # not a contrast between two cell types
        "largest gap between them is 0.213, in neural crest",
        # the mouse example is one of the 34 studies, not an extra cohort
        "its 68,417 cells are part of the 923,389 reported in the main text",
        "### Case study: a monocyte-versus-T-cell junction-usage difference",
        "### Cell-type junction-usage analysis",
        "13 of 13 donors",
        "across the nine loci",
        "chr6:41,940,586-41,941,451",
        "Five loci passed at *q* < 0.05",
        "positive in 4 of 4 donors",
        "the three-locus pan-myeloid model was not supported",
        "no shared regulation or mechanism is claimed",
        # the replicated junction program: numbers and the prospective negative
    ):
        require(phrase in normalized, f"Required manuscript fact missing: {phrase}")

    # Table 1 and the figure legends sit after the reference list, so anything
    # that rewrites the tail can silently truncate them.
    # the tail was duplicated once during a restore; a section must appear once
    for section in ("## Table 1", "## Figure legends", "## Supplementary Note 1"):
        require(text.count(section) == 1,
                f"{section} appears {text.count(section)} times, expected once")
    for section in ("## Table 1", "## Figure legends",
                    "### Figure 1.", "### Figure 2.", "### Figure 3.", "### Figure 4."):
        require(section in text, f"Section lost from the manuscript tail: {section}")

    for token in ("2,369,193", "860,127", "710,535"):
        require(
            normalized.count(token) == 1,
            f"{token} is reported more than once in the body text",
        )

    for phrase in (
        # claims that failed their own orthogonal tests upstream and must never
        # reach the paper: CUTA NMD (three independent no-support results, one
        # in the opposite direction) and the falsified reciprocal-isoform set
        "nonsense-mediated decay",
        # 17 runs across 5 datasets are no_celltype_by_design in the release
        "every cell carries a cell-type assignment",
        # run-as-donor gene universes, superseded by the biological-unit rerun
        "10,494",
        "13,214",
        # CRA044500 is one of the 34 studies: 6 runs, 68,417 of the 923,389
        # cells. Calling it outside the snapshot invents a second denominator.
        "additional to the counts in the main text",
        "outside the frozen manuscript snapshot",
        "escapes NMD",
        "junction-usage program",
        "junction-usage signature",
        "coordinated program",
        "reflects shared regulation",
        "demonstrates causal",
        "reciprocal isoform program",
        "S100A6",
        "EMP3",
        "469 IsoQuant-complete",
        "neither a long-read isoform collection",
        "not first collection",
        "author-provided LAN deployment",
        "frozen candidate manifest",
        "corrected analysis-status",
        "Cell-type-specific Malat1 DTU",
        "gene-level FDR",
        "The labelled release currently",
        "A gene was called only",
        "no FDR discoveries",
        "internally generated mouse",
        "internally produced mouse",
        "REPOSITORY AND ACCESSION TO BE SUPPLIED",
        "D291–D303",
        # live-registry counts must never leak back into the frozen snapshot text
        "941,733",
        "478 run",
        "33 studies",
        # pre-audit denominators must not survive anywhere in the text
        "850,938",
        "434 run",
        "30 studies",
        "469-record",
        "469-row",
        # pre-correction cell totals must not survive either
        "845,746",
        "845,448",
        "161,933",
        "135 of 136",
        # 850,640 differs from 850,938 by exactly the 298 one-well cells and
        # only confused readers; the tier is stated in records instead
        "850,640",
        # the ANCHOR embedding size belongs to Supplementary Table S8, not the main text
        "25,621-cell",
        # the comparison table is now main text, not supplementary-only
        "comparison is provided in Supplementary Table S10",
        # the paper must not describe the database as tiered or partially processed
        "membership metadata",
        "membership-only",
        "completeness tier",
        "transcript-feature-complete",
        "one-well records",
        "separately identified subset",
        "labelled subset",
        # and must not name the resource before introducing it
        "Among the resources compared, scTHREAD",
        # Internal identifiers are not reader-facing. A column name, a dataset key,
        # a command-line flag or a file name tells the reader nothing they can act
        # on, and their repeated appearance is what made earlier drafts read like
        # engineering documentation rather than a paper.
        "isoquant_cells",
        "isoquant_transcripts",
        "annotation_cells",
        "smartseq2_one_well",
        "study::source",
        "hapA",
        "hapB",
        "OWN_ASE_scONT",
        "IsoQuant-core-complete",
        "--data_type",
        "--read_group",
        "skera",
        "pbmm2",
        "isoseq groupdedup",
        ".tsv",
        ".gtf",
        ".xlsx",
        "/gpfs/",
        # a paper states inclusion criteria; it does not narrate the engineering audit
        "modality audit",
        "PRJNA591948",
        "LibrarySource=GENOMIC",
        "file-contract audit",
        "The file recount",
        "rather than from screenshots",
        "only such correction",
        # the superseded portal embedding must not reappear
        "8,367",
        "shared 25,621-cell embedding",
    ):
        require(phrase not in normalized, f"Forbidden legacy wording remains: {phrase}")

    require("## Table 1" in text, "Main-text Table 1 is missing")

    # Nothing that ships may carry CJK text. The release manifest is the file
    # Data availability offers for download, and it carried Chinese in two
    # columns for 333 of its 453 rows.
    cjk = re.compile(r"[㐀-鿿　-〿＀-￯]")
    for shipped in (
        MANUSCRIPT,
        PROJECT / "tables/release_20260803/release_manifest.tsv",
        PROJECT / "tables/release_20260803/release_study_cells.tsv",
        PROJECT / "docs/NAR_SF_LEGENDS.md",
    ):
        if not shipped.is_file():
            continue
        hits = [
            f"line {n}: {line.strip()[:60]}"
            for n, line in enumerate(shipped.read_text(encoding="utf-8").splitlines(), 1)
            if cjk.search(line)
        ]
        require(
            not hits,
            f"{shipped.name} ships CJK text ({len(hits)} lines): {hits[0] if hits else ''}",
        )

    # A whole Table 1 - caption, table and footnote - once sat between the last
    # reference and "## Figure legends", carrying no heading of its own, so a
    # heading count of 1 was satisfied while the docx printed the table twice.
    # Count the captions, not the headings.
    for caption in (
        "**Table 1. Feature comparison between scTHREAD",
        "### Figure 1. Content, construction and access routes",
    ):
        require(
            text.count(caption) == 1,
            f"{caption[:44]}... appears {text.count(caption)} times, expected once",
        )
    # and hold the reference list to nothing but references, which is the shape
    # every truncation and duplication bug here has broken
    refs_block = text.split("\n## References", 1)[1]
    refs_block = re.split(r"^## ", refs_block, maxsplit=1, flags=re.M)[0]
    stray = [
        line for line in refs_block.splitlines()
        if line.strip() and not re.match(r"^\d+\. ", line)
    ]
    require(
        not stray,
        "the reference list contains non-reference lines: "
        + "; ".join(line[:60] for line in stray[:3]),
    )

    # The supplementary pack was compacted from ten figures to four on
    # 2026-08-04, so a citation above S4 points at a figure that no longer
    # exists. Ranges count too: "S1-S10" is the wording that was in the
    # Supplementary data section. The draft hard-wraps, so match across newlines.
    stale = set()
    for match in re.finditer(
        r"Supplementary\s+Figures?\s+S(\d+)(?:\s*[-–—]\s*S?(\d+))?", text
    ):
        low = int(match.group(1))
        high = int(match.group(2)) if match.group(2) else low
        stale.update(n for n in range(low, high + 1) if n > SUPPLEMENTARY_FIGURES)
    require(
        not stale,
        "supplementary figures above S%d are cited but do not exist: %s"
        % (SUPPLEMENTARY_FIGURES, ", ".join(f"S{n}" for n in sorted(stale))),
    )
    present = sorted(
        int(p.stem.removeprefix("NAR_SF")) for p in (PROJECT / "figures").glob("NAR_SF*.pdf")
    )
    require(
        present == list(range(1, SUPPLEMENTARY_FIGURES + 1)),
        f"expected figures/NAR_SF1..SF{SUPPLEMENTARY_FIGURES}.pdf, found {present}",
    )

    main_text, references = text.split("## References\n", 1)
    cited = cited_numbers(text)  # main text plus the tables/legends printed after the list
    listed = {int(value) for value in re.findall(r"^(\d+)\.", references, re.MULTILINE)}
    require(len(listed) >= 30, f"corpus expects >=30 references, found {len(listed)}")
    require(
        cited == listed == set(range(1, len(listed) + 1)),
        "Citation/reference numbering mismatch",
    )
    dois = re.findall(
        r"doi:(10\.\d{4,9}/[^\s]+)\.$", references, re.I | re.MULTILINE
    )
    require(
        len(dois) == len(set(map(str.lower, dois))) == len(listed),
        f"every reference needs a unique DOI ({len(dois)} DOIs for {len(listed)} entries)",
    )
    # Cold Spring Harbor uses 10.1101 for both bioRxiv and its journals; bioRxiv
    # suffixes are numeric (10.1101/2023.05.17.541248, 10.1101/068809) whereas
    # journal suffixes start with a letter (10.1101/gr.222976.117).
    require(
        not [doi for doi in dois if re.match(r"10\.1101/\d", doi.lower())],
        "bioRxiv DOI used where a journal version exists",
    )
    for phrase in (
        "Zhu S, Lian Q, Ye W, et al.",
        "Moon Y, Herrmann CJ, Mironov A, Zavolan M.",
        "McInnes L, Healy J, Saul N, Großberger L.",
    ):
        require(phrase in references, f"Crossref-corrected reference missing: {phrase}")

    print(
        "NAR MANUSCRIPT NARRATIVE QA PASS",
        {
            "abstract_words": word_count(abstract),
            "references": len(listed),
            "main_sections": len(ordered_headings),
        },
    )


if __name__ == "__main__":
    main()
