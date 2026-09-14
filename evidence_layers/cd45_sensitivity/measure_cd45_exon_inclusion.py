#!/usr/bin/env python3
"""Measure CD45 variable-exon (A/B/C) usage per cell type at the junction level.

Why not at the isoform level
----------------------------
The two annotated models that differ by exactly the three variable exons
(ENST00000442510 retains A/B/C, ENST00000348564 lacks all three) are seen in
only 1.7% and 1.0% of the portal atlas cells, so an isoform-resolved view of the
CD45RA/RO switch is under-powered. Worse, PTPRC isoform assignment as a whole is
confounded by 3' truncation: only 17-27% of long reads over PTPRC reach the
canonical 3' end, and read 3' ends in the proximal zone carry no poly(A) support
above background (6.2-7.7% within 50 bp of a PolyASite 2.0 cluster, against a
4.8% background), so proximal-terminating models act as truncation sinks.

The variable exons sit at chr1:198.696-198.703 Mb while the median read 3' end
is ~198.707 Mb, so most reads span the decision point regardless of where they
stop. Measuring the splice choice directly is therefore both better powered and
immune to the truncation artefact.

Definitions
-----------
informative read : has a splice junction whose donor is the constitutive exon-3
                   donor (chr1:198,692,373), i.e. the read demonstrably made the
                   include-or-skip choice.
RO               : that junction's acceptor is the constitutive exon-7 acceptor
                   (chr1:198,703,298), skipping all three variable exons.
otherwise        : classified by which of A/B/C the read's aligned blocks cover,
                   giving the high-molecular-weight side (RA/RB/RC).

Counts are reads, not UMI-collapsed molecules, and are reported per run so that
run-to-run reproducibility is visible rather than pooled away.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pysam


# GENCODE v44, GRCh38. PTPRC is on the + strand.
E3_DONOR = 198_692_373
E7_ACCEPTOR = 198_703_298
VARIABLE_EXONS = {
    "A": (198_696_712, 198_696_909),   # 198 bp
    "B": (198_699_564, 198_699_704),   # 141 bp
    "C": (198_702_387, 198_702_530),   # 144 bp
}
SCAN = ("chr1", 198_690_000, 198_706_000)
TOLERANCE = 10
MIN_EXON_OVERLAP = 20
MIN_INTRON = 20


def classify(read: pysam.AlignedSegment) -> str | None:
    """Return the variable-exon class of one read, or None if uninformative."""
    blocks = read.get_blocks()
    if len(blocks) < 2:
        return None
    junctions = [
        (end, start)
        for (_, end), (start, _) in zip(blocks, blocks[1:])
        if start - end > MIN_INTRON
    ]
    from_e3 = [j for j in junctions if abs(j[0] - E3_DONOR) <= TOLERANCE]
    if not from_e3:
        return None
    if any(abs(acceptor - E7_ACCEPTOR) <= TOLERANCE for _, acceptor in from_e3):
        return "RO"
    covered = [
        name
        for name, (start, end) in VARIABLE_EXONS.items()
        if any(
            min(end, b_end) - max(start, b_start) >= MIN_EXON_OVERLAP
            for b_start, b_end in blocks
        )
    ]
    return "".join(sorted(covered)) if covered else "other"


def count_run(bam: Path, barcodes: dict[str, str]) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    with pysam.AlignmentFile(bam) as handle:
        for read in handle.fetch(*SCAN):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            try:
                barcode = read.get_tag("CB")
            except KeyError:
                continue
            cell_type = barcodes.get(str(barcode).split("-")[0])
            if cell_type is None:
                continue
            label = classify(read)
            if label is None:
                continue
            key = (cell_type, label)
            counts[key] = counts.get(key, 0) + 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam-root", type=Path, required=True)
    parser.add_argument(
        "--barcode-map",
        type=Path,
        required=True,
        help="JSON: {run: {barcode: cell_type}} from the portal cell map",
    )
    parser.add_argument("--study", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    barcode_map = json.loads(args.barcode_map.read_text())
    # rglob can surface the same BAM twice through nested/symlinked work dirs;
    # counting a run twice would silently double its weight.
    bams = {}
    for path in sorted(args.bam_root.rglob("*.tagged.bam")):
        bams.setdefault(path.name.replace(".tagged.bam", ""), path)
    if not bams:
        raise SystemExit(f"No tagged BAMs under {args.bam_root}")
    bams = [bams[run] for run in sorted(bams)]

    rows = []
    for bam in bams:
        run = bam.name.replace(".tagged.bam", "")
        barcodes = barcode_map.get(run)
        if not barcodes:
            print(f"  skip {run}: no annotated barcodes in the cell map")
            continue
        counts = count_run(bam, barcodes)
        total = sum(counts.values())
        print(f"  {run}: {total:,} informative reads over {len(barcodes):,} annotated cells")
        for (cell_type, label), n in sorted(counts.items()):
            rows.append(
                {
                    "study": args.study,
                    "run": run,
                    "cell_type": cell_type,
                    "class": label,
                    "reads": n,
                }
            )
    if not rows:
        raise SystemExit("No informative reads found; check barcode-map run naming")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["study", "run", "cell_type", "class", "reads"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
