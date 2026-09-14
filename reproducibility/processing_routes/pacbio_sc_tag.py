#!/usr/bin/env python3
# Read minimap2 BAM whose read names carry |CB:xxx|UB:yyy; set CB/UB tags, strip suffix from name.
# usage: pacbio_sc_tag.py <in.bam> <out.bam>
import pysam, sys
inb, outb = sys.argv[1:3]
bf = pysam.AlignmentFile(inb, "rb")
of = pysam.AlignmentFile(outb, "wb", template=bf)
n = 0
for r in bf:
    nm = r.query_name
    if "|CB:" in nm:
        base, rest = nm.split("|CB:", 1)
        cb, ub = (rest.split("|UB:", 1) + [""])[:2]
        r.query_name = base
        r.set_tag("CB", cb, "Z")
        if ub: r.set_tag("UB", ub, "Z")
    of.write(r); n += 1
bf.close(); of.close()
sys.stderr.write(f"tagged {n} alignments\n")
