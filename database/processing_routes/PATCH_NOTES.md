# Patch notes — public copies vs in-house sources

The three PacBio files are byte-identical copies of the in-house scripts used
for the frozen release (SHA-256):

```
6694c916573a14f951d696ce30f69808ca15f568ed4f3708ee18eae22e25ae2a  pacbio_sc_extract.py
5865949f8bc146099e5abab70b3f4497cf54efaae8f245dc09f0c3c0bc6ac399  pacbio_sc_extract_kinnex.py
1d41ec93fa4d3ada97ed1f1ca9eeaaac43c1a045b3178074d4b63fef5b0af8e6  pacbio_sc_tag.py
```

The two shell scripts differ from the in-house originals in exactly one
respect: hard-coded internal data roots are resolved through environment
variables (`SCTHREAD_SCLONG_ROOT`, `SCTHREAD_EXTSSD_ROOT`) whose defaults
reproduce the original machine layout, so behaviour at the default location is
unchanged. No flags, parameters, thresholds or pipeline steps were altered.

```diff
--- run_smartseq2_isoquant.sh (in-house)
+++ run_smartseq2_isoquant.sh (public)
@@ -16,9 +16,14 @@
 #   run_smartseq2_isoquant.sh <GSE> <GENOME.fa> <GENEDB.db> [ALIGN_JOBS] [SMOKE_N]
 #     ALIGN_JOBS : concurrent minimap2 cells (default 16)
 #     SMOKE_N    : if set, only process the first N cells (smoke test)
+#
+# Environment:
+#   SCTHREAD_SCLONG_ROOT : data root hosting data/fastq/<GSE> and the output
+#                          tree (default is the internal HPC layout; see
+#                          ENVIRONMENT.md (repository root) to override)
 set -uo pipefail
 GSE="$1"; GENOME="$2"; GENEDB="$3"; JOBS="${4:-16}"; SMOKE="${5:-0}"
-ROOT=/gpfs/home/fuzc/Seq_Database/scLong
+ROOT="${SCTHREAD_SCLONG_ROOT:-/gpfs/home/fuzc/Seq_Database/scLong}"
 FQDIR="$ROOT/data/fastq/$GSE"
 OUT="$ROOT/data/isoquant_smartseq2/$GSE"
 BAMDIR="$OUT/bam"; LOGDIR="$OUT/logs"
```

```diff
--- run_fluidigm_c1_mouse.sh (in-house)
+++ run_fluidigm_c1_mouse.sh (public)
@@ -1,12 +1,16 @@
 #!/bin/bash
 # Fluidigm C1 PacBio (非10x, 每SRR=1 cell) → minimap2 splice:hq → IsoQuant file_name. mouse GRCm39.
+# SCTHREAD_EXTSSD_ROOT: root hosting isoquant_env/, ref/GRCm39, fluidigm_input/
+#                      (default is the internal workstation layout; override
+#                       to point at your own tool env, reference and FASTQs)
 set -eo pipefail
 GSE=$1; shift
-E=/mnt/extssd2/isoquant_env/bin; REF=/mnt/extssd2/ref/GRCm39
-W=/mnt/extssd2/fluidigm/$GSE; rm -rf "$W"; mkdir -p "$W/bam"; cd "$W"
+X="${SCTHREAD_EXTSSD_ROOT:-/mnt/extssd2}"
+E=$X/isoquant_env/bin; REF=$X/ref/GRCm39
+W=$X/fluidigm/$GSE; rm -rf "$W"; mkdir -p "$W/bam"; cd "$W"
 for s in "$@"; do
   [ -f bam/$s.bam ] && continue
-  $E/minimap2 -t 24 -ax splice:hq -uf --MD $REF/genome.fa /mnt/extssd2/fluidigm_input/${s}.fastq.gz 2>/dev/null | $E/samtools sort -@8 -o bam/$s.bam -
+  $E/minimap2 -t 24 -ax splice:hq -uf --MD $REF/genome.fa $X/fluidigm_input/${s}.fastq.gz 2>/dev/null | $E/samtools sort -@8 -o bam/$s.bam -
   $E/samtools index bam/$s.bam
 done
 $E/isoquant --bam bam/*.bam --reference $REF/genome.fa --genedb $REF/genes.db \
```

In-house sources (authoritative, unchanged):

- `Seq_Database/scLong/pacbio_pilot/pacbio_sc_extract.py`
- `Seq_Database/scLong/pacbio_pilot/pacbio_sc_extract_kinnex.py`
- `Seq_Database/scLong/pacbio_pilot/pacbio_sc_tag.py`
- `Seq_Database/scLong/scripts/smartseq2_longread/run_smartseq2_isoquant.sh`
- `Seq_Database/scLong/scripts/smartseq2_longread/run_fluidigm_c1_mouse.sh`
