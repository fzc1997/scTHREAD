#!/bin/bash
# Fluidigm C1 PacBio (非10x, 每SRR=1 cell) → minimap2 splice:hq → IsoQuant file_name. mouse GRCm39.
# SCTHREAD_EXTSSD_ROOT: root hosting isoquant_env/, ref/GRCm39, fluidigm_input/
#                      (default is the internal workstation layout; override
#                       to point at your own tool env, reference and FASTQs)
set -eo pipefail
GSE=$1; shift
X="${SCTHREAD_EXTSSD_ROOT:-/mnt/extssd2}"
E=$X/isoquant_env/bin; REF=$X/ref/GRCm39
W=$X/fluidigm/$GSE; rm -rf "$W"; mkdir -p "$W/bam"; cd "$W"
for s in "$@"; do
  [ -f bam/$s.bam ] && continue
  $E/minimap2 -t 24 -ax splice:hq -uf --MD $REF/genome.fa $X/fluidigm_input/${s}.fastq.gz 2>/dev/null | $E/samtools sort -@8 -o bam/$s.bam -
  $E/samtools index bam/$s.bam
done
$E/isoquant --bam bam/*.bam --reference $REF/genome.fa --genedb $REF/genes.db \
  --data_type pacbio_ccs --read_group file_name --fl_data --threads 32 --prefix $GSE -o out 2>&1 | tail -3
echo "FLUIDIGM_DONE_$GSE $(head -1 out/$GSE/$GSE.gene_grouped_file_name_counts.tsv 2>/dev/null|cut -f1)"
