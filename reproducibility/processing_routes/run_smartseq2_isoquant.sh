#!/usr/bin/env bash
# Long-read Smart-seq2 (scONT plate) per-cell isoform pipeline.
#
# WHY this differs from the 10x route (wf-single-cell -> IsoQuant):
#   Smart-seq2 is plate-based: one CELL = one WELL = one library = one SRR.
#   Reads are full-length cDNA with NO cell barcode and NO UMI. So there is
#   nothing for wf-single-cell to demultiplex. Instead each SRR is aligned on
#   its own, and IsoQuant groups counts by INPUT FILE (--read_group file_name),
#   i.e. one column per cell. No UMI dedup (Smart-seq2 has none) -> read counts.
#
# Flow:  per cell  minimap2 -ax splice (ONT cDNA) -> sorted+indexed BAM
#        once      IsoQuant --bam <all cells> --read_group file_name --fl_data
#                           --data_type nanopore  -> per-cell isoform matrix
#
# Usage:
#   run_smartseq2_isoquant.sh <GSE> <GENOME.fa> <GENEDB.db> [ALIGN_JOBS] [SMOKE_N]
#     ALIGN_JOBS : concurrent minimap2 cells (default 16)
#     SMOKE_N    : if set, only process the first N cells (smoke test)
#
# Environment:
#   SCTHREAD_SCLONG_ROOT : data root hosting data/fastq/<GSE> and the output
#                          tree (default is the internal HPC layout; see
#                          reproducibility/ENVIRONMENT.md to override)
set -uo pipefail
GSE="$1"; GENOME="$2"; GENEDB="$3"; JOBS="${4:-16}"; SMOKE="${5:-0}"
ROOT="${SCTHREAD_SCLONG_ROOT:-/gpfs/home/fuzc/Seq_Database/scLong}"
FQDIR="$ROOT/data/fastq/$GSE"
OUT="$ROOT/data/isoquant_smartseq2/$GSE"
BAMDIR="$OUT/bam"; LOGDIR="$OUT/logs"
MM2="$ROOT/../../anaconda3/envs/isoquant/bin/minimap2"; [ -x "$MM2" ] || MM2=minimap2
ST="$ROOT/../../anaconda3/envs/isoquant/bin/samtools";  [ -x "$ST" ] || ST=samtools
IQ="$ROOT/../../anaconda3/envs/isoquant/bin/isoquant";  [ -x "$IQ" ] || IQ=isoquant
ALIGN_THREADS=6
mkdir -p "$BAMDIR" "$LOGDIR"

# --- collect per-cell single-end fastqs (each = one cell) ---
mapfile -t FQS < <(ls "$FQDIR"/*_1.fastq.gz 2>/dev/null | sort)
[ "${#FQS[@]}" -eq 0 ] && mapfile -t FQS < <(ls "$FQDIR"/*.fastq.gz 2>/dev/null | grep -vE "_2\.fastq" | sort)
[ "$SMOKE" -gt 0 ] && FQS=("${FQS[@]:0:$SMOKE}")
echo "[$(date '+%F %T')] $GSE: ${#FQS[@]} cells; align_jobs=$JOBS smoke=$SMOKE"

# --- 1. per-cell spliced alignment (parallel, capped) ---
[ -f "$GENOME.fai" ] || "$ST" faidx "$GENOME"
align_one() {
  local fq="$1"; local cell; cell=$(basename "$fq"); cell="${cell%%_1.fastq.gz}"; cell="${cell%%.fastq.gz}"
  local bam="$BAMDIR/$cell.bam"
  [ -s "$bam" ] && [ -s "$bam.bai" ] && return 0
  "$MM2" -ax splice -uf --MD -t "$ALIGN_THREADS" "$GENOME" "$fq" 2>"$LOGDIR/$cell.mm2.log" \
    | "$ST" sort -@ 2 -o "$bam" - 2>>"$LOGDIR/$cell.mm2.log" \
    && "$ST" index "$bam"
}
export -f align_one; export MM2 ST BAMDIR LOGDIR GENOME ALIGN_THREADS
printf '%s\n' "${FQS[@]}" | xargs -P "$JOBS" -I{} bash -c 'align_one "$@"' _ {}

# --- 2. verify all cells aligned, build bam list ---
BAMS=(); miss=0
for fq in "${FQS[@]}"; do
  cell=$(basename "$fq"); cell="${cell%%_1.fastq.gz}"; cell="${cell%%.fastq.gz}"
  if [ -s "$BAMDIR/$cell.bam.bai" ]; then BAMS+=("$BAMDIR/$cell.bam"); else echo "  MISS align: $cell"; miss=$((miss+1)); fi
done
echo "[$(date '+%F %T')] aligned ${#BAMS[@]}/${#FQS[@]} cells (miss=$miss)"
[ "${#BAMS[@]}" -eq 0 ] && { echo "no BAMs, abort"; exit 1; }

# --- 3. IsoQuant: one run, group counts by input file (= per cell) ---
echo "[$(date '+%F %T')] IsoQuant $GSE over ${#BAMS[@]} cells (file_name grouping, FL, no UMI)"
"$IQ" --bam "${BAMS[@]}" --reference "$GENOME" --genedb "$GENEDB" \
  --data_type nanopore --read_group file_name --fl_data \
  --output "$OUT/isoquant" --prefix "$GSE" --threads 32 >"$LOGDIR/isoquant.$GSE.log" 2>&1 \
  && touch "$OUT/ISOQ_OK" && echo "[$(date '+%F %T')] DONE $GSE" \
  || { echo "[$(date '+%F %T')] IsoQuant FAILED $GSE (see $LOGDIR/isoquant.$GSE.log)"; exit 1; }
