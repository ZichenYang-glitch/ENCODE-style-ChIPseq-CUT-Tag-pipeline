#!/usr/bin/env bash
set -euo pipefail

############################################
# ENCODE-style ChIP/CUT&Tag pipeline (no input)
# Author: yangzichen-glitch
# Usage examples:
#   Single-end:
#     ./chipseq.sh -s SAMPLE -x /path/to/bowtie2/index -g mm --se -r1 sample.fastq.gz
#
#   Paired-end:
#     ./chipseq.sh -s SAMPLE -x /path/to/bowtie2/index -g mm --pe -r1 R1.fastq.gz -r2 R2.fastq.gz
#
# Notes:
# - Default: NO input/control (works for CUT&Tag/CUT&RUN, and acceptable for some ChIP-seq cases)
# - You can optionally provide --control_bam CONTROL.bam to do with-control MACS2.
############################################

# ---------- Defaults ----------
THREADS=8
GENOME="mm"          # mm (mouse) or hs (human)
MAPQ=30
OUTDIR="./chip_out"
BINSIZE=10
REMOVE_DUP="auto"    # auto | yes | no
PEAK_MODE="narrow"   # narrow | broad
TECH="chip"          # chip | cuttag (Switch for MACS3 parameters)
EXTEND_READS="auto"  # auto | yes | no
TRIM="yes"           # yes | no

# Input/control
SAMPLE=""
BOWTIE2_INDEX=""
SE_MODE="no"
PE_MODE="no"
R1=""
R2=""
CONTROL_BAM=""

# ---------- Helpers ----------
# Define log function: prints messages with a timestamp
log() {
    echo -e "[$(date '+%F %T')] $*"
}

# Define die function: prints an error message and exits the script
die() {
    echo -e "[$(date '+%F %T')] [ERROR] $*" >&2
    exit 1
}

# ---------- Parse args ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample) SAMPLE="$2"; shift 2;;
    -x|--index) BOWTIE2_INDEX="$2"; shift 2;;
    -g|--genome) GENOME="$2"; shift 2;;
    -p|--threads) THREADS="$2"; shift 2;;
    --mapq) MAPQ="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    --binSize) BINSIZE="$2"; shift 2;;
    --trim) TRIM="yes"; shift 1;;
    --tech) TECH="$2"; shift 2;; # Added: chip or cuttag

    --se) SE_MODE="yes"; PE_MODE="no"; shift 1;;
    --pe) PE_MODE="yes"; SE_MODE="no"; shift 1;;
    -r1) R1="$2"; shift 2;;
    -r2) R2="$2"; shift 2;;

    --removeDup) REMOVE_DUP="$2"; shift 2;;
    --peakMode) PEAK_MODE="$2"; shift 2;;
    --extendReads) EXTEND_READS="$2"; shift 2;;

    --control_bam) CONTROL_BAM="$2"; shift 2;;
    -h|--help)
      sed -n '1,120p' "$0"
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

# ---------- Validate ----------
[[ -n "$SAMPLE" ]] || die "Missing --sample"
[[ -n "$BOWTIE2_INDEX" ]] || die "Missing --index (bowtie2 index prefix)"
[[ "$SE_MODE" == "yes" || "$PE_MODE" == "yes" ]] || die "Choose one: --se or --pe"
[[ -n "$R1" ]] || die "Missing -r1 FASTQ"
if [[ "$PE_MODE" == "yes" ]]; then
  [[ -n "$R2" ]] || die "Paired-end needs -r2 FASTQ"
fi
[[ "$PEAK_MODE" == "narrow" || "$PEAK_MODE" == "broad" ]] || die "--peakMode must be narrow or broad"

# ---------- Genome size for MACS2 ----------
# MACS3 -g: Accepts shortcuts like 'mm' or 'hs', or a specific number (e.g., 1.87e9).
# This mapping ensures compatibility between common genome build names and MACS shortcuts.
MACS_G="$GENOME"
if [[ "$GENOME" == "mm10" || "$GENOME" == "mm39" ]]; then
    MACS_G="mm"
elif [[ "$GENOME" == "hg19" || "$GENOME" == "hg38" ]]; then
    MACS_G="hs"
fi

log "Genome size shortcut used: $MACS_G"

# ---------- Directories ----------
RAWDIR="${OUTDIR}/00_raw"
QCDIR="${OUTDIR}/01_qc"
ALNDIR="${OUTDIR}/02_align"
BWDIR="${OUTDIR}/03_bigwig"
PEAKDIR="${OUTDIR}/04_peaks"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$RAWDIR" "$QCDIR" "$ALNDIR" "$BWDIR" "$PEAKDIR" "$LOGDIR"

log "Sample            : $SAMPLE"
log "Genome (-g)       : $GENOME (MACS uses: $MACS_G)"  
log "Threads          : $THREADS"
log "Mode             : $([[ "$PE_MODE" == "yes" ]] && echo PE || echo SE)"
log "Peak mode        : $PEAK_MODE"
log "Protocol (Tech)  : $TECH"                          
log "Outdir           : $OUTDIR"
log "Control BAM      : ${CONTROL_BAM:-NONE (no-input)}"

# ---------- Step 1: FastQC (optional but nice) ----------
if command -v fastqc >/dev/null 2>&1; then
  log "FastQC..."
  fastqc -t "$THREADS" -o "$QCDIR" "$R1" ${R2:+$R2} |& tee "${LOGDIR}/${SAMPLE}.fastqc.log" || true
else
  log "FastQC not found, skip."
fi

# ---------- Step 2: Trim (optional) ----------
TRIM_R1="$R1"
TRIM_R2="$R2"
if [[ "$TRIM" == "yes" ]]; then
  command -v trim_galore >/dev/null 2>&1 || die "trim_galore not found but --trim set"
  log "Trimming with Trim Galore..."
  if [[ "$SE_MODE" == "yes" ]]; then
    trim_galore --cores "$THREADS" -o "$RAWDIR" "$R1" |& tee "${LOGDIR}/${SAMPLE}.trim.log"
    TRIM_R1="$(ls -1 "$RAWDIR"/*_trimmed.fq.gz | tail -n 1)"
  else
    trim_galore --paired --cores "$THREADS" -o "$RAWDIR" "$R1" "$R2" |& tee "${LOGDIR}/${SAMPLE}.trim.log"
    TRIM_R1="$(ls -1 "$RAWDIR"/*_val_1.fq.gz | tail -n 1)"
    TRIM_R2="$(ls -1 "$RAWDIR"/*_val_2.fq.gz | tail -n 1)"
  fi
  log "Trimmed R1: $TRIM_R1"
  [[ "$PE_MODE" == "yes" ]] && log "Trimmed R2: $TRIM_R2"
fi

# ---------- Step 3: Align (Bowtie2) ----------
# ---------- Step 3: Align (Bowtie2) ----------
log "Aligning reads with Bowtie2..."
BAM_SORT="${ALNDIR}/${SAMPLE}.sorted.bam"

# Define Read Group (RG) information to ensure compatibility with Picard MarkDuplicates
# RG ID: Unique ID, SM: Sample Name, LB: Library, PL: Platform
RG_FLAGS="--rg-id ${SAMPLE} --rg SM:${SAMPLE} --rg LB:${SAMPLE} --rg PL:ILLUMINA"

if [[ "$SE_MODE" == "yes" ]]; then
    # Single-End (SE) Alignment mode
    bowtie2 $RG_FLAGS -x "$BOWTIE2_INDEX" -U "$TRIM_R1" -p "$THREADS" 2> "${LOGDIR}/${SAMPLE}.bowtie2.log" \
      | samtools view -b - \
      | samtools sort -@ "$THREADS" -o "$BAM_SORT" -
else
    # Paired-End (PE) Alignment mode
    bowtie2 $RG_FLAGS -x "$BOWTIE2_INDEX" -1 "$TRIM_R1" -2 "$TRIM_R2" -p "$THREADS" 2> "${LOGDIR}/${SAMPLE}.bowtie2.log" \
      | samtools view -b - \
      | samtools sort -@ "$THREADS" -o "$BAM_SORT" -
fi

# Index the sorted BAM file for fast random access
samtools index "$BAM_SORT"

# Generate basic alignment statistics for Quality Control
log "Generating alignment stats..."
samtools flagstat "$BAM_SORT" > "${QCDIR}/${SAMPLE}.flagstat.txt"
samtools idxstats "$BAM_SORT" > "${QCDIR}/${SAMPLE}.idxstats.txt"

# ---------- Step 4: Filter MAPQ + remove secondary/supplementary ----------
log "Filtering MAPQ >= $MAPQ, keep primary alignments..."
BAM_FILT="${ALNDIR}/${SAMPLE}.mapq${MAPQ}.bam"
# -F 0x904 filters: 0x100 secondary, 0x800 supplementary, 0x4 unmapped
samtools view -@ "$THREADS" -b -q "$MAPQ" -F 0x904 "$BAM_SORT" > "$BAM_FILT"
samtools index "$BAM_FILT"

# ---------- Step 5: Mark/Remove duplicates (policy depends) ----------
# For CUT&Tag/TF narrow peaks: removing dups is often OK.
# For broad histone marks: sometimes keep dups (or mark only).
# We'll do:
#   auto: narrow -> yes ; broad -> no
if [[ "$REMOVE_DUP" == "auto" ]]; then
  if [[ "$PEAK_MODE" == "narrow" ]]; then REMOVE_DUP="yes"; else REMOVE_DUP="no"; fi
fi

BAM_FINAL="$BAM_FILT"
if [[ "$REMOVE_DUP" == "yes" ]]; then
  log "Removing duplicates (REMOVE_DUP=yes)..."
  BAM_DEDUP="${ALNDIR}/${SAMPLE}.mapq${MAPQ}.dedup.bam"

  if command -v picard >/dev/null 2>&1; then
    picard MarkDuplicates I="$BAM_FILT" O="$BAM_DEDUP" M="${QCDIR}/${SAMPLE}.dup_metrics.txt" REMOVE_DUPLICATES=true \
      |& tee "${LOGDIR}/${SAMPLE}.picard.dedup.log"
  else
    # samtools markdup requires name-sorted then fixmate
    TMP_NAME="${ALNDIR}/${SAMPLE}.name.bam"
    TMP_FIX="${ALNDIR}/${SAMPLE}.fixmate.bam"
    TMP_POS="${ALNDIR}/${SAMPLE}.pos.bam"
    samtools sort -n -@ "$THREADS" -o "$TMP_NAME" "$BAM_FILT"
    samtools fixmate -m -@ "$THREADS" "$TMP_NAME" "$TMP_FIX"
    samtools sort -@ "$THREADS" -o "$TMP_POS" "$TMP_FIX"
    samtools markdup -r -@ "$THREADS" "$TMP_POS" "$BAM_DEDUP" |& tee "${LOGDIR}/${SAMPLE}.samtools.markdup.log"
    rm -f "$TMP_NAME" "$TMP_FIX" "$TMP_POS"
    echo "samtools markdup used; metrics file not identical to picard." > "${QCDIR}/${SAMPLE}.dup_metrics.txt"
  fi

  samtools index "$BAM_DEDUP"
  BAM_FINAL="$BAM_DEDUP"
else
  log "Duplicate removal skipped (REMOVE_DUP=no)."
  # still useful to estimate duplicates if picard exists
  if command -v picard >/dev/null 2>&1; then
    picard MarkDuplicates I="$BAM_FILT" O=/dev/null M="${QCDIR}/${SAMPLE}.dup_metrics.txt" REMOVE_DUPLICATES=false \
      |& tee "${LOGDIR}/${SAMPLE}.picard.markdup.log" || true
  fi
fi

# Final stats
samtools flagstat "$BAM_FINAL" > "${QCDIR}/${SAMPLE}.final.flagstat.txt"

# ---------- Step 6: Peak Calling (MACS3 Adaptive Logic) ----------
# Identify enriched regions. 
# Switch between Standard ChIP (auto-model) and CUT&Tag (manual shift).
COUNT=0
log "Calling peaks with MACS3 (Protocol: $TECH)..."
mkdir -p "$PEAKDIR/$SAMPLE"

# 1. Determine Input Format
FORMAT="BAM"
[[ "$PE_MODE" == "yes" ]] && FORMAT="BAMPE"

# 2. Define Base Parameters
MACS3_ARGS=(-t "$BAM_FINAL" -f "$FORMAT" -g "$MACS_G" -n "$SAMPLE" --outdir "$PEAKDIR/$SAMPLE" -q 0.01)

# Add control/input if provided
[[ -n "${CONTROL_BAM}" ]] && MACS3_ARGS+=(-c "$CONTROL_BAM")

# 3. Execute Peak Calling
if [[ "$PEAK_MODE" == "narrow" ]]; then
    log "Peak Type: Narrow"
    
    if [[ "$TECH" == "cuttag" ]]; then
        # CUT&Tag Mode: High-resolution centering on Tn5 insertion sites
        log "Applying CUT&Tag-specific parameters: --nomodel --shift -100 --extsize 200"
        macs3 callpeak "${MACS3_ARGS[@]}" --nomodel --shift -100 --extsize 200 |& tee "${LOGDIR}/${SAMPLE}.macs3.log"
    else
        # Standard ChIP-seq Mode: 
        # For PE, BAMPE handles fragment size. For SE, MACS3 builds a model.
        log "Applying standard ChIP-seq logic (BAMPE/Modeling)."
        macs3 callpeak "${MACS3_ARGS[@]}" |& tee "${LOGDIR}/${SAMPLE}.macs3.log"
    fi
else
    # Broad Peak Mode: For repressive marks (H3K27me3, H3K9me3)
    log "Peak Type: Broad (using --broad and --broad-cutoff 0.1)"
    macs3 callpeak "${MACS3_ARGS[@]}" --broad --broad-cutoff 0.1 |& tee "${LOGDIR}/${SAMPLE}.macs3.log"
fi

# 4. Result Verification & Peak Counting
PEAK_FILE="${PEAKDIR}/${SAMPLE}/${SAMPLE}_peaks.narrowPeak"
[[ "$PEAK_MODE" == "broad" ]] && PEAK_FILE="${PEAKDIR}/${SAMPLE}/${SAMPLE}_peaks.broadPeak"

if [[ -f "$PEAK_FILE" ]]; then
    COUNT=$(wc -l < "$PEAK_FILE")
    log "Success: Found $COUNT peaks."
else
    log "ERROR: Peak file not found. Check ${LOGDIR}/${SAMPLE}.macs3.log"
fi

# ---------- Step 7: bigWig (deepTools bamCoverage) ----------
command -v bamCoverage >/dev/null 2>&1 || die "bamCoverage (deepTools) not found"

log "Generating bigWig (CPM-normalized) ..."
BW="${BWDIR}/${SAMPLE}.CPM.bw"

# Try to extract fragment size from MACS3 log if SE mode and standard ChIP
MACS_FRAGSIZE=""
if [[ "$SE_MODE" == "yes" && "$TECH" == "chip" && -f "${LOGDIR}/${SAMPLE}.macs3.log" ]]; then
    # MACS3 log usually contains: "# predicted fragment length is 123 bps"
    MACS_FRAGSIZE=$(grep "predicted fragment length is" "${LOGDIR}/${SAMPLE}.macs3.log" | tail -n 1 | sed -E 's/.*predicted fragment length is ([0-9]+) bps.*/\1/')
    if [[ -n "$MACS_FRAGSIZE" ]]; then
        log "Extracted fragment size from MACS3: $MACS_FRAGSIZE"
    fi
fi

# extendReads:
# - For PE BAM, deepTools can infer fragment length; for SE, extendReads is usually helpful.
EXT_ARG=()
if [[ "$EXTEND_READS" == "auto" ]]; then
  if [[ "$SE_MODE" == "yes" ]]; then EXTEND_READS="yes"; else EXTEND_READS="no"; fi
fi

if [[ "$EXTEND_READS" == "yes" ]]; then
  if [[ -n "$MACS_FRAGSIZE" ]]; then
      EXT_ARG=(--extendReads "$MACS_FRAGSIZE")
  else
      EXT_ARG=(--extendReads)
  fi
elif [[ "$EXTEND_READS" != "no" ]]; then
  # Assume user provided a specific integer length
  EXT_ARG=(--extendReads "$EXTEND_READS")
fi

bamCoverage \
  -b "$BAM_FINAL" \
  -o "$BW" \
  --normalizeUsing CPM \
  --binSize "$BINSIZE" \
  "${EXT_ARG[@]}" \
  --numberOfProcessors "$THREADS" \
  |& tee "${LOGDIR}/${SAMPLE}.bamCoverage.log"

# ---------- Optional: fingerprint plot (useful QC, but needs multiple samples usually) ----------
if command -v plotFingerprint >/dev/null 2>&1; then
  log "Optional QC: plotFingerprint (single sample still OK, more meaningful with many samples)"
  plotFingerprint -b "$BAM_FINAL" --labels "$SAMPLE" -plot "${QCDIR}/${SAMPLE}.fingerprint.png" \
    --numberOfProcessors "$THREADS" |& tee "${LOGDIR}/${SAMPLE}.fingerprint.log" || true
fi

log "DONE."
log "Final BAM : $BAM_FINAL"
log "bigWig     : $BW"
log "Peaks dir  : $PEAKDIR/$SAMPLE"

# ---------- Final Report ----------
echo "========================================================="
echo "        ChIP-seq/CUT&Tag Pipeline Summary"
echo "========================================================="
echo " Sample ID      : $SAMPLE"
echo " Protocol       : $TECH"
echo " Final BAM      : $BAM_FINAL"
echo " Peak File      : $PEAK_FILE"
echo " BigWig File    : $BW"
echo " Total Peaks    : $COUNT"
echo "========================================================="
log "Pipeline completed successfully."
