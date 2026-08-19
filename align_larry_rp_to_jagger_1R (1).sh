#!/bin/bash
#SBATCH --job-name=Jag1R_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00
#SBATCH --array=0
#SBATCH --output=Jag1R_align_%A_%a.out
#SBATCH --error=Jag1R_align_%A_%a.err

set -euo pipefail

# ============================================================================
# Align Larry (known 1RS.1BL carrier) and KS061278M-4 to the same competitive
# Jagger + Weining 1R reference.
#
# Align both samples:
#   sbatch --array=0-1%2 align_larry_rp_to_jagger_1R.sh
# ============================================================================

# ------------------------------ EDIT HERE -----------------------------------
ANALYSIS_ROOT="/bulk/akf/Addy/Variant_Calling/1RS_1BL_analysis"

REFERENCE="$ANALYSIS_ROOT/combined_reference/Jagger_plus_Weining_1R.fa"
OUTPUT_ROOT="$ANALYSIS_ROOT/alignments"
TEMP_ROOT="$ANALYSIS_ROOT/alignment_tmp"

LARRY_R1="$ANALYSIS_ROOT/reference_downloads/larry/fastq/SRR13572423_1.fastq.gz"
LARRY_R2="$ANALYSIS_ROOT/reference_downloads/larry/fastq/SRR13572423_2.fastq.gz"

RP_R1="/bulk/akf/Addy/fastq_QC/SRR15101982_1.fq.gz"
RP_R2="/bulk/akf/Addy/fastq_QC/SRR15101982_2.fq.gz"
# ---------------------------------------------------------------------------

BWA_THREADS=16
VIEW_THREADS=2
SORT_THREADS=6
SORT_MEMORY_PER_THREAD="4G"
QC_THREADS=6

JAGGER_1B_ID="LR862509.1"
JAGGER_1B_LENGTH="705338699"
RYE_1R_ID="Weining_1R"
RYE_1R_LENGTH="940967580"

SAMPLES=("Larry" "KS061278M-4")
R1_FILES=("$LARRY_R1" "$RP_R1")
R2_FILES=("$LARRY_R2" "$RP_R2")

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if ! [[ "$TASK_ID" =~ ^[0-9]+$ ]] || (( TASK_ID >= ${#SAMPLES[@]} )); then
    echo "ERROR: array index must be 0 (Larry) or 1 (KS061278M-4)." >&2
    exit 1
fi

SAMPLE="${SAMPLES[$TASK_ID]}"
R1="${R1_FILES[$TASK_ID]}"
R2="${R2_FILES[$TASK_ID]}"

SAMPLE_DIR="$OUTPUT_ROOT/$SAMPLE"
FINAL_BAM="$SAMPLE_DIR/${SAMPLE}.Jagger_plus_Weining_1R.mapped.sorted.bam"
PARTIAL_BAM="$SAMPLE_DIR/.${SAMPLE}.Jagger_plus_Weining_1R.mapped.sorted.partial.bam"
ALIGNMENT_LOG="$SAMPLE_DIR/${SAMPLE}.bwa_mem.log"
RUN_INFO="$SAMPLE_DIR/${SAMPLE}.run_info.txt"
TMP_DIR="$TEMP_ROOT/${SAMPLE}.${SLURM_JOB_ID:-manual}.${TASK_ID}"

mkdir -p "$SAMPLE_DIR" "$TEMP_ROOT"

module purge
module load SAMtools
module load BWA/0.7.17-GCCcore-11.3.0

for program in bwa samtools awk grep sort; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: required program is unavailable: $program" >&2
        exit 1
    fi
done

if [[ "$SAMPLE" == "KS061278M-4" ]] &&
   { [[ "$R1" == REPLACE_WITH_* ]] || [[ "$R2" == REPLACE_WITH_* ]]; }; then
    echo "ERROR: enter the two original KS061278M-4 FASTQ paths at the top" >&2
    echo "       of this script before submitting array task 1." >&2
    exit 1
fi

for input_file in "$REFERENCE" "${REFERENCE}.fai" "$R1" "$R2"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: required input is missing or empty:" >&2
        echo "  $input_file" >&2
        exit 1
    fi
done

BWA_INDEX_FILES=(
    "${REFERENCE}.amb"
    "${REFERENCE}.ann"
    "${REFERENCE}.bwt"
    "${REFERENCE}.pac"
    "${REFERENCE}.sa"
)
for index_file in "${BWA_INDEX_FILES[@]}"; do
    if [[ ! -s "$index_file" ]]; then
        echo "ERROR: BWA index file is missing or empty:" >&2
        echo "  $index_file" >&2
        exit 1
    fi
done

observed_1b_length="$(
    awk -F '\t' -v id="$JAGGER_1B_ID" '$1 == id {print $2}' "${REFERENCE}.fai"
)"
observed_1r_length="$(
    awk -F '\t' -v id="$RYE_1R_ID" '$1 == id {print $2}' "${REFERENCE}.fai"
)"

if [[ "$observed_1b_length" != "$JAGGER_1B_LENGTH" ]]; then
    echo "ERROR: $JAGGER_1B_ID is absent or has an unexpected length." >&2
    exit 1
fi
if [[ "$observed_1r_length" != "$RYE_1R_LENGTH" ]]; then
    echo "ERROR: $RYE_1R_ID is absent or has an unexpected length." >&2
    exit 1
fi

cleanup_temp() {
    if [[ -n "${TMP_DIR:-}" && "$TMP_DIR" == "$TEMP_ROOT/"* &&
          "$TMP_DIR" != "$TEMP_ROOT" ]]; then
        rm -rf -- "$TMP_DIR"
    fi
}
trap cleanup_temp EXIT

publish_partial_bam_if_complete() {
    if [[ ! -s "$PARTIAL_BAM" ]]; then
        return 1
    fi

    if samtools quickcheck -v "$PARTIAL_BAM"; then
        if ! samtools view -H "$PARTIAL_BAM" |
             awk '$1 == "@HD" && $0 ~ /SO:coordinate/ {found=1}
                  END {exit !found}'; then
            echo "ERROR: completed partial BAM is not coordinate sorted:" >&2
            echo "  $PARTIAL_BAM" >&2
            exit 1
        fi
        echo "A complete partial BAM was found; publishing it without realignment."
        mv "$PARTIAL_BAM" "$FINAL_BAM"
        return 0
    fi

    echo "Removing an incomplete partial BAM from an interrupted run."
    rm -f -- "$PARTIAL_BAM"
    return 1
}

echo "Alignment job started: $(date)"
echo "Sample: $SAMPLE"
echo "Reference: $REFERENCE"
echo "R1: $R1"
echo "R2: $R2"
echo "Output BAM: $FINAL_BAM"

NEED_ALIGNMENT=true
if [[ -e "$FINAL_BAM" ]]; then
    if [[ ! -s "$FINAL_BAM" ]] || ! samtools quickcheck -v "$FINAL_BAM"; then
        echo "ERROR: final BAM exists but is empty or invalid; it was not overwritten:" >&2
        echo "  $FINAL_BAM" >&2
        exit 1
    fi
    echo "A complete final BAM already exists; alignment will be skipped."
    NEED_ALIGNMENT=false
elif publish_partial_bam_if_complete; then
    NEED_ALIGNMENT=false
fi

if [[ "$NEED_ALIGNMENT" == true ]]; then
    mkdir -p "$TMP_DIR"

    echo "Running BWA-MEM and coordinate sorting: $(date)"
    bwa mem \
        -t "$BWA_THREADS" \
        -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA" \
        "$REFERENCE" \
        "$R1" \
        "$R2" \
        2> "$ALIGNMENT_LOG" |
        samtools view \
            -@ "$VIEW_THREADS" \
            -u \
            -F 4 \
            - |
        samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEMORY_PER_THREAD" \
            -T "$TMP_DIR/sort" \
            -o "$PARTIAL_BAM" \
            -

    samtools quickcheck -v "$PARTIAL_BAM"
    if ! samtools view -H "$PARTIAL_BAM" |
         awk '$1 == "@HD" && $0 ~ /SO:coordinate/ {found=1}
              END {exit !found}'; then
        echo "ERROR: output BAM is not marked as coordinate sorted." >&2
        exit 1
    fi

    mv "$PARTIAL_BAM" "$FINAL_BAM"
fi

if [[ ! -s "${FINAL_BAM}.csi" ]] ||
   ! samtools view -c "$FINAL_BAM" "$RYE_1R_ID" >/dev/null 2>&1; then
    rm -f -- "${FINAL_BAM}.csi"
    echo "Building CSI index: $(date)"
    samtools index -@ "$QC_THREADS" -c "$FINAL_BAM"
fi

samtools quickcheck -v "$FINAL_BAM"
samtools idxstats "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE}.idxstats.tsv"
samtools flagstat -@ "$QC_THREADS" "$FINAL_BAM" \
    > "$SAMPLE_DIR/${SAMPLE}.flagstat.txt"
samtools coverage -r "$JAGGER_1B_ID" "$FINAL_BAM" \
    > "$SAMPLE_DIR/${SAMPLE}.${JAGGER_1B_ID}.coverage.tsv"
samtools coverage -r "$RYE_1R_ID" "$FINAL_BAM" \
    > "$SAMPLE_DIR/${SAMPLE}.${RYE_1R_ID}.coverage.tsv"

# Record software versions, exact inputs, and primary filtering decisions.
{
    echo "sample=$SAMPLE"
    echo "completed=$(date --iso-8601=seconds)"
    echo "reference=$REFERENCE"
    echo "read1=$R1"
    echo "read2=$R2"
    echo "bam=$FINAL_BAM"
    echo "jagger_1b_id=$JAGGER_1B_ID"
    echo "weining_1r_id=$RYE_1R_ID"
    echo "bwa_threads=$BWA_THREADS"
    echo "samtools_view_filter=-F 4"
    echo "mapq_filter=none"
    echo "duplicate_filter=none"
    echo "proper_pair_filter=none"
    echo "bwa_version=$({ bwa 2>&1 || true; } | awk '/Version:/ {version=$2} END {print version}')"
    echo "samtools_version=$(samtools --version | awk 'NR == 1 {print $2}')"
} > "$RUN_INFO"

echo "Alignment and QC completed successfully: $(date)"
echo "Final BAM: $FINAL_BAM"
echo "CSI index: ${FINAL_BAM}.csi"
echo "QC directory: $SAMPLE_DIR"
