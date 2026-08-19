#!/bin/bash
#SBATCH --job-name=Jag1R_controls
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-7
#SBATCH --output=Jag1R_controls_%A_%a.out
#SBATCH --error=Jag1R_controls_%A_%a.err

set -euo pipefail

# Align KS061406LN~26 and seven wild-emmer donors to the same complete Jagger + Weining chromosome 1R competitive reference used for Larry and KS061278M-4
#
# Submit all samples:
#   sbatch align_rp2_emmer_to_jagger_1R.sh

# ------------------------------ EDIT HERE -----------------------------------
ANALYSIS_ROOT="/bulk/akf/Addy/Variant_Calling/1RS_1BL_analysis"
FASTQ_ROOT="/bulk/akf/Addy/fastq_QC"
REFERENCE="$ANALYSIS_ROOT/combined_reference/Jagger_plus_Weining_1R.fa"
# ---------------------------------------------------------------------------

OUTPUT_ROOT="$ANALYSIS_ROOT/alignments"
TEMP_ROOT="$ANALYSIS_ROOT/alignment_tmp"

BWA_THREADS=16
VIEW_THREADS=2
SORT_THREADS=6
SORT_MEMORY_PER_THREAD="4G"
QC_THREADS=6

JAGGER_1B_ID="LR862509.1"
JAGGER_1B_LENGTH="705338699"
RYE_1R_ID="Weining_1R"
RYE_1R_LENGTH="940967580"

SAMPLES=(
    "KS061406LN-26"
    "PI487264"
    "PI428088"
    "PI471021"
    "PI466984"
    "PI471815"
    "PI470736"
    "PI414718"
)

DATA_TYPES=(
    "WGS"
    "EXOME"
    "EXOME"
    "EXOME"
    "WGS"
    "WGS"
    "WGS"
    "WGS"
)

RUN_PREFIX_LISTS=(
    "SRR15129972 SRR15129973 SRR15129974 SRR15129975"
    "SRR8157374"
    "SRR8157382"
    "SRR8157418"
    "WE-1"
    "WE-2"
    "WE-3"
    "WE-4"
)

declare -A ALLOW_LATE_PAIR_MISMATCH=(
    [SRR8157374]=1
    [SRR8157418]=1
)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if ! [[ "$TASK_ID" =~ ^[0-9]+$ ]] || (( TASK_ID >= ${#SAMPLES[@]} )); then
    echo "ERROR: array index must be between 0 and 7." >&2
    exit 1
fi

SAMPLE="${SAMPLES[$TASK_ID]}"
DATA_TYPE="${DATA_TYPES[$TASK_ID]}"
read -r -a RUN_PREFIXES <<< "${RUN_PREFIX_LISTS[$TASK_ID]}"

SAMPLE_DIR="$OUTPUT_ROOT/$SAMPLE"
RUN_BAM_DIR="$SAMPLE_DIR/run_bams"
FINAL_BAM="$SAMPLE_DIR/${SAMPLE}.Jagger_plus_Weining_1R.mapped.sorted.bam"
PARTIAL_FINAL_BAM="$SAMPLE_DIR/.${SAMPLE}.Jagger_plus_Weining_1R.mapped.sorted.partial.bam"
RUN_INFO="$SAMPLE_DIR/${SAMPLE}.run_info.txt"
TMP_DIR="$TEMP_ROOT/${SAMPLE}.${SLURM_JOB_ID:-manual}.${TASK_ID}"

mkdir -p "$SAMPLE_DIR" "$RUN_BAM_DIR" "$TMP_DIR"

module purge
module load SAMtools
module load BWA/0.7.17-GCCcore-11.3.0

for program in bwa samtools awk; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: required program is unavailable: $program" >&2
        exit 1
    fi
done

for input_file in "$REFERENCE" "${REFERENCE}.fai"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: required reference file is missing or empty: $input_file" >&2
        exit 1
    fi
done

for suffix in amb ann bwt pac sa; do
    if [[ ! -s "${REFERENCE}.${suffix}" ]]; then
        echo "ERROR: BWA index file is missing: ${REFERENCE}.${suffix}" >&2
        exit 1
    fi
done

observed_1b_length="$(awk -F '\t' -v id="$JAGGER_1B_ID" '$1 == id {print $2}' "${REFERENCE}.fai")"
observed_1r_length="$(awk -F '\t' -v id="$RYE_1R_ID" '$1 == id {print $2}' "${REFERENCE}.fai")"
if [[ "$observed_1b_length" != "$JAGGER_1B_LENGTH" ||
      "$observed_1r_length" != "$RYE_1R_LENGTH" ]]; then
    echo "ERROR: the reference does not contain the expected Jagger 1B and Weining 1R records." >&2
    exit 1
fi

cleanup_temp() {
    if [[ -n "${TMP_DIR:-}" && "$TMP_DIR" == "$TEMP_ROOT/"* &&
          "$TMP_DIR" != "$TEMP_ROOT" ]]; then
        rm -rf -- "$TMP_DIR"
    fi
}
trap cleanup_temp EXIT

align_one_run() {
    local prefix="$1"
    local r1="$FASTQ_ROOT/${prefix}_1.fq.gz"
    local r2="$FASTQ_ROOT/${prefix}_2.fq.gz"
    local run_bam="$RUN_BAM_DIR/${prefix}.mapped.sorted.bam"
    local partial_bam="$RUN_BAM_DIR/.${prefix}.mapped.sorted.partial.bam"
    local log="$RUN_BAM_DIR/${prefix}.bwa_mem.log"
    local run_tmp="$TMP_DIR/${prefix}"
    local -a pipeline_status
    local bwa_status view_status sort_status

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "ERROR: paired FASTQs are missing or empty for $prefix:" >&2
        echo "  $r1" >&2
        echo "  $r2" >&2
        exit 1
    fi

    if [[ -s "$run_bam" ]] && samtools quickcheck -v "$run_bam"; then
        echo "Reusing completed run BAM: $run_bam"
        return
    fi
    if [[ -e "$run_bam" ]]; then
        echo "ERROR: existing run BAM is invalid and was not overwritten: $run_bam" >&2
        exit 1
    fi

    if [[ -s "$partial_bam" ]] && samtools quickcheck -v "$partial_bam"; then
        if [[ -s "$log" ]] &&
           grep -Eq 'paired reads have different names|paired reads have different name|mem_sam_pe' "$log"; then
            if [[ "${ALLOW_LATE_PAIR_MISMATCH[$prefix]:-0}" != 1 ]]; then
                echo "ERROR: an existing partial BAM ended at an unapproved mate-name mismatch: $prefix" >&2
                exit 1
            fi
            echo "WARNING: Publishing the valid BAM retained at the verified late mismatch for $prefix." >&2
        else
            echo "Publishing completed partial BAM for $prefix."
        fi
        mv "$partial_bam" "$run_bam"
        return
    fi

    rm -f -- "$partial_bam"
    mkdir -p "$run_tmp"
    echo "Aligning $SAMPLE run $prefix: $(date)"

    set +e
    bwa mem \
        -t "$BWA_THREADS" \
        -R "@RG\\tID:${prefix}\\tSM:${SAMPLE}\\tPL:ILLUMINA" \
        "$REFERENCE" "$r1" "$r2" \
        2> "$log" |
        samtools view -@ "$VIEW_THREADS" -u -F 4 - |
        samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEMORY_PER_THREAD" \
            -T "$run_tmp/sort" \
            -o "$partial_bam" -
    pipeline_status=("${PIPESTATUS[@]}")
    set -e

    bwa_status="${pipeline_status[0]}"
    view_status="${pipeline_status[1]}"
    sort_status="${pipeline_status[2]}"

    if (( view_status != 0 || sort_status != 0 )); then
        echo "ERROR: SAMtools failed for $prefix (view=$view_status, sort=$sort_status)." >&2
        rm -f -- "$partial_bam"
        exit 1
    fi

    if (( bwa_status != 0 )); then
        if [[ "${ALLOW_LATE_PAIR_MISMATCH[$prefix]:-0}" == 1 ]] &&
           grep -Eq 'paired reads have different names|paired reads have different name|mem_sam_pe' "$log"; then
            echo "WARNING: BWA stopped at the verified late mate-name mismatch for $prefix." >&2
            echo "WARNING: Retaining the valid mapped alignment produced before the mismatch." >&2
        else
            echo "ERROR: BWA failed for $prefix with exit status $bwa_status." >&2
            echo "See: $log" >&2
            rm -f -- "$partial_bam"
            exit 1
        fi
    fi

    samtools quickcheck -v "$partial_bam"
    mv "$partial_bam" "$run_bam"
}

echo "Alignment job started: $(date)"
echo "Sample: $SAMPLE"
echo "Data type: $DATA_TYPE"
echo "Runs: ${RUN_PREFIXES[*]}"
echo "Reference: $REFERENCE"

if [[ -s "$FINAL_BAM" ]] && samtools quickcheck -v "$FINAL_BAM"; then
    echo "Reusing completed biological-sample BAM: $FINAL_BAM"
else
    if [[ -e "$FINAL_BAM" ]]; then
        echo "ERROR: existing final BAM is invalid and was not overwritten: $FINAL_BAM" >&2
        exit 1
    fi

    if [[ -s "$PARTIAL_FINAL_BAM" ]] &&
       samtools quickcheck -v "$PARTIAL_FINAL_BAM"; then
        echo "Publishing completed partial biological-sample BAM."
        mv "$PARTIAL_FINAL_BAM" "$FINAL_BAM"
    else
        rm -f -- "$PARTIAL_FINAL_BAM"

        RUN_BAMS=()
        for prefix in "${RUN_PREFIXES[@]}"; do
            align_one_run "$prefix"
            RUN_BAMS+=("$RUN_BAM_DIR/${prefix}.mapped.sorted.bam")
        done

        if (( ${#RUN_BAMS[@]} == 1 )); then
            if ! ln "${RUN_BAMS[0]}" "$PARTIAL_FINAL_BAM" 2>/dev/null; then
                cp --reflink=auto "${RUN_BAMS[0]}" "$PARTIAL_FINAL_BAM"
            fi
        else
            echo "Merging ${#RUN_BAMS[@]} runs for $SAMPLE: $(date)"
            samtools merge -@ "$QC_THREADS" -o "$PARTIAL_FINAL_BAM" "${RUN_BAMS[@]}"
        fi

        samtools quickcheck -v "$PARTIAL_FINAL_BAM"
        mv "$PARTIAL_FINAL_BAM" "$FINAL_BAM"
    fi
fi

if [[ ! -s "${FINAL_BAM}.csi" ]] ||
   ! samtools view -c "$FINAL_BAM" "$RYE_1R_ID" >/dev/null 2>&1; then
    rm -f -- "${FINAL_BAM}.csi"
    echo "Building CSI index: $(date)"
    samtools index -@ "$QC_THREADS" -c "$FINAL_BAM"
fi

samtools quickcheck -v "$FINAL_BAM"
samtools idxstats "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE}.idxstats.tsv"
samtools flagstat -@ "$QC_THREADS" "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE}.flagstat.txt"
samtools coverage -r "$JAGGER_1B_ID" "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE}.${JAGGER_1B_ID}.coverage.tsv"
samtools coverage -r "$RYE_1R_ID" "$FINAL_BAM" > "$SAMPLE_DIR/${SAMPLE}.${RYE_1R_ID}.coverage.tsv"

{
    echo "sample=$SAMPLE"
    echo "display_name=$([[ "$SAMPLE" == "KS061406LN-26" ]] && echo 'KS061406LN~26' || echo "$SAMPLE")"
    echo "data_type=$DATA_TYPE"
    echo "completed=$(date --iso-8601=seconds)"
    echo "reference=$REFERENCE"
    echo "fastq_root=$FASTQ_ROOT"
    echo "runs=${RUN_PREFIXES[*]}"
    echo "bam=$FINAL_BAM"
    echo "samtools_view_filter=-F 4"
    echo "mapq_filter=none"
    echo "duplicate_filter=none"
    echo "proper_pair_filter=none"
    for prefix in "${RUN_PREFIXES[@]}"; do
        echo "allow_late_pair_mismatch_${prefix}=${ALLOW_LATE_PAIR_MISMATCH[$prefix]:-0}"
    done
    echo "bwa_version=$({ bwa 2>&1 || true; } | awk '/Version:/ {version=$2} END {print version}')"
    echo "samtools_version=$(samtools --version | awk 'NR == 1 {print $2}')"
} > "$RUN_INFO"

echo "Alignment and QC completed successfully: $(date)"
echo "Final BAM: $FINAL_BAM"
echo "CSI index: ${FINAL_BAM}.csi"
