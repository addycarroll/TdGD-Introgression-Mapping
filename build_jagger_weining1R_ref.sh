#!/bin/bash
#SBATCH --job-name=build_Jag_1R_ref
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=build_Jag_1R_ref_%j.out
#SBATCH --error=build_Jag_1R_ref_%j.err

set -euo pipefail

# ============================================================
# Build and index the competitive reference for 1RS.1BL testing
#
# Reference contents:
#   - complete Jagger wheat assembly, GCA_903993795.1
#   - Weining rye chromosome 1R only, GCA_016097815.1
#
# The rye record is renamed "Weining_1R" so it is easy to identify in BAM files, coverage tables, and plots.
# ============================================================

# -------------------------- EDIT IF NEEDED -------------------
ANALYSIS_ROOT="/bulk/akf/Addy/Variant_Calling/1RS_1BL_analysis"

JAGGER_FASTA_GZ="$ANALYSIS_ROOT/reference_downloads/references/GCA_903993795.1_10wheat_assembly_jagger_genomic.fna.gz"
JAGGER_REPORT="$ANALYSIS_ROOT/reference_downloads/references/GCA_903993795.1_10wheat_assembly_jagger_assembly_report.txt"

RYE_1R_FASTA_GZ="$ANALYSIS_ROOT/reference_downloads/references/Weining_chr1R.fna.gz"
RYE_REPORT="$ANALYSIS_ROOT/reference_downloads/references/GCA_016097815.1_HAU_Weining_v1.0_assembly_report.txt"

OUTPUT_DIR="$ANALYSIS_ROOT/combined_reference"
COMBINED_FASTA="$OUTPUT_DIR/Jagger_plus_Weining_1R.fa"
# --------------------------------------------------------------

THREADS="${SLURM_CPUS_PER_TASK:-4}"
RYE_OUTPUT_ID="Weining_1R"

mkdir -p "$OUTPUT_DIR"

module purge
module load SAMtools
module load BWA/0.7.17-GCCcore-11.3.0

for program in gzip awk grep sort uniq samtools bwa; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: required program not found after loading modules: $program" >&2
        exit 1
    fi
done

for input_file in \
    "$JAGGER_FASTA_GZ" \
    "$JAGGER_REPORT" \
    "$RYE_1R_FASTA_GZ" \
    "$RYE_REPORT"
do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: required input is missing or empty:" >&2
        echo "  $input_file" >&2
        exit 1
    fi
done

echo "Build started: $(date)"
echo "Combined reference: $COMBINED_FASTA"
echo "Threads: $THREADS"

echo "Testing compressed input FASTAs."
gzip -t "$JAGGER_FASTA_GZ"
gzip -t "$RYE_1R_FASTA_GZ"

# Obtain chromosome accessions and lengths from the reports.
JAGGER_1B_ACCESSION="$(
    awk -F '\t' '
        $0 !~ /^#/ && $2 == "assembled-molecule" && $3 == "1B" {
            print $5
        }
    ' "$JAGGER_REPORT"
)"

JAGGER_1B_LENGTH="$(
    awk -F '\t' '
        $0 !~ /^#/ && $2 == "assembled-molecule" && $3 == "1B" {
            print $9
        }
    ' "$JAGGER_REPORT"
)"

RYE_1R_ACCESSION="$(
    awk -F '\t' '
        $0 !~ /^#/ && $2 == "assembled-molecule" && $3 == "1R" {
            print $5
        }
    ' "$RYE_REPORT"
)"

RYE_1R_LENGTH="$(
    awk -F '\t' '
        $0 !~ /^#/ && $2 == "assembled-molecule" && $3 == "1R" {
            print $9
        }
    ' "$RYE_REPORT"
)"

JAGGER_REPORT_RECORDS="$(
    awk -F '\t' '$0 !~ /^#/ && NF >= 9 {count++} END {print count + 0}' \
        "$JAGGER_REPORT"
)"

if [[ "$JAGGER_1B_ACCESSION" != "LR862509.1" ||
      "$JAGGER_1B_LENGTH" != "705338699" ]]; then
    echo "ERROR: unexpected Jagger chromosome 1B identity or length." >&2
    echo "  observed: $JAGGER_1B_ACCESSION  $JAGGER_1B_LENGTH" >&2
    exit 1
fi

if [[ "$RYE_1R_ACCESSION" != "CM027771.1" ||
      "$RYE_1R_LENGTH" != "940967580" ]]; then
    echo "ERROR: unexpected Weining chromosome 1R identity or length." >&2
    echo "  observed: $RYE_1R_ACCESSION  $RYE_1R_LENGTH" >&2
    exit 1
fi

read -r RYE_FASTA_RECORDS RYE_FASTA_LENGTH < <(
    gzip -cd "$RYE_1R_FASTA_GZ" |
        awk '
            /^>/ {
                records++
                next
            }
            {
                length_sum += length($0)
            }
            END {
                printf "%.0f %.0f\n", records + 0, length_sum + 0
            }
        '
)

if [[ "$RYE_FASTA_RECORDS" -ne 1 ||
      "$RYE_FASTA_LENGTH" -ne "$RYE_1R_LENGTH" ]]; then
    echo "ERROR: extracted rye 1R FASTA failed validation." >&2
    echo "  records: $RYE_FASTA_RECORDS (expected 1)" >&2
    echo "  length:  $RYE_FASTA_LENGTH (expected $RYE_1R_LENGTH)" >&2
    exit 1
fi

EXPECTED_COMBINED_RECORDS=$((JAGGER_REPORT_RECORDS + 1))

validate_combined_reference() {
    local fasta="$1"
    local fai="${fasta}.fai"
    local observed_records
    local observed_jagger_1b_length
    local observed_rye_1r_length
    local duplicate_ids

    if [[ ! -s "$fasta" ]]; then
        return 1
    fi

    rm -f "$fai"
    samtools faidx "$fasta"

    observed_records="$(awk 'END {print NR + 0}' "$fai")"
    observed_jagger_1b_length="$(
        awk -F '\t' -v id="$JAGGER_1B_ACCESSION" '$1 == id {print $2}' "$fai"
    )"
    observed_rye_1r_length="$(
        awk -F '\t' -v id="$RYE_OUTPUT_ID" '$1 == id {print $2}' "$fai"
    )"
    duplicate_ids="$(
        cut -f 1 "$fai" | sort | uniq -d
    )"

    if [[ "$observed_records" -ne "$EXPECTED_COMBINED_RECORDS" ]]; then
        echo "ERROR: combined FASTA has $observed_records records;" >&2
        echo "       expected $EXPECTED_COMBINED_RECORDS." >&2
        return 1
    fi

    if [[ "$observed_jagger_1b_length" != "$JAGGER_1B_LENGTH" ]]; then
        echo "ERROR: Jagger 1B is absent or has the wrong length." >&2
        return 1
    fi

    if [[ "$observed_rye_1r_length" != "$RYE_1R_LENGTH" ]]; then
        echo "ERROR: renamed Weining 1R is absent or has the wrong length." >&2
        return 1
    fi

    if [[ -n "$duplicate_ids" ]]; then
        echo "ERROR: duplicate sequence IDs found in combined FASTA:" >&2
        printf '%s\n' "$duplicate_ids" >&2
        return 1
    fi
}

# ============================================================
# 1. Create the combined reference atomically
# ============================================================

if [[ -s "$COMBINED_FASTA" ]] &&
   validate_combined_reference "$COMBINED_FASTA"; then
    echo "Validated combined FASTA already exists; reusing it."
else
    TEMP_FASTA="$OUTPUT_DIR/.Jagger_plus_Weining_1R.${SLURM_JOB_ID:-manual}.tmp.fa"
    TEMP_FAI="${TEMP_FASTA}.fai"

    cleanup_temp() {
        rm -f "$TEMP_FASTA" "$TEMP_FAI"
    }
    trap cleanup_temp EXIT

    echo "Creating combined FASTA."
    echo "  Jagger: $JAGGER_FASTA_GZ"
    echo "  rye 1R: $RYE_1R_FASTA_GZ"

    gzip -cd "$JAGGER_FASTA_GZ" > "$TEMP_FASTA"

    gzip -cd "$RYE_1R_FASTA_GZ" |
        awk -v new_id="$RYE_OUTPUT_ID" '
            /^>/ {
                print ">" new_id
                next
            }
            {
                print
            }
        ' >> "$TEMP_FASTA"

    validate_combined_reference "$TEMP_FASTA"

    mv "$TEMP_FASTA" "$COMBINED_FASTA"
    mv "$TEMP_FAI" "${COMBINED_FASTA}.fai"
    trap - EXIT

    echo "Combined FASTA created and validated."
fi

if [[ ! -s "${COMBINED_FASTA}.fai" ]]; then
    samtools faidx "$COMBINED_FASTA"
fi

# ============================================================
# 2. Build the BWA index
# ============================================================

BWA_INDEX_FILES=(
    "${COMBINED_FASTA}.amb"
    "${COMBINED_FASTA}.ann"
    "${COMBINED_FASTA}.bwt"
    "${COMBINED_FASTA}.pac"
    "${COMBINED_FASTA}.sa"
)

BWA_INDEX_COMPLETE=true
for index_file in "${BWA_INDEX_FILES[@]}"; do
    if [[ ! -s "$index_file" ]]; then
        BWA_INDEX_COMPLETE=false
        break
    fi
done

if [[ "$BWA_INDEX_COMPLETE" == true ]]; then
    echo "Complete BWA index already exists; skipping indexing."
else
    echo "Removing any incomplete BWA index files."
    rm -f \
        "${COMBINED_FASTA}.amb" \
        "${COMBINED_FASTA}.ann" \
        "${COMBINED_FASTA}.bwt" \
        "${COMBINED_FASTA}.pac" \
        "${COMBINED_FASTA}.sa"

    echo "Building BWA index with the bwtsw algorithm."
    echo "The bwtsw algorithm is required for references larger than 2 GB."
    bwa index -a bwtsw "$COMBINED_FASTA"

    for index_file in "${BWA_INDEX_FILES[@]}"; do
        if [[ ! -s "$index_file" ]]; then
            echo "ERROR: expected BWA index file was not created:" >&2
            echo "  $index_file" >&2
            exit 1
        fi
    done
fi

# ============================================================
# 3. Final summary
# ============================================================

echo
echo "Reference build completed successfully: $(date)"
echo
echo "Combined reference:"
echo "  $COMBINED_FASTA"
echo "samtools index:"
echo "  ${COMBINED_FASTA}.fai"
echo "Jagger chromosome 1B ID:"
echo "  $JAGGER_1B_ACCESSION"
echo "Weining chromosome 1R ID:"
echo "  $RYE_OUTPUT_ID"
echo "Total reference records:"
awk 'END {print "  " NR}' "${COMBINED_FASTA}.fai"
echo
echo "Use this FASTA for every Larry and KS061278M-4 alignment:"
echo "  $COMBINED_FASTA"
