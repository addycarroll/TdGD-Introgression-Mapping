#!/bin/bash
#SBATCH --job-name=get_refs_larry
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=get_refs_larry_%j.out
#SBATCH --error=get_refs_larry_%j.err

set -euo pipefail

# ============================================================
# Download inputs for a Parker et al.-style 1RS.1BL analysis
#
# Downloads:
#   1. Complete Jagger wheat assembly (GCA_903993795.1)
#   2. Complete Weining rye assembly (GCA_016097815.1)
#   3. Larry paired-end reads (SRR13572423)
#
# Also extracts:
#   - Weining rye chromosome 1R
#   - Jagger wheat chromosome 1B (for inspection only)
#
# IMPORTANT:
# Use the COMPLETE Jagger assembly, plus rye chromosome 1R, for
# competitive alignment. Do not map against only Jagger 1B + rye 1R.
# ============================================================

# -------------------------- EDIT THIS -------------------------
DOWNLOAD_ROOT="/bulk/akf/Addy/Variant_Calling/1RS_1BL_analysis/reference_downloads"
# --------------------------------------------------------------

THREADS="${SLURM_CPUS_PER_TASK:-8}"

JAGGER_ACCESSION="GCA_903993795.1"
RYE_ACCESSION="GCA_016097815.1"
LARRY_RUN="SRR13572423"

JAGGER_STEM="GCA_903993795.1_10wheat_assembly_jagger"
RYE_STEM="GCA_016097815.1_HAU_Weining_v1.0"

JAGGER_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/993/795/${JAGGER_STEM}/${JAGGER_STEM}_genomic.fna.gz"
JAGGER_REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/993/795/${JAGGER_STEM}/${JAGGER_STEM}_assembly_report.txt"

RYE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/${RYE_STEM}/${RYE_STEM}_genomic.fna.gz"
RYE_REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/097/815/${RYE_STEM}/${RYE_STEM}_assembly_report.txt"

REFERENCE_DIR="$DOWNLOAD_ROOT/references"
LARRY_DIR="$DOWNLOAD_ROOT/larry"
SRA_CACHE="$LARRY_DIR/sra"
FASTQ_DIR="$LARRY_DIR/fastq"
TEMP_DIR="$DOWNLOAD_ROOT/fasterq_tmp_${SLURM_JOB_ID:-manual}"

JAGGER_FASTA_GZ="$REFERENCE_DIR/${JAGGER_STEM}_genomic.fna.gz"
JAGGER_REPORT="$REFERENCE_DIR/${JAGGER_STEM}_assembly_report.txt"
JAGGER_1B_GZ="$REFERENCE_DIR/Jagger_chr1B.fna.gz"

RYE_FASTA_GZ="$REFERENCE_DIR/${RYE_STEM}_genomic.fna.gz"
RYE_REPORT="$REFERENCE_DIR/${RYE_STEM}_assembly_report.txt"
RYE_1R_GZ="$REFERENCE_DIR/Weining_chr1R.fna.gz"

mkdir -p "$REFERENCE_DIR" "$SRA_CACHE" "$FASTQ_DIR" "$TEMP_DIR"

echo "Download started: $(date)"
echo "Destination: $DOWNLOAD_ROOT"
echo "Threads: $THREADS"

# Load the cluster's SRA Toolkit and compression program.
module purge
module load SRA-Toolkit
module load pigz

for program in wget gzip awk grep wc prefetch fasterq-dump vdb-validate pigz; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    if [[ -s "$destination" ]]; then
        echo "Already present: $destination"
        return
    fi

    echo "Downloading: $url"
    wget --continue --output-document="$partial" "$url"
    mv "$partial" "$destination"
}

extract_chromosome() {
    local source_fasta_gz="$1"
    local assembly_report="$2"
    local chromosome_label="$3"
    local destination_fasta_gz="$4"
    local chromosome_accession

    chromosome_accession="$(
        awk -F '\t' -v chromosome="$chromosome_label" '
            $0 !~ /^#/ &&
            $2 == "assembled-molecule" &&
            $3 == chromosome {
                print $5
            }
        ' "$assembly_report"
    )"

    if [[ -z "$chromosome_accession" ]]; then
        echo "ERROR: chromosome $chromosome_label was not found in:" >&2
        echo "  $assembly_report" >&2
        exit 1
    fi

    if [[ "$(printf '%s\n' "$chromosome_accession" | wc -l)" -ne 1 ]]; then
        echo "ERROR: chromosome $chromosome_label has multiple accessions in:" >&2
        echo "  $assembly_report" >&2
        printf '%s\n' "$chromosome_accession" >&2
        exit 1
    fi

    echo "Extracting chromosome $chromosome_label ($chromosome_accession) from:"
    echo "  $source_fasta_gz"

    gzip -cd "$source_fasta_gz" |
        awk -v accession="$chromosome_accession" -v chromosome="$chromosome_label" '
            BEGIN {
                keep = 0
                count = 0
            }
            /^>/ {
                header_accession = $1
                sub(/^>/, "", header_accession)
                keep = (header_accession == accession)
                if (keep) {
                    count++
                }
            }
            keep {
                print
            }
            END {
                if (count != 1) {
                    print "ERROR: expected one " chromosome \
                          " FASTA record but found " count > "/dev/stderr"
                    exit 1
                }
            }
        ' |
        gzip -c > "$destination_fasta_gz"

    gzip -t "$destination_fasta_gz"

    local record_count
    record_count="$(gzip -cd "$destination_fasta_gz" | grep -c '^>')"

    if [[ "$record_count" -ne 1 ]]; then
        echo "ERROR: $destination_fasta_gz contains $record_count records." >&2
        exit 1
    fi
}

# ============================================================
# 1. Download the two genome assemblies and assembly reports
# ============================================================

download_file "$JAGGER_URL" "$JAGGER_FASTA_GZ"
download_file "$JAGGER_REPORT_URL" "$JAGGER_REPORT"
download_file "$RYE_URL" "$RYE_FASTA_GZ"
download_file "$RYE_REPORT_URL" "$RYE_REPORT"

echo "Testing compressed reference FASTA files."
gzip -t "$JAGGER_FASTA_GZ"
gzip -t "$RYE_FASTA_GZ"

# ============================================================
# 2. Extract the chromosomes useful for inspection
#
# Rye 1R will later be appended to the complete Jagger assembly.
# Jagger 1B is extracted only as a convenient inspection copy.
# ============================================================

extract_chromosome "$RYE_FASTA_GZ" "$RYE_REPORT" "1R" "$RYE_1R_GZ"
extract_chromosome "$JAGGER_FASTA_GZ" "$JAGGER_REPORT" "1B" "$JAGGER_1B_GZ"

echo "Extracted FASTA headers:"
gzip -cd "$RYE_1R_GZ" | grep '^>'
gzip -cd "$JAGGER_1B_GZ" | grep '^>'

# ============================================================
# 3. Download and validate the Larry SRA run
# ============================================================

echo "Downloading Larry run $LARRY_RUN."
prefetch \
    --max-size u \
    --output-directory "$SRA_CACHE" \
    "$LARRY_RUN"

echo "Validating downloaded SRA data."
vdb-validate "$SRA_CACHE/$LARRY_RUN"

# ============================================================
# 4. Convert Larry to paired FASTQ and compress the files
# ============================================================

LARRY_R1_FASTQ="$FASTQ_DIR/${LARRY_RUN}_1.fastq"
LARRY_R2_FASTQ="$FASTQ_DIR/${LARRY_RUN}_2.fastq"
LARRY_R1_GZ="${LARRY_R1_FASTQ}.gz"
LARRY_R2_GZ="${LARRY_R2_FASTQ}.gz"

if [[ -s "$LARRY_R1_GZ" && -s "$LARRY_R2_GZ" ]]; then
    echo "Compressed Larry FASTQ files already exist; skipping conversion."
else
    echo "Converting Larry SRA data to paired FASTQ."
    fasterq-dump \
        --split-files \
        --threads "$THREADS" \
        --temp "$TEMP_DIR" \
        --outdir "$FASTQ_DIR" \
        --progress \
        "$SRA_CACHE/$LARRY_RUN"

    if [[ ! -s "$LARRY_R1_FASTQ" || ! -s "$LARRY_R2_FASTQ" ]]; then
        echo "ERROR: Expected paired Larry FASTQ files were not created." >&2
        echo "Inspect the SRA layout with:" >&2
        echo "  sra-stat --quick $SRA_CACHE/$LARRY_RUN" >&2
        exit 1
    fi

    echo "Compressing Larry FASTQ files."
    pigz --processes "$THREADS" "$LARRY_R1_FASTQ" "$LARRY_R2_FASTQ"
fi

echo "Testing compressed Larry FASTQ files."
gzip -t "$LARRY_R1_GZ"
gzip -t "$LARRY_R2_GZ"

# ============================================================
# 5. Final summary
# ============================================================

echo
echo "Download completed successfully: $(date)"
echo
echo "Complete Jagger assembly:"
echo "  $JAGGER_FASTA_GZ"
echo "Jagger chromosome 1B inspection copy:"
echo "  $JAGGER_1B_GZ"
echo "Complete Weining rye assembly:"
echo "  $RYE_FASTA_GZ"
echo "Weining rye chromosome 1R:"
echo "  $RYE_1R_GZ"
echo "Larry paired FASTQ files:"
echo "  $LARRY_R1_GZ"
echo "  $LARRY_R2_GZ"
echo
echo "Keep the complete Jagger assembly and Weining chromosome 1R"
echo "for construction of the later competitive-alignment reference."
