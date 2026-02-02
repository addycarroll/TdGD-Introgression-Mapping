#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=124G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --job-name=align

# Load necessary modules
module load BWA
module load picard
module load SAMtools

# Set directories
REF="/bulk/akf/Addy/Reference_Genomes/data/Indexing/iwgsc_refseqv1.0_genomic"  # Adjust the file name if needed
FASTQ_DIR="/bulk/akf/Addy/fastq_QC"
OUT_DIR="/bulk/akf/Addy/Alignment"

# Command-line arguments:
#   $1: Base name (also used as sample name)
#   $2...: One or more mapping quality thresholds (optional; default: 40)
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <base_name> [mapq_threshold1 mapq_threshold2 ...]"
  exit 1
fi

BASE=$1
SAMPLE=$1

# Build an array of MQ thresholds. If none are given, default to 40.
if [ "$#" -eq 1 ]; then
  MQ_ARRAY=(40)
else
  shift  # Remove the base name argument so that $@ now contains only MQ thresholds.
  MQ_ARRAY=("$@")
fi
echo "Running alignment for base/sample: $BASE"
echo "Mapping quality thresholds to test: ${MQ_ARRAY[*]}"

# Change to the directory containing fastq files
cd "$FASTQ_DIR" || { echo "Failed to cd to $FASTQ_DIR"; exit 1; }

# Run BWA mem
echo "Running BWA mem..."
bwa mem -t 24 -a -M -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" "$REF" "${BASE}_1.fq.gz" "${BASE}_2.fq.gz" > "${OUT_DIR}/${BASE}.sam"

# Sort the SAM file using Picard SortSam
echo "Sorting SAM file..."
java -Xmx120g -jar $EBROOTPICARD/picard.jar SortSam --INPUT "${OUT_DIR}/${BASE}.sam" --OUTPUT "${OUT_DIR}/${BASE}_sorted.bam" --SORT_ORDER coordinate

# Mark duplicates using Picard MarkDuplicates
echo "Marking duplicates..."
java -Xmx120g -jar $EBROOTPICARD/picard.jar MarkDuplicates --INPUT "${OUT_DIR}/${BASE}_sorted.bam" --OUTPUT "${OUT_DIR}/${BASE}_sorted_marked.bam" --METRICS_FILE "${OUT_DIR}/${BASE}_dup_metrics.txt"
echo "Initial alignment, sorting, and duplicate marking completed"

# Iterate over each mapping quality threshold and run quality filtering with SAMtools
for MQ in "${MQ_ARRAY[@]}"; do
    echo "Filtering with mapping quality threshold: $MQ"
    OUTPUT_BAM="${OUT_DIR}/${BASE}_sorted_marked_q${MQ}.bam"
    samtools view -f 0x2 -F 0x808 -F 0x400 -q "$MQ" -@ 12 "${OUT_DIR}/${BASE}_sorted_marked.bam" -o "$OUTPUT_BAM"
    echo "Created filtered BAM: $OUTPUT_BAM"
done
echo "Quality filtering iterations completed"
