#!/bin/bash
#SBATCH --job-name=MQ_comparison
#SBATCH --output=MQ_comparison_%j.out
#SBATCH --error=MQ_comparison_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=124G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24

# Load necessary modules (adjust as needed)
module load SAMtools
module load BamTools

# Check for sample ID input
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 SAMPLE_ID"
    exit 1
fi

SAMPLE_ID=$1

# Define input and output directories
INPUT_DIR="/bulk/akf/Addy/Alignment"
OUTPUT_DIR="/bulk/akf/Addy/Alignment/MQ_comparison"

# Use shell globbing to find BAM files matching the pattern in the input directory
# Expected pattern: SAMPLEID_sorted_marked_q*.bam
BAM_FILES=( ${INPUT_DIR}/${SAMPLE_ID}_sorted_marked_q*.bam )

# Check if any files were found
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found matching: ${INPUT_DIR}/${SAMPLE_ID}_sorted_marked_q*.bam"
    exit 1
fi

# Loop through each BAM file and process
for BAM in "${BAM_FILES[@]}"; do
    echo "Processing ${BAM}..."
    # Extract the base filename without directory and extension
    BASENAME=$(basename "${BAM}" .bam)
    # Define output filenames in the output directory
    STATS_OUT="${OUTPUT_DIR}/${BASENAME}_stats.txt"
    DEPTH_OUT="${OUTPUT_DIR}/${BASENAME}_depth.txt"
    AVG_COVERAGE_OUT="${OUTPUT_DIR}/${BASENAME}_avg_coverage.txt"
    DISTRO_OUT="${OUTPUT_DIR}/${BASENAME}_coverage_distribution.txt"
    # Run bamtools stats for alignment statistics
    echo "Running bamtools stats on ${BAM}..."
    bamtools stats -in "${BAM}" > "${STATS_OUT}"
    # Use samtools depth to calculate per-base coverage
    echo "Running samtools depth on ${BAM}..."
    samtools depth "${BAM}" > "${DEPTH_OUT}"
    # Calculate average coverage from the samtools depth output
    echo "Calculating average coverage for ${BAM}..."
    awk '{total+=$3; count++} END {if (count > 0) print "Average coverage =", total/count; else print "No coverage data"}' "${DEPTH_OUT}" > "${AVG_COVERAGE_OUT}"
    # Calculate coverage distribution (number of bases at each coverage level)
    echo "Calculating coverage distribution for ${BAM}..."
    awk '{cov[$3]++} END {for (level in cov) print level, cov[level]}' "${DEPTH_OUT}" | sort -n > "${DISTRO_OUT}"
done
echo "Analysis complete! Output files are located in ${OUTPUT_DIR}"
