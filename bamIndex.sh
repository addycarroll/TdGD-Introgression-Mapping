#!/bin/bash
#SBATCH --time=72:00:00           # Maximum walltime (72 hours)
#SBATCH --mem=64G                 # Memory allocation
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks per node
#SBATCH --job-name=bamindex       # Job name
#SBATCH --output=slurm-%j.out      # Output file (job ID will replace %j)
#SBATCH --error=slurm-%j.err       # Error file

# Usage:
#   sbatch indexBamFile.sh sample.bam
# This script uses fixed directories:
#   Input Directory:  /bulk/akf/Addy/Variant_Calling/bam_files
#   Output Directory: /bulk/akf/Addy/Variant_Calling/bam_files/Indexing
module load SAMtools

# Get the BAM file name from the command-line argument.
bam_filename="$1"
if [ -z "$bam_filename" ]; then
    echo "Usage: sbatch indexBamFile.sh <bam_file_name>"
    exit 1
fi

# Fixed directories.
INPUT_DIR="/bulk/akf/Addy/Variant_Calling/bam_files/WholeGenome_LowCoverage"
OUTPUT_DIR="/bulk/akf/Addy/Variant_Calling/bam_files/WholeGenome_LowCoverage"


# Construct the full path to the BAM file.
bam_file="${INPUT_DIR}/${bam_filename}"
if [ ! -f "$bam_file" ]; then
    echo "Error: File '${bam_file}' does not exist."
    exit 1
fi

echo "Indexing: ${bam_file}"
# Run the Samtools index command.
samtools index -c "$bam_file" "${OUTPUT_DIR}/${bam_filename}.csi"
