#!/bin/bash
#SBATCH --job-name=convert_sra_single       # Job name
#SBATCH --output=convert_sra_single_%j.out    # Standard output file (with job ID)
#SBATCH --error=convert_sra_single_%j.err     # Standard error file
#SBATCH --time=24:00:00                       # Maximum wall time (adjust if needed)
#SBATCH --mem=8G                              # Memory per node (adjust if needed)

# Source the conda initialization script so that 'conda' is available
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the sratoolkit environment
conda activate sratoolkit

# Get the SRA file from the command-line argument
SRA_FILE="$1"

# Check if an SRA file was provided
if [ -z "$SRA_FILE" ]; then
    echo "Error: No SRA file specified."
    echo "Usage: sbatch convert_single_sra.sh /path/to/file.sra"
    exit 1
fi
echo "Starting conversion of SRA file '$SRA_FILE'..."

# Run fasterq-dump on the specified SRA file, splitting paired-end reads if present
fasterq-dump "$SRA_FILE" --split-files
echo "Conversion complete for $SRA_FILE."
