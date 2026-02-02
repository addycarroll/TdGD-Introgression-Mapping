#!/bin/bash
#SBATCH --job-name=download_ncbi       # Job name
#SBATCH --output=download_ncbi_%j.out   # Standard output file
#SBATCH --error=download_ncbi_%j.err    # Standard error file
#SBATCH --time=06:00:00                 # Maximum wall time (adjust as needed)
#SBATCH --mem=8G                        # Memory per node (adjust as needed)

# Load shell configuration and activate conda environment
source ~/.bashrc
conda activate sratoolkit

# Array of SRA accessions
accessions=(
    "SRR15129972" "SRR15129973" "SRR15129974" "SRR15129975"
    "SRR15101982" "SRR15101996" "SRR15101990"
    "SRR8157374" "SRR8157382" "SRR8157418"
)

for acc in "${accessions[@]}"; do
    echo "Processing ${acc}..."

    # Remove any existing directory for this accession to avoid conflicts
    if [ -d "${acc}" ]; then
        echo "Removing existing directory ${acc}..."
        rm -rf "${acc}"
    fi
    echo "Downloading ${acc}..."

    # Use --max-size to allow for large files and --output-directory . to save in current directory
    prefetch --max-size 70G --output-directory . ${acc}

    # Check if the .sra file exists in the current directory
    if [ -f "${acc}.sra" ]; then
        echo "Converting ${acc}.sra to FASTQ format..."
        fasterq-dump ${acc}.sra --split-files
    else
        echo "Error: ${acc}.sra not found. Skipping conversion."
    fi
done
echo "All downloads and conversions are complete."
