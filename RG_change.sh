#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=RG
#SBATCH --out=RG-%j.out
#SBATCH --error=RG-%j.err

module load picard
# Usage:
#   sbatch process_RG.sh <file_identifier> <uniform_RG>
# Example:
#   sbatch process_RG.sh SRR15129972 KS061406LN~26

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file_identifier> <uniform_RG>"
    exit 1
fi

# Get command-line arguments.
file_id=$1
uniform_rg=$2

# Set fixed input and output directories.
input_dir="/bulk/akf/Addy/Alignment"
output_dir="/bulk/akf/Addy/Variant_Calling/bam_files"

echo "Processing files with identifier: $file_id"
echo "Uniform read group value: $uniform_rg"
echo "Input directory: $input_dir"
echo "Output directory: $output_dir"

# Loop over all BAM files matching the pattern in the input directory.
for bam in "$input_dir"/${file_id}_sorted_marked_q50.bam; do
    if [ ! -f "$bam" ]; then
        echo "No file matching pattern ${file_id}_sorted_marked_q*.bam found in $input_dir."
        exit 1
    fi
    echo "Processing file: $bam"
    # Get the base name (remove .bam extension)
    base_name=$(basename "$bam" .bam)
    # Construct the output filename (placed in the output directory)
    output_bam="${output_dir}/${base_name}_rg.bam"
    echo "Updating read group to uniform value: $uniform_rg"
    java -Xmx4g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups --INPUT "$bam" --OUTPUT "$output_bam" --RGID "$uniform_rg" --RGLB "$uniform_rg" --RGPL ILLUMINA --RGPU "${uniform_rg}_unit" --RGSM "$uniform_rg" --VALIDATION_STRINGENCY SILENT
    echo "Created updated file: $output_bam"
done
echo "All files processed."

