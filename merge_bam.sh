#!/bin/bash
#SBATCH --job-name=merge_srr
#SBATCH --output=merge_srr_%j.out
#SBATCH --error=merge_srr_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load SAMtools
INPUT_DIR="/bulk/akf/Addy/Alignment/RG"
OUTPUT_DIR="/bulk/akf/Addy/Variant_Calling/bam_files"

samtools merge \
    "$OUTPUT_DIR"/SRR1512997x_merged.bam \
    "$INPUT_DIR"/SRR15129972_sorted_marked_q50_rg.bam \
    "$INPUT_DIR"/SRR15129973_sorted_marked_q50_rg.bam \
    "$INPUT_DIR"/SRR15129974_sorted_marked_q50_rg.bam \
    "$INPUT_DIR"/SRR15129975_sorted_marked_q50_rg.bam

echo "Merging complete. Output file:"
echo "$OUTPUT_DIR/SRR1512997x_merged.bam"
