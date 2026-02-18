#!/usr/bin/env bash
#SBATCH --job-name=extract_gt
#SBATCH --output=extract_gt_%A_%a.out
#SBATCH --error=extract_gt_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21

# Load BCFtools
module load BCFtools

# File list with full paths to the VCFs
VCF_LIST="/bulk/akf/Addy/Variant_Calling/AcrossRun_10DP60_QUAL30_MISS_VCF/vcf_files_all_10DP60_QUAL30_MISS.txt"
OUT_DIR="/bulk/akf/Addy/Introgression_mapping/Genotype_matrices"

# Get file name for this task
VCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$VCF_LIST")
BASE=$(basename "$VCF" .vcf.gz)
CHR=$(echo "$BASE" | grep -oE 'chr[1-7][ABD]')

# Define output file
OUTFILE="$OUT_DIR/${BASE}_alleleGenos.tsv"

# Get sample names for header
HEADER=$(bcftools query -l "$VCF" | paste -sd '\t' -)
echo -e "Marker\tREF\tALT\t$HEADER" > "$OUTFILE"

# Extract genotypes and reformat the first column to Marker format
bcftools query \
  -f '%POS\t%REF\t%ALT[\t%TGT]\n' "$VCF" | \
  awk -v chr="$CHR" 'BEGIN {OFS="\t"} { $1 = chr "_" $1; print }' \
  >> "$OUTFILE"
