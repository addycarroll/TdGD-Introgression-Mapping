#!/usr/bin/env bash
#SBATCH --job-name=qual_miss_filter
#SBATCH --output=qual_miss_%A_%a.out
#SBATCH --error=qual_miss_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-21

module load BCFtools

# Directories
IN_DIR="/bulk/akf/Addy/Variant_Calling/AcrossRun_10DP60_VCF"
OUT_DIR="/bulk/akf/Addy/Variant_Calling/AcrossRun_10DP60_QUAL30_MISS_VCF"

# List of chromosomes
mapfile -t VCFS < <(ls "$IN_DIR"/Parents_chr*_D250_DP10_60.vcf.gz | sort)
IN_VCF="${VCFS[$SLURM_ARRAY_TASK_ID-1]}"
base=$(basename "$IN_VCF" .vcf.gz)

# Extract chromosome name from file
chrom=$(echo "$base" | awk -F'_' '{print $2}')

# Define genomes based on chromosome
if [[ "$chrom" =~ D$ ]]; then
    MISS_THRESH=0.50  # Retain sites with ≥2 out of 4 samples (50% missingness)
else
    MISS_THRESH=0.555  # Retain sites with ≥4 out of 9 samples (~55.5% missingness)
fi

# Apply QUAL and missingness filtering
bcftools view -i "QUAL>=30 && F_MISSING<=$MISS_THRESH" "$IN_VCF" -Oz -o "$OUT_DIR/${base}_qual_miss_filt.vcf.gz"
bcftools index "$OUT_DIR/${base}_qual_miss_filt.vcf.gz"

# Completion message
echo "[$(date)] Completed filtering for chromosome $chrom"
