#!/usr/bin/env bash
#SBATCH --job-name=filter_DP10_60_array
#SBATCH --output=slurm_array_%A_%a.out
#SBATCH --error=slurm_array_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21

module load BCFtools

VCF_LIST="/bulk/akf/Addy/Variant_Calling/vcf_files/vcf_files_all.txt"
OUT_DIR="/bulk/akf/Addy/Variant_Calling/AcrossRun_10DP60_VCF"


# Pick the VCF for this task
in_vcf=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$VCF_LIST")
base=$(basename "$in_vcf" .vcf.gz)
out_vcf="$OUT_DIR/${base}_DP10_60.vcf.gz"
echo "→ [Task ${SLURM_ARRAY_TASK_ID}] Filtering $base (DP ≥10 & ≤60)…"

bcftools view \
  -v snps \
  -m2 -M2 \
  -i 'INFO/DP>=10 && INFO/DP<=60' \
  -Oz \
  -o "$out_vcf" \
  "$in_vcf"


bcftools index "$out_vcf"

