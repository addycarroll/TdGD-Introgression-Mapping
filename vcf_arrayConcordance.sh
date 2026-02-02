#!/usr/bin/env bash
#SBATCH --job-name=arraysites_call_D250
#SBATCH --output=arraysites_call_D250_%j.out
#SBATCH --error=arraysites_call_D250_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load BCFtools

module load SAMtools

SITES_DIR="/bulk/akf/Addy/Variant_Calling/Array_concordance"

TARGETS="${SITES_DIR}/array_sites.pos.tsv"   # lines like LS992080.1:1207866-1207866

TMP_BAM_LIST="${SITES_DIR}/bam_list.txt"

REF_GENOME="/bulk/akf/Addy/Reference_Genomes/data_CSv1/Indexing/iwgsc_refseqv1.0_genomic.fna"

MAX_DEPTH=250

OUT_PATH="${SITES_DIR}"

OUT_VCF="${OUT_PATH}/Parents_arraySites_D${MAX_DEPTH}.vcf.gz"

bcftools mpileup \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --annotate FORMAT/AD,FORMAT/DP \
  --skip-indels \
  -f "$REF_GENOME" \
  -R "$TARGETS" \
  -b "$TMP_BAM_LIST" \
  -B \
  -d "$MAX_DEPTH" \
  -Ou \
| bcftools call \
    --threads "${SLURM_CPUS_PER_TASK:-8}" \
    -m \
    --variants-only \
    --skip-variants indels \
    -Oz \
    -o "$OUT_VCF"

bcftools index -t "$OUT_VCF"
echo "[DONE] $OUT_VCF"
