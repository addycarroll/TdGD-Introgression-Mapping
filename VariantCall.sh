#!/bin/bash
#SBATCH --job-name=variant_call       # Job name
#SBATCH --output=variant_call_%A_%a.out   # Standard output file
#SBATCH --error=variant_call_%A_%a.err    # Standard error file
#SBATCH --time=72:00:00                # Maximum wall time (adjust if needed)
#SBATCH --mem=8G                        # Memory per node (adjust if needed)
#SBATCH --cpus-per-task=4               # Use 4 CPUs
#SBATCH --array=1-21                    # Job array for 21 chromosomes

module load BCFtools
module load SAMtools

# Usage check: Provide subgenome as argument (A, B, or D)
if [ -z "$1" ]; then
    echo "Usage: sbatch variant_calling.sh <subgenome>"
    echo "Example: sbatch variant_calling.sh A"
    exit 1
fi
SUBGENOME=$1

# Define mapping arrays for each subgenome
refseqs_A=( "LS992080.1" "LS992083.1" "LS992086.1" "LS992089.1" "LS992092.1" "LS992095.1" "LS992098.1" )
refseqs_B=( "LS992081.1" "LS992084.1" "LS992087.1" "LS992090.1" "LS992093.1" "LS992096.1" "LS992099.1" )
refseqs_D=( "LS992082.1" "LS992085.1" "LS992088.1" "LS992091.1" "LS992094.1" "LS992097.1" "LS992100.1" )

# Define depth thresholds to test
DEPTHS=(250 500 1000)

# Determine chromosome index and depth index from SLURM_ARRAY_TASK_ID
CHROM_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % 7 ))  # 7 chromosomes per subgenome
DEPTH_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / 7 ))  # Each chromosome runs for each depth

# Select chromosome based on subgenome
if [ "$SUBGENOME" == "A" ]; then
    REGION=${refseqs_A[$CHROM_INDEX]}
elif [ "$SUBGENOME" == "B" ]; then
    REGION=${refseqs_B[$CHROM_INDEX]}
elif [ "$SUBGENOME" == "D" ]; then
    REGION=${refseqs_D[$CHROM_INDEX]}
else
    echo "Invalid subgenome specified. Use A, B, or D."
    exit 1
fi

# Select depth threshold
MAX_DEPTH=${DEPTHS[$DEPTH_INDEX]}
echo "Processing region: ${REGION} with depth threshold: ${MAX_DEPTH}"

# Directory where your BAM files are stored.
BAM_DIR="/bulk/akf/Addy/Variant_Calling/bam_files"

# Set the appropriate BAM list file.
if [ "$SUBGENOME" == "A" ] || [ "$SUBGENOME" == "B" ]; then
    BAM_LIST="/bulk/akf/Addy/Variant_Calling/bam_list_AB.txt"
elif [ "$SUBGENOME" == "D" ]; then
    BAM_LIST="/bulk/akf/Addy/Variant_Calling/bam_list_D.txt"
fi

# Set the output directory for VCF files.
OUT_PATH="/bulk/akf/Addy/Variant_Calling/vcf_files"

# Path to your reference genome.
REF_GENOME="/bulk/akf/Addy/Reference_Genomes/data/Indexing/iwgsc_refseqv1.0_genomic.fna"

echo "Using BAM list file: ${BAM_LIST}"
echo "Assuming BAM files are located in: ${BAM_DIR}"

# Prepare a temporary BAM list file with full paths.
TMP_BAM_LIST=$(mktemp)
while IFS= read -r line; do
    if [[ $line == /* ]]; then
        echo "$line" >> "$TMP_BAM_LIST"
    else
        echo "${BAM_DIR}/${line}" >> "$TMP_BAM_LIST"
    fi
done < "$BAM_LIST"

# Run bcftools mpileup and call variants with different depth thresholds
bcftools mpileup --threads 10 --annotate AD,DP,INFO/AD --skip-indels -f "$REF_GENOME" -r "$REGION" -b "$TMP_BAM_LIST" -B -d "$MAX_DEPTH" | bcftools call --threads 10 -m --variants-only --skip-variants indels --output-type v -o "${OUT_PATH}/Parents_${REGION}_D${MAX_DEPTH}.vcf" --group-samples -

# Clean up temporary file.
rm "$TMP_BAM_LIST"
