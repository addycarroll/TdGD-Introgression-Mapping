#!/bin/bash
#SBATCH --job-name=fastp_QC_single            # Job name
#SBATCH --output=fastp_QC_single_%j.out         # Output file (with job ID)
#SBATCH --error=fastp_QC_single_%j.err          # Error file
#SBATCH --time=24:00:00                        # Wall time (adjust if needed)
#SBATCH --mem=64G                              # Memory (adjust depending on node resources)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16                   # Number of tasks (threads) per node; adjust based on your node

# Load the conda environment that has fastp installed.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastp

# Define the input directory containing the FASTQ files.
INPUT_DIR="/bulk/akf/Addy/WE_fastq_fromEduard"

# Define the output directory where processed files will be written.
OUTPUT_DIR="/bulk/akf/Addy/fastq_QC"

# Check that a sample prefix was provided as a command-line argument.
if [ -z "$1" ]; then
  echo "Usage: sbatch $0 <sample_prefix>"
  echo "Example: sbatch $0 SRR15101990"
  exit 1
fi

# Set the sample prefix from the first command-line argument.
sample="$1"
echo "Processing sample: $sample"

# Construct full paths to input FASTQ files.
IN1="$INPUT_DIR/${sample}_1.fastq.gz"
IN2="$INPUT_DIR/${sample}_2.fastq.gz"

# Construct output file names.
OUT1="$OUTPUT_DIR/${sample}_1.fq.gz"
OUT2="$OUTPUT_DIR/${sample}_2.fq.gz"
UNP1="$OUTPUT_DIR/${sample}_unpaired_1.fq.gz"
UNP2="$OUTPUT_DIR/${sample}_unpaired_2.fq.gz"
HTML_REPORT="$OUTPUT_DIR/${sample}_fastp.html"

# Check that both input files exist before processing.
if [[ -f "$IN1" && -f "$IN2" ]]; then echo "Running fastp on $sample..."; fastp --in1 "$IN1" --in2 "$IN2" -f 6 -q 30 --detect_adapter_for_pe --out1 "$OUT1" --out2 "$OUT2" --unpaired1 "$UNP1" --unpaired2 "$UNP2" -h "$HTML_REPORT" --thread 16; echo "Finished processing sample: $sample"; else echo "Warning: Input files for sample $sample not found. Skipping."; fi
echo "Processing complete."
