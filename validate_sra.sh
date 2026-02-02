#!/bin/bash
#SBATCH --job-name=validate_SRA             # Job name
#SBATCH --output=validate_SRA_%j.out        # Standard output file (includes job ID)
#SBATCH --error=validate_SRA_%j.err         # Standard error file
#SBATCH --time=06:00:00                     # Time limit (adjust as needed)
#SBATCH --mem=4G                            # Memory allocation (adjust as needed)

# Source the conda initialization script so that 'conda' is available.
# Adjust the path if your conda installation is elsewhere.
source ~/miniconda3/etc/profile.d/conda.sh

# Activate the environment that contains the SRA Toolkit.
conda activate sratoolkit

# Define the SRA file you want to validate.
# You can either hard-code the path or pass it as an argument.
SRA_FILE="/bulk/akf/Addy/NCBI_SRA_files/SRR8157382.sra"

# Run vdb-validate on the specified SRA file.
vdb-validate "$SRA_FILE"
echo "Validation complete for $SRA_FILE."
