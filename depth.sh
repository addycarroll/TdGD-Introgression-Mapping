#!/bin/bash
#SBATCH --job-name=depth_max
#SBATCH --output=depth_max_%A_%a.out
#SBATCH --error=depth_max_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --array=0-13    # Adjust the upper bound to (number_of_depth_files - 1)

# Directory containing your depth files
DEPTH_DIR="/bulk/akf/Addy/Alignment/MQ_comparison/MQ50"

# Create an array of all depth files in the directory.
# This example assumes your depth files end with "_depth.txt"
depth_files=($(find "$DEPTH_DIR" -maxdepth 1 -type f -name "*_depth.txt" | sort))

# Check that the array is not empty
if [ ${#depth_files[@]} -eq 0 ]; then
    echo "No depth files found in ${DEPTH_DIR}"
    exit 1
fi

# Select the file corresponding to this job array index
file=${depth_files[$SLURM_ARRAY_TASK_ID]}
echo "Processing file: ${file}"

# Compute the maximum depth from the third column
max_depth=$(awk '{if($3>max) max=$3} END {print max}' "$file")
echo "File: ${file}"
echo "Max depth: ${max_depth}"
