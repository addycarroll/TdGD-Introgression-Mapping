#!/bin/bash
#SBATCH --time=72:00:00           # Maximum walltime (72 hours)
#SBATCH --mem=64G                 # Memory allocation
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks-per-node=1       # Number of tasks per node
#SBATCH --job-name=bwaindex       # Job name
#SBATCH --output=slurm-%j.out      # Output file (job ID will replace %j)
#SBATCH --error=slurm-%j.err       # Error file

# Example usage:
#   sbatch indexbwaAuto.sh /path/to/ReferenceFile.fna /path/to/output_directory
# Note: Adjust the file extension in the basename command if necessary.

# Get the reference genome file path from the first argument.
RefFastaPath=$1

# Get the output directory from the second argument.
OutDir=$2

# If no output directory is provided, default to the current directory.
if [ -z "$OutDir" ]; then
  OutDir="."
fi

# Create the output directory if it doesn't exist.
mkdir -p "$OutDir"

# Extract the base name of the reference file without the extension.
# Change '.fna' to your file's extension if necessary.
refBase=$(basename -s .fna "$RefFastaPath")

# Load the BWA module (if your cluster uses modules)
module load BWA
module load SAMtools

# Run the BWA indexing command:
# - '-a bwtsw' selects the indexing algorithm (recommended for large genomes)
# - '-p' sets the prefix for the output index files, combining the output directory and base name.
bwa index -a bwtsw -p "${OutDir}/${refBase}" "$RefFastaPath"

# FASTA index (.fai) with samtools
samtools faidx "$RefFastaPath"

# (optional) copy fai into OutDir if you want everything together
cp "${RefFastaPath}.fai" "${OutDir}/${refBase}.fai"
