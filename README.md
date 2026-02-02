# TdGD Introgression Mapping - Pipeline Description
### Retrieve sequencing data
> 
**1. download_ncbi.sh:** Download SRA files from NCBI
- INPUT: Hard-coded list of SRA ids
- OUTPUT: SRA files
- Uses sratoolkit activated within a conda environment to run prefetch 
