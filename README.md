# TdGD Introgression Mapping - Pipeline Description
## Retrieve sequencing data
**1. download_ncbi.sh:** Download SRA files from NCBI
- INPUT: Hard-coded list of SRA ids
- OUTPUT: SRA files
- Uses sratoolkit activated within a conda environment to run prefetch
***
**2. validate_sra.sh:** Validate SRA file integrity
- INPUT: A single SRA file
- OUTPUT: Console output/logs confirming file validation
- Uses sratoolkit activated within a conda environment to run vdb-validate
## Convert to FASTQ and prepare for alignment
**3. convert_single_sra.sh:** Convert a single SRA file to FASTQ format
- INPUT: A single SRA file specified as a command-line argument
- OUTPUT: Paired-end FASTQ files
- Uses sratoolkit activated within a conda environment to run fasterq-dump
***
**4. fastp_QC_single.sh:** Run quality control on paired-end FASTQ files
- INPUT: Paired-end FASTQ files from the input directory specified via a sample prefix command-line argument
- OUTPUT:
  - Quality-controlled FASTQ files (trimmed/filtered paired-end files and unpaired reads) in the output directory
  - An HTML report summarizing the QC results
- Uses fastp installed in a conda environment
***
**5. bwaIndex.sh:** Generate BWA index files for a reference genome
- INPUT: Reference genome file (.fna)
- OUTPUT: BWA index files with various extensions stored in the specified output directory
- Uses BWA to build the reference index using the bwtsw indexing algorithm
