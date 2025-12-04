#!/bin/bash
# Bash script for downloading genome sequencing reads from NCBI Read Sequencing Archive (SRA)
#
# Author: Marcus Vinicius Canário Viana
# Date: 04/12/2025
# Repository: https://github.com/canarioviana/bacterial_genome_assembly
# More info: see README.md in the repository
#
# Instructions:
#
# **A. Create the input file**
#
# 1. Create a **tab-separated file** named **"1_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"1_reads_accessions.tsv"** file in the working directory.
#
# **B. Execution**
#
# Place this script (**bga_download_reads_end2end.sh**) in the working directory and execute it **using the following commands**:
# chmod +x bga_download_reads_end2end.sh
# ./bga_download_reads_end2end.sh


############################################################
## 0) Error handling and checking Conda installation
############################################################

############################################################
## Error handling
# Exit immediately if a command fails (returns a non-zero exit code)
set -e
# Ensure that a pipeline (command1 | command2) fails if any command in the pipe fails
set -o pipefail

############################################################
## Check if the 'conda' command is available on the system PATH
if command -v conda &> /dev/null; then
    # Locate the base Conda installation path
    CONDA_BASE=$(conda info --base)
    # Makes 'conda activate' available in the current subshell
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "ERROR: The 'conda' command was not found" >&2
    echo "Ensure Conda or Miniconda is installed and configured in your PATH" >&2
    exit 1
fi

############################################################
## 1) Sequencing reads directory and files
############################################################

############################################################
## Reads from NCBI SRA

# Delete previous file of not used reads
rm -f 1_reads_single-end.tsv
rm -f 1_reads_paired-end.tsv

# Verify the presence of the file 1_reads_accessions.tsv with a list of accessions
if [ -f 1_reads_accessions.tsv ]; then
    echo "The file 1_reads_accessions.tsv was found. The sequencing reads will be downloaded."

    # Create output directory 
    mkdir -p 1_reads

    # SRA Tools
    # Activate Conda environment
    conda activate sra-tools

    # Loop through file lines
    while IFS=$'\t' read -r accession sample others; do
        if [ -f "1_reads/${sample}_1.fq.gz" ] && [ -f "1_reads/${sample}_2.fq.gz" ]; then
            echo "Sample $sample paired files found. Skipping download."
        elif [ -f "1_reads/${sample}.fq.gz" ]; then
            echo "Sample $sample single-end file found. Skipping download."
        else
            echo "Downloading sample: $sample (accession: $accession)"

            # Remove any incomplete files
            rm -f "1_reads/${sample}_1.fq.gz" "1_reads/${sample}_2.fq.gz" "1_reads/${sample}.fq.gz"    
            
            # Run prefetch
            prefetch -p -O 1_reads "${accession}"

            # Run fasterq-dump
            cd 1_reads
            fasterq-dump \
            --threads $(nproc --ignore=1) \
            -p \
            --split-files "${accession}" \
            -O .
            # Delete temporary directories
            rm -r "${accession}"
            # Compress files
            echo "Compressing fastq files."
            pigz -p $(nproc --ignore=1) ${accession}*.fastq

            # Check if the reads are paired-end
            if [ -f "${accession}_1.fastq.gz" ] && [ -f "${accession}_2.fastq.gz" ];  then
                echo "Sample ${sample} has paired-end reads."
                # Rename files
                echo "Renaming files."
                mv "${accession}_1.fastq.gz" "${sample}_1.fq.gz"
                mv "${accession}_2.fastq.gz" "${sample}_2.fq.gz"
                # Add file to the list of single-end libraries
                echo -e "${accession}\t${sample}_1.fq.gz\t${sample}_2.fq.gz" >> ../1_reads_paired-end.tsv
                # Go back to main directory
                cd ..
            # Check if the reads are not paired-end
            elif [ -f "${accession}.fastq.gz" ]; then
                echo "Sample ${sample} has single-end reads."
                # Renaming file
                echo "Renaming file."
                mv "${accession}.fastq.gz" "${sample}.fq.gz"
                # Add file to the list of single-end libraries
                echo -e "${accession}\t${sample}.fq.gz" >> ../1_reads_single-end.tsv
                # Go back to main directory
                cd ..
            fi
        fi
    done < 1_reads_accessions.tsv
    echo "Download process complete. Deactivating the environment."
    # Deactivate Conda environment
    conda deactivate
else
    echo "The file 1_reads_accessions.tsv was not found."
    exit 1
fi
