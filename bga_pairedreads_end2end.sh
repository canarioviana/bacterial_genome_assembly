#!/bin/bash
# Bash script for bacterial genome assembly from short-read sequencing data (end-to-end worflow)
#
# Author: Marcus Vinicius Canário Viana
# Date: 14/09/2026
# Repository: https://github.com/canarioviana/bacterial_genome_assembly
# More info: see README.md in the repository
#
# Instructions:
#
# **A. Using Local Read Files**
#
# 1. Standardize the paired-end file names of each sample to **samplename_1.fq.gz** and **samplename_2.fq.gz**.
# 2. In the working directory, create the directory **1_reads** and place the read files **inside it**.
#
# **B. Downloading Reads from NCBI SRA**
#
# 1. Create a **tab-separated file** named **"0_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"0_reads_accessions.tsv"** file in the working directory.
#
# **C. Execution**
#
# Place this script (**bga_pairedreads_end2end.sh**) in the working directory and execute it **using the following commands**:
# chmod +x bga_pairedreads_end2end.sh
# ./bga_pairedreads_end2end.sh


############################################################
## SUMMARY OF END-TO-END GENOME ASSEMBLY WORKFLOW FROM SHORT-READS
############################################################

## 0) Error handling and checking Conda installation
## 1) Sequencing reads directory and files
    # Download reads from NCBI SRA (sra-tools)
    # Check local reads files
## 2) Raw reads quality assessment
    # FastQC
    # MultiQC
## 3) Raw reads trimming, estimation of genome size and downsampling 
    # Fastp
    # Estimation of genome size (KMC and GenomeScope) and downsampling (Rasusa)
## 4) Trimmed reads quality assessment
    # FastQC
    # MultiQC
## 5) De novo assembly
    # Unicycler
## 6) Organization of de novo assembly files
## 7) Assembly quality assessment
    # CheckM2
    # GUNC
    # QUAST
    # Barrnap
    # Calculation of vertical sequencing coverage
## 8) Taxonomic assignment
    # GTDB-Tk
## 9) Plasmid identification
    # MOB-suite 
## 10) Assignment of contigs to molecules


############################################################
## Preparing input files
############################################################

# ## Metadata for sequencing reads from NCBI SRA

# 1. Create a **tab-separated file** named **"0_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"0_reads_accessions.tsv"** file in the working directory.

# ---

# ## Sequencing reads as local files

# 1. The sequencing reads must be in FASTQ format and compressed, with the suffixes `_1.fq.gz` and `_2.fq.gz`, or `_1.fastq.gz` and `_2.fastq.gz` or `_R1_001.fastq.gz` and `_R2_001.fastq.gz`
# 2. In the working directory, create the directory `1_reads/` and place the read files inside it.

# ---

############################################################
# 0) Beginning of the script - Error handling and others
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
    echo "✗  ERROR: The 'conda' command was not found." >&2
    echo "Ensure Conda or Miniconda is installed and configured in your PATH." >&2
    # Exit the script
    exit 1
fi


############################################################
##  Mark the start of a full pipeline run (fresh start or restart)

echo -e "════════════════════════════════════════════════════════
WORKFLOW STARTED @ $(date +'%Y-%m-%d %H:%M:%S')
════════════════════════════════════════════════════════\n" | tee -a 0_workflow_progress.txt

############################################################
## Check for metadata files and validate consistency

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="0) Metadata check"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt

# Check if reads metadata file '0_reads_accessions.tsv' exists and is not empty
if [ -f "0_reads_accessions.tsv" ]; then
    if [ ! -s "0_reads_accessions.tsv" ]; then
        echo "✗  ERROR: Metadata file for reads from NCBI SRA '0_reads_accessions.tsv' exists but is empty!" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    echo "✔  Metadata file for reads from NCBI SRA '0_reads_accessions.tsv' found and non-empty." | tee -a 0_workflow_progress.txt
    # Remove Windows CRLF line endings
    sed -i 's/\r$//' 0_reads_accessions.tsv
else
    echo "⚠️ Warning: The file 0_reads_accessions.tsv was not found. No attempt will be made to download sequencing reads in the next step." | tee -a 0_workflow_progress.txt
fi

# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
# 1) Reads files and renaming
############################################################

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="1) Reads files and renaming"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

############################################################
## 1.1) Reads from NCBI SRA (SRA Tools)

# Avoid literal glob pattern
shopt -s nullglob

# Delete previous file of not used reads
rm -f 1_reads_not_pe.tsv

if [ -f 0_reads_accessions.tsv ]; then
    echo "✔  The file 0_reads_accessions.tsv was found. The sequencing reads will be downloaded."

    # Create output directory
    mkdir -p 1_reads

    # Activate Conda environment
    conda activate sra-tools

    # Skip blank lines and an optional header line
    tr -d '\r' < 0_reads_accessions.tsv | \
    while IFS=$'\t' read -r accession sample others; do

        # Check if the first column is "accession"
        [ -z "$accession" ] && continue
        [ "$accession" = "accession" ] && continue

        # Declare variables to check file pairing
        paired_ok=0
        single_ok=0

        # Check if paired read files are already present
        if [ -f "1_reads/${sample}_1.fq.gz" ] && [ -f "1_reads/${sample}_2.fq.gz" ] \
            && gzip -t "1_reads/${sample}_1.fq.gz" 2>/dev/null && gzip -t "1_reads/${sample}_2.fq.gz" 2>/dev/null; then
            # Pair found
            paired_ok=1
        elif [ -f "1_reads/${sample}.fq.gz" ] && gzip -t "1_reads/${sample}.fq.gz" 2>/dev/null; then
            # Single read file found
            single_ok=1
        fi
        if [ "$paired_ok" -eq 1 ]; then
            # Inform that a valid read file pair was found
            echo "✔  Sample $sample paired files found and valid. Skipping download."
        elif [ "$single_ok" -eq 1 ]; then
            # Inform that a valid single read file was found
            echo "✔  Sample $sample single-end file found and valid. Skipping download."
            echo -e "${accession}\t${sample}.fq.gz" >> 1_reads_not_pe.tsv
        else
            # For any other case delete files and start download
            echo "▶  Downloading sample: $sample (accession: $accession)"
            rm -f "1_reads/${sample}_1.fq.gz" "1_reads/${sample}_2.fq.gz" "1_reads/${sample}.fq.gz"
            rm -rf "1_reads/${accession}"
            rm -f "1_reads/${accession}"*.fastq "1_reads/${accession}"*.fastq.gz

            # Download SRA file
            prefetch -p -O 1_reads "${accession}"

            (
                # Go to 1_reads
                cd 1_reads || exit 1

                # Create fastq files from SRA file
                fasterq-dump \
                    --threads $(nproc --ignore=1) \
                    -p \
                    --split-files "${accession}" \
                    -O .

                # Delete SRA file
                rm -rf "${accession}"

                # Compresss fastq files
                echo "Compressing fastq files."
                pigz -p $(nproc --ignore=1) ${accession}*.fastq

                if [ -f "${accession}_1.fastq.gz" ] && [ -f "${accession}_2.fastq.gz" ] \
                    && gzip -t "${accession}_1.fastq.gz" 2>/dev/null && gzip -t "${accession}_2.fastq.gz" 2>/dev/null; then
                    # Inform that a valid read file pair was downloaded
                    echo "Sample ${sample} has paired-end reads."
                    echo "Renaming files."
                    # Rename files
                    mv "${accession}_1.fastq.gz" "${sample}_1.fq.gz"
                    mv "${accession}_2.fastq.gz" "${sample}_2.fq.gz"
                    # Remove unpaired
                    rm -f "${accession}.fastq.gz"

                elif [ -f "${accession}.fastq.gz" ] && gzip -t "${accession}.fastq.gz" 2>/dev/null; then
                    # Inform that a valid single read file was downloaded
                    echo "Sample ${sample} has single-end reads."
                    echo "Renaming file."
                    # Rename file
                    mv "${accession}.fastq.gz" "${sample}.fq.gz"
                    echo "Warning: this script only uses paired-end reads. This file will not be used."
                    echo -e "${accession}\t${sample}.fq.gz" >> ../1_reads_not_pe.tsv

                else
                    # Inform that the download was not successful
                    echo "✗  ERROR: No valid output produced for sample ${sample} (accession: ${accession}). Check the accession ID"
                    rm -f "${accession}"*.fastq.gz
                    exit 1
                fi
            )
        fi
    done
    echo "✔  Download process complete."
    # Deactivate Conda environment
    conda deactivate
else
    echo "⚠️  Warning: The file 0_reads_accessions.tsv was not found. Proceeding using local files."
fi

############################################################
## 1.2) Reads stored as local files

# Checking wether the directory 1_reads exists
if [ ! -d 1_reads ]; then
    echo "✗  ERROR: The reads directory '1_reads' was not found" | tee -a 0_workflow_progress.txt
    echo "Please create it and put the files in it"
    exit 1
fi

# Avoid literal glob pattern
shopt -s nullglob

# Check and rename Illumina standard format (*_R1_001.fastq.gz to *_1.fq.gz)
r1_files=(1_reads/*_R1_001.fastq.gz)
if [ ${#r1_files[@]} -gt 0 ]; then
    echo "Found files in the format *_R1_001.fastq.gz and *_R2_001.fastq.gz" | tee -a 0_workflow_progress.txt
    echo "Renaming them to the format *_1.fq.gz and *_2.fq.gz..." | tee -a 0_workflow_progress.txt

    # Rename files
    rename 's/_R1_001\.fastq\.gz/_1.fq.gz/; s/_R2_001\.fastq\.gz/_2.fq.gz/' 1_reads/*.fastq.gz 2>/dev/null

    # Rename files inside md5 file
    r1_md5_files=(1_reads/*_R1_001.fastq.gz.md5)
    if [ ${#r1_md5_files[@]} -gt 0 ]; then
        # .md5 content references the data filename, not the .md5 filename itself
        sed -i 's/_R1_001\.fastq\.gz/_1.fq.gz/; s/_R2_001\.fastq\.gz/_2.fq.gz/' 1_reads/*.md5 2>/dev/null
        rename 's/_R1_001\.fastq\.gz\.md5/_1.fq.gz\.md5/; s/_R2_001\.fastq\.gz\.md5/_2.fq.gz\.md5/' 1_reads/*.fastq.gz.md5 2>/dev/null
    fi
fi

# Check and rename *.fastq.gz to *.fq.gz
fastq_files=(1_reads/*_1.fastq.gz)
if [ ${#fastq_files[@]} -gt 0 ]; then
    echo "Found files in the format *_1.fastq.gz and *_2.fastq.gz" | tee -a 0_workflow_progress.txt
    echo "Renaming them and their pairs to the format *_1.fq.gz and *_2.fq.gz..." | tee -a 0_workflow_progress.txt

    # Rename files
    rename 's/\.fastq\.gz$/.fq.gz/' 1_reads/*.fastq.gz 2>/dev/null

    # Rename files inside md5 file
    md5_files=(1_reads/*_1.fastq.gz.md5)
    if [ ${#md5_files[@]} -gt 0 ]; then
        # 
        sed -i 's/\_1.fastq\.gz$/\_1.fq.gz/; s/\_2.fastq\.gz$/\_2.fq.gz/' 1_reads/*.md5 2>/dev/null
        rename 's/\_1.fastq\.gz\.md5$/\_1.fq.gz\.md5/; s/\_2.fastq\.gz\.md5$/\_2.fq.gz\.md5/' 1_reads/*.fastq.gz.md5 2>/dev/null
    fi
fi

# Verify the presence, names, and pairing of input files
echo "Checking the presence, names, and pairing of input files..." | tee -a 0_workflow_progress.txt

# Declare variables
files_found=0
missing_pairs=0
fq_files=(1_reads/*_1.fq.gz)

if [ ${#fq_files[@]} -gt 0 ]; then
    files_found=1
    for r1 in "${fq_files[@]}"; do
        # Obtain r2 path
        r2="${r1/_1.fq.gz/_2.fq.gz}"
        if [ ! -f "$r2" ]; then
            echo "✗  Expected pair file '${r2}' not found!"
            missing_pairs=1
        fi
    done
fi

if [ "$files_found" -eq 0 ]; then
    echo "✗  ERROR: No read files matching 'samplename_1.fq.gz' were found in '1_reads'." | tee -a 0_workflow_progress.txt
    exit 1
fi
if [ "$missing_pairs" -eq 1 ]; then
    echo "✗  ERROR: The script requires matching read pairs (*_1.fq.gz and *_2.fq.gz) for all samples." | tee -a 0_workflow_progress.txt
    echo "Please rename or add the missing pair files into '1_reads'."
    exit 1
fi
echo "✔  All read files verified successfully and paired properly." | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
# 2) Raw reads quality assessment
############################################################

############################################################
## 2.1) FastQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) FastQC"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="2_fastqc.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (2_fastqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ -d "2_fastqc" ]; then
        echo "${workflow_step}: Found incomplete 2_fastqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "2_fastqc"
    fi
    rm -f "2_fastqc.tar.gz" "2_fastqc.tar.gz.md5"

    fq_files=(1_reads/*.fq.gz)
    if [ ${#fq_files[@]} -eq 0 ]; then
        echo "✗  ERROR: No .fq.gz files found in 1_reads/ for FastQC." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create an output directory
    mkdir -p 2_fastqc

    # Activate Conda environment
    conda activate fastqc
    # Run main software
    fastqc -t $(nproc --ignore=1) "${fq_files[@]}" -o 2_fastqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    itens_to_compress=(2_fastqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 2.2) FastQC -> MultiQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="2) MultiQC"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="2_fastqc_multiqc.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (2_fastqc_multiqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
    rm -rf "2_fastqc" 2>/dev/null
else
    if [ -d "2_fastqc_multiqc" ]; then
        echo "${workflow_step}: Found incomplete 2_fastqc_multiqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "2_fastqc_multiqc"
    fi
    rm -f "2_fastqc_multiqc.tar.gz" "2_fastqc_multiqc.tar.gz.md5"

    # Check for the presence of the input files
    fastqc_zips=(2_fastqc/*_fastqc.zip)
    if [ ${#fastqc_zips[@]} -eq 0 ]; then
        echo "✗  ERROR: No *_fastqc.zip files found in 2_fastqc/ for MultiQC." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Activate Conda environment
    conda activate multiqc
    # Run main software
    multiqc "${fastqc_zips[@]}" -o 2_fastqc_multiqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    itens_to_compress=(2_fastqc_multiqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete output directory after the compressed file is verified valid.
    rm -r 2_fastqc 2_fastqc_multiqc
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
# 3) Raw reads trimming
############################################################

############################################################
## 3.1) Fastp

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="3) Fastp"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create an output directory
mkdir -p 3_fastp    

# Calculate sample size
i=1
r1_files=(1_reads/*_1.fq.gz)
sample_count=${#r1_files[@]}
if [ "$sample_count" -eq 0 ]; then
    echo "✗  ${workflow_step}: No input found matching 1_reads/*_1.fq.gz. Check that the previous step completed successfully." | tee -a 0_workflow_progress.txt
    exit 1
fi

# Track whether any sample was actually processed
work_done=false

# Activate Conda environment
conda activate fastp
# Loop through a list of sample files
for r1 in "${r1_files[@]}"; do

    # Obtain r2 path
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%_1.fq.gz}

    # Skip sample if trimmed and downsampled output already exist and is valid
    out1="3_fastp/${sample}_trimmed_1.fq.gz"
    out2="3_fastp/${sample}_trimmed_2.fq.gz"
    ds1="3_fastp_downsampling/${sample}_trimmed_ds_1.fq.gz"
    ds2="3_fastp_downsampling/${sample}_trimmed_ds_2.fq.gz"
    if [ -f "$out1" ] && [ -f "$out2" ] && gzip -t "$out1" 2>/dev/null && gzip -t "$out2" 2>/dev/null; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$ds1" ] && [ -f "$ds2" ] && gzip -t "$ds1" 2>/dev/null && gzip -t "$ds2" 2>/dev/null; then
        echo "${workflow_step} output already consumed by downsampling step for sample: $sample. Skipping sample."
        i=$((i + 1))
        continue
    elif [ -f "$out1" ] || [ -f "$out2" ] || [ -f "3_fastp/${sample}_trimmed_fastp.html" ] || [ -f "3_fastp/${sample}_trimmed_fastp.json" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$out1" "$out2" "3_fastp/${sample}_trimmed_fastp.html" "3_fastp/${sample}_trimmed_fastp.json"
    fi

    # The sample will actually be (re)processed
    work_done=true

    # Inform current sample
    echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    echo "${workflow_step} - Read 1 file: ${r1}"
    echo "${workflow_step} - Read 2 file: ${r2}"
    # Start counting the loop running time
    loop_start_time=$SECONDS

    # Run main software
    fastp \
        --thread $(nproc --ignore=1) \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        --overrepresentation_analysis \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "$out1" \
        --out2 "$out2" \
        --html "3_fastp/${sample}_trimmed_fastp.html" \
        --json "3_fastp/${sample}_trimmed_fastp.json"
    
    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time}" | tee -a 0_workflow_progress.txt

    # Increate sample count
    i=$((i + 1))

done
# Deactivate Conda environment
conda deactivate

compressed_file="3_fastp.tar.gz"
itens_to_compress=(3_fastp/*.json 3_fastp/*.html)

# Skip recompressing/re-verifying if nothing changed and a valid archive from a previous run already exists
if [ "$work_done" = false ] && [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] \
    && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (${compressed_file} verified). Skipping compression." | tee -a 0_workflow_progress.txt
else
    # Compress the output directory
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
    fi
    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    # Check file integrity
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt
fi

# Generate checksum files for the reads
(cd 3_fastp && for file in *.gz; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity of reads" | tee -a 0_workflow_progress.txt
(cd 3_fastp && md5sum -c *.md5) | tee -a 0_workflow_progress.txt

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 3.2) Downsampling (KMC, GenomeScope and Rasusa)

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="3) Downsampling (KMC, GenomeScope and Rasusa)"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Create output directory
mkdir -p 3_fastp_downsampling

# Count and verify sample files
i=1
r1_files=(3_fastp/*_trimmed_1.fq.gz)
sample_count=${#r1_files[@]}
r1_ds_files=(3_fastp_downsampling/*_trimmed_ds_1.fq.gz)
if [ "$sample_count" -eq 0 ] && [ ${#r1_ds_files[@]} -eq 0 ]; then
    echo "✗  ${workflow_step}: No input found matching 3_fastp/*_trimmed_1.fq.gz nor existing output in 3_fastp_downsampling/. Check that the previous step completed successfully." | tee -a 0_workflow_progress.txt
    exit 1
fi

# Track whether any sample was actually processed
work_done=false

for r1 in "${r1_files[@]}"; do
    # Obtain r2 path
    r2="${r1/_trimmed_1.fq.gz/_trimmed_2.fq.gz}"
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%_trimmed_1.fq.gz}

    # Obtain paths to output files
    out1="3_fastp_downsampling/${sample}_trimmed_ds_1.fq.gz"
    out2="3_fastp_downsampling/${sample}_trimmed_ds_2.fq.gz"
    genomesizedir="3_fastp_downsampling/${sample}_genomesize"

    if [ -f "$out1" ] && [ -f "$out2" ] && gzip -t "$out1" 2>/dev/null && gzip -t "$out2" 2>/dev/null; then
        echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
        # Increate sample count
        i=$((i + 1))
        continue
    elif [ -f "$out1" ] || [ -f "$out2" ] || [ -f "3_fastp_downsampling/${sample}_trimmed_fastp.html" ] || [ -f "3_fastp_downsampling/${sample}_trimmed_fastp.json" ] \
         || [ -d "$genomesizedir" ]; then
        echo "${workflow_step} found incomplete/corrupted output for sample: $sample. Removing partial files and reprocessing."
        rm -f "$out1" "$out2" "3_fastp_downsampling/${sample}_trimmed_fastp.html" "3_fastp_downsampling/${sample}_trimmed_fastp.json"
        rm -rf "$genomesizedir" kmc_tmp kmc_count* kmc_histogram.tsv kmc_input_reads.txt
    fi

    # The sample will actually be (re)processed
    work_done=true

    echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
    loop_start_time=$SECONDS

    # Create output directory
    mkdir -p "${genomesizedir}"

    # Activate Conda environment
    conda activate kmc
    echo "Counting k-mers from the sequencing reads"  | tee -a 0_workflow_progress.txt
    # Create kmc temporary directory
    mkdir kmc_tmp
    # Create kmc list of input read files
    ls -1 ${r1} ${r2} > kmc_input_reads.txt
    # Run software
    kmc \
        -t$(nproc --ignore=1) \
        -k21 \
        -m32 \
        -ci1 \
        -cs10000 \
        @kmc_input_reads.txt \
        kmc_count \
        kmc_tmp
    # Generate histogram for GenomeScope
    echo "Generating k-mers histogram from the sequencing reads"
    kmc_tools transform kmc_count histogram kmc_histogram.tsv -cx10000
    # Deactivate Conda environment
    conda deactivate

    # Estimate genome size using GenomeScope
    # Activate Conda environment
    conda activate genomescope
    # Run software
    echo "Estimating genome size"
    genomescope2 \
    -k 21 \
    -i kmc_histogram.tsv \
    -o "${genomesizedir}/genomescope"
    # Deactivate Conda environment
    conda deactivate

    # Verify if summary.txt exists
    summary_file="${genomesizedir}/genomescope/summary.txt"
    if [[ ! -f "$summary_file" ]]; then
        echo "✗  ERROR: GenomeScope failed for sample ${sample} (no summary.txt). Aborting." | tee -a 0_workflow_progress.txt
        echo -e "${sample}\tNA" >> 3_genomesize.tsv
        # Clean up temp files before exiting
        mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
        rm -r kmc_tmp
        exit 1
    fi

    # Extract genome size estimation
    genomesize_bp=$(grep "Genome Haploid Length" "$summary_file" | awk '{print $(NF-1)}' | tr -d ',')
    # Validate the value
    if [[ "$genomesize_bp" =~ ^[0-9]+$ ]]; then
        genomesize_mb=$(echo "scale=2; $genomesize_bp / 1000000" | bc)
        echo "Estimated genome size of sample ${sample}: ${genomesize_mb} Mb"
        # 
        [ -f 3_genomesize.tsv ] && sed -i "/^${sample}\t/d" 3_genomesize.tsv
        echo -e "${sample}\t${genomesize_bp}" >> 3_genomesize.tsv
    else
        # Unreliable GenomeScope model fit : skip this sample and keep genomescope directory for manual review.
        echo "✗  WARNING: Invalid genome size estimate (${genomesize_bp}) for sample ${sample}. Skipping sample (kept for manual review in ${genomesizedir}/genomescope)." | tee -a 0_workflow_progress.txt
        [ -f 3_genomesize.tsv ] && sed -i "/^${sample}\t/d" 3_genomesize.tsv
        echo -e "${sample}\tNA" >> 3_genomesize.tsv
        rm -r kmc_tmp kmc_count* kmc_histogram.tsv kmc_input_reads.txt
        # Increase sample count
        i=$((i + 1))
        continue
    fi

    # Move kmc temporary files to the genome size directory
    mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
    # Delete the directory kmc_tmp
    rm -r kmc_tmp

    # Declare the desired coverage
    coverage=100
    # Activate Conda environment
    conda activate rasusa
    output_r1="${r1filename/_trimmed_1.fq.gz/_trimmed_ds_1.fq.gz}"
    output_r2="${r1filename/_trimmed_1.fq.gz/_trimmed_ds_2.fq.gz}"
    echo "Downsampling the sequencing reads"  | tee -a 0_workflow_progress.txt
    # Run software
    if ! rasusa reads \
        --coverage ${coverage} \
        --genome-size ${genomesize_mb}mb \
        -s 100 \
        -o "3_fastp_downsampling/${output_r1}" \
        -o "3_fastp_downsampling/${output_r2}" \
        "${r1}" "${r2}"; then
        echo "✗  ERROR: Rasusa failed for sample ${sample}. Aborting without deleting original trimmed reads." | tee -a 0_workflow_progress.txt
        exit 1
    fi
    # Deactivate Conda environment
    conda deactivate

    # Delete the original trimmed reads after confirming rasusa actually produced valid output
    if [ -f "3_fastp_downsampling/${output_r1}" ] && [ -f "3_fastp_downsampling/${output_r2}" ] \
        && gzip -t "3_fastp_downsampling/${output_r1}" 2>/dev/null && gzip -t "3_fastp_downsampling/${output_r2}" 2>/dev/null; then
        rm "${r1}" "${r2}"
    else
        echo "✗  ERROR: Rasusa did not produce valid output for sample ${sample}. Aborting without deleting original trimmed reads." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Stop counting the running time
    loop_elapsed_time=$((SECONDS - $loop_start_time))
    # Calculate the running time
    loop_hours=$((loop_elapsed_time / 3600))
    loop_minutes=$(((loop_elapsed_time % 3600) / 60))
    loop_seconds=$((loop_elapsed_time % 60))
    loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
    # Show the running time
    echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time} " | tee -a 0_workflow_progress.txt

    # Increate sample count
    i=$((i + 1))

done

# Generate checksum file
(cd 3_fastp_downsampling && for file in *.gz; do
    [ -f "${file}.md5" ] && md5sum -c "${file}.md5" >/dev/null 2>&1 && continue
    echo "Processing checksum of file: ${file}"
    md5sum ${file} > ${file}.md5
done) | tee -a 0_workflow_progress.txt
# Check file integrity
echo "${workflow_step}: Checking file integrity of reads" | tee -a 0_workflow_progress.txt
if ! (cd 3_fastp_downsampling && md5sum -c *.md5) | tee -a 0_workflow_progress.txt; then
    echo "✗  ERROR: ${workflow_step}: read checksum verification failed" | tee -a 0_workflow_progress.txt
    exit 1
fi

# Compress the output directory
compressed_file="3_fastp_downsampling.tar.gz"
itens_to_compress=(3_fastp_downsampling/*_genomesize 3_genomesize.tsv)
if [ ${#itens_to_compress[@]} -eq 0 ]; then
    echo "✗  ERROR: ${workflow_step}: No genomesize output found to compress (3_fastp_downsampling/*_genomesize)." | tee -a 0_workflow_progress.txt
    exit 1
fi
if [ "$work_done" = false ] && [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] \
    && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (${compressed_file} verified). Skipping compression." | tee -a 0_workflow_progress.txt
else
    # Compress the output directory
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    # Generate checksum files for the reads
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
# 4) Trimmed reads quality assessment
############################################################

############################################################
## 4.1) FastQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) FastQC"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="4_fastqc.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (4_fastqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ -d "4_fastqc" ]; then
        echo "${workflow_step}: Found incomplete 4_fastqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "4_fastqc"
    fi
    rm -f "4_fastqc.tar.gz" "4_fastqc.tar.gz.md5"

    # Guard against an empty glob before calling fastqc with no input.
    fq_files=(3_fastp_downsampling/*.fq.gz)
    if [ ${#fq_files[@]} -eq 0 ]; then
        echo "✗  ERROR: No .fq.gz files found in 3_fastp_downsampling/ for FastQC." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create an output directory
    mkdir -p 4_fastqc

    # Activate Conda environment
    conda activate fastqc
    # Run main software
    fastqc -t $(nproc --ignore=1) 3_fastp_downsampling/*.gz -o 4_fastqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    itens_to_compress=(4_fastqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 4.2) FastQC -> MultiQC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="4) MultiQC"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="4_fastqc_multiqc.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (4_fastqc_multiqc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
    rm -rf "4_fastqc" 2>/dev/null
else
    if [ -d "4_fastqc_multiqc" ]; then
        echo "${workflow_step}: Found incomplete 4_fastqc_multiqc directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "4_fastqc_multiqc"
    fi
    # Delete incomplete files
    rm -f "4_fastqc_multiqc.tar.gz" "4_fastqc_multiqc.tar.gz.md5"

    # Guard against an empty glob before calling multiqc with no input.
    fastqc_zips=(4_fastqc/*_fastqc.zip)
    if [ ${#fastqc_zips[@]} -eq 0 ]; then
        echo "✗  ERROR: No *_fastqc.zip files found in 4_fastqc/ for MultiQC." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Activate Conda environment
    conda activate multiqc
    # Run main software
    multiqc "${fastqc_zips[@]}" -o 4_fastqc_multiqc
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    itens_to_compress=(4_fastqc_multiqc)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Generate checksum file of compressed directory file
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    # Check file integrity
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete output directory after the compressed file is verified valid
    rm -r 4_fastqc 4_fastqc_multiqc
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
# 5) De novo assembly
############################################################

############################################################
## 5) Unicycler

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="5) Unicycler"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Track whether any sample was actually processed
work_done=false

# Inform the output file
compressed_file="5_unicycler.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (5_unicycler.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else

    # Count and verify sample files
    i=1
    r1_files=(3_fastp_downsampling/*_trimmed_ds_1.fq.gz)
    sample_count=${#r1_files[@]}
    if [ "$sample_count" -eq 0 ]; then
        echo "✗  ${workflow_step}: No input found matching 3_fastp_downsampling/*_trimmed_ds_1.fq.gz. Check that the previous step completed successfully." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create output directory
    mkdir -p 5_unicycler

    # Activate Conda environment
    conda activate unicycler

    # Loop through a list of files
    for r1 in "${r1_files[@]}"; do
        r2=${r1/_trimmed_ds_1.fq.gz/_trimmed_ds_2.fq.gz}
        # Extract file name
        filename=${r1##*/}
        # Extract sample name
        sample=${filename%_trimmed_ds_1.fq.gz}
        out_dir="5_unicycler/${sample}_unicycler"

        # Verify if the output files exists
        if [ -s "${out_dir}/assembly.fasta" ]; then
            echo "${workflow_step} output already exists and is non-empty for sample: $sample. Skipping sample."
            # Increate sample count
            i=$((i + 1))
            continue
        elif [ -d "${out_dir}" ]; then
            echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
            rm -rf "${out_dir}"
        fi

        # Verify in the output file is missing
        if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
            echo "✗  ERROR: The input files of sample ${sample} are empty or missing. Aborting." | tee -a 0_workflow_progress.txt
            exit 1
        fi

        # The sample will actually be (re)processed  
        work_done=true

        # Inform current sample
        echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Create output directory
        mkdir -p "${out_dir}"

        # Run software
        unicycler \
            -t $(nproc --ignore=1) \
            --spades_options "--cov-cutoff auto" \
            --min_fasta_length 200 \
            -1 "${r1}" \
            -2 "${r2}" \
            -o "${out_dir}"

        # Verify if Unicycler actually produced a valid assembly before counting the sample as done
        if [ ! -s "${out_dir}/assembly.fasta" ]; then
            echo "✗  ERROR: Unicycler did not produce a valid assembly.fasta for sample ${sample}. Aborting." | tee -a 0_workflow_progress.txt
            exit 1
        fi

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time} " | tee -a 0_workflow_progress.txt

        # Increate sample count
        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Compress the output directory
    itens_to_compress=(5_unicycler)
    # Skip recompression if nothing changed and a valid archive already exists
    if [ "$work_done" = false ] && [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] \
        && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
        echo "✔  ${workflow_step}: nothing changed, ${compressed_file} already valid. Skipping compression." | tee -a 0_workflow_progress.txt
    else
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: Failed to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: Integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
## 6) Organization of de novo assembly files
############################################################

############################################################
## 6.1) Organizing assemblies

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="6) Organizing assemblies"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="6_assemblies.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (${compressed_file} verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Check if Unicycler directory or compressed file are present
    if [ ! -d "5_unicycler" ]; then
        echo "✗  ERROR: ${workflow_step}: '5_unicycler' directory not found and no valid ${compressed_file} exists. Check that step 5 completed successfully." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Check if any Unicycler sample directory is present
    unicycler_dirs=(5_unicycler/*/)
    if [ ${#unicycler_dirs[@]} -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No sample directories found in 5_unicycler/." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create output directory
    mkdir -p 6_assemblies

    for dir in "${unicycler_dirs[@]}"; do
        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_unicycler*}
        # Check if the assembly file exists and is non-empty before copying
        if [ ! -s "${dir}assembly.fasta" ]; then
            echo "✗  ERROR: ${workflow_step}: Missing or empty assembly.fasta for sample ${sample} (${dir}assembly.fasta)." | tee -a 0_workflow_progress.txt
            exit 1
        fi
        cp "${dir}assembly.fasta" "6_assemblies/${sample}.fasta"
    done

    # Compress the output directory
    itens_to_compress=(6_assemblies)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: Failed to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: Integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete output directory after the compressed file is verified valid
    rm -r 5_unicycler
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


############################################################
## 7) Assembly quality assessment
############################################################

############################################################
## 7.1) CheckM2

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) CheckM2"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="7_checkm2.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (7_checkm2.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Check for incomplete results
    output_file="7_checkm2.tsv"
    if [ -f "$output_file" ] || [ -d "7_checkm2" ]; then
        echo "${workflow_step}: found partial output from a previous interrupted run. Removing it and reprocessing." | tee -a 0_workflow_progress.txt
        rm -f "$output_file"
        rm -rf "7_checkm2"
    fi

    # Guard against empty input
    fasta_files=(6_assemblies/*.fasta)
    if [ ${#fasta_files[@]} -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No .fasta files found in 6_assemblies/." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create output directory
    mkdir -p 7_checkm2

    # Activate conda environment
    conda activate checkm2
    checkm2 predict \
        --threads $(nproc --ignore=1) \
        -x fasta \
        --input "6_assemblies" \
        --output-directory "7_checkm2"
    # Deactivate Conda environment
    conda deactivate

    # Verify if Checkm2 actually produced its report before trusting it.
    if [ ! -s "7_checkm2/quality_report.tsv" ]; then
        echo "✗  ERROR: ${workflow_step}: CheckM2 did not produce 7_checkm2/quality_report.tsv." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Copy and rename the output file
    cp "7_checkm2/quality_report.tsv" "7_checkm2.tsv"

    # Create list of filtered genomes
    echo "${workflow_step}: Filtering genomes from 7_checkm2.tsv (Completeness >= 70%, Contamination <= 5%, N50 >= 5000 bp)" | tee -a 0_workflow_progress.txt
    awk -F'\t' 'NR>1 { if ($2 >= 70 && $3 <= 5 && $7 >= 5000) print $1 }' 7_checkm2.tsv > 7_checkm2_passed.txt

    # Compress the output directory
    itens_to_compress=(7_checkm2 7_checkm2.tsv 7_checkm2_passed.txt)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete output directory after the compressed file is verified valid
    rm -r 7_checkm2
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 7.2) GUNC

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) GUNC"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="7_gunc.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (7_gunc.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ ! -f "7_checkm2_passed.txt" ]; then
        echo "${workflow_step}: 7_checkm2_passed.txt not found — attempting to extract it from 7_checkm2.tar.gz" | tee -a 0_workflow_progress.txt
        if [ -f "7_checkm2.tar.gz" ]; then
            tar -xzf 7_checkm2.tar.gz 7_checkm2_passed.txt
        fi
        if [ ! -f "7_checkm2_passed.txt" ]; then
            echo "✗ ERROR: ${workflow_step}: 7_checkm2_passed.txt could not be found or recovered. Run CheckM2 first." | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi

    # Check for incomplete results
    output_file="7_gunc.tsv"
    if [ -f "$output_file" ] || [ -d 7_gunc ] || [ -d 7_gunc_temp ]; then
        echo "${workflow_step}: found partial output from a previous interrupted run. Removing it and reprocessing." | tee -a 0_workflow_progress.txt
        rm -f "$output_file"
        rm -rf 7_gunc 7_gunc_temp
    fi

    # Build the sample list from genomes that passed CheckM2
    # Create input directory
    filtered_input_dir="7_gunc/gunc_input"
    mkdir -p "$filtered_input_dir"
    for file in 6_assemblies/*.fasta; do
        [ -f "$file" ] || continue
        filename=$(basename "$file")
        samplename="${filename%.fasta}"
        if grep -qxF "$samplename" 7_checkm2_passed.txt; then
            ln -sf "$(readlink -f "$file")" "${filtered_input_dir}/${filename}"
        fi
    done

    # Check if there are input files to analyze
    if [ -z "$(ls -A "$filtered_input_dir" 2>/dev/null)" ]; then
        echo "✗  ERROR: ${workflow_step}: No genomes passed the CheckM2 quality filter." | tee -a 0_workflow_progress.txt
        rm -rf 7_gunc 7_gunc_temp
        exit 1
    fi

    # Create output directory
    mkdir -p 7_gunc 7_gunc_temp

    # Activate conda environment
    conda activate gunc
    gunc run \
        --threads $(nproc --ignore=1) \
        --contig_taxonomy_output \
        --file_suffix .fasta \
        --input_dir "$filtered_input_dir" \
        --temp_dir "7_gunc_temp" \
        --out_dir "7_gunc"
    # Deactivate Conda environment
    conda deactivate

    # Verify if gunc actually produced output before trusting the glob/cp.
    gunc_reports=(7_gunc/*maxCSS_level.tsv)
    if [ ${#gunc_reports[@]} -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: GUNC did not produce a *maxCSS_level.tsv report." | tee -a 0_workflow_progress.txt
        exit 1
    fi
    # Copy and rename the output file
    cp "${gunc_reports[0]}" "7_gunc.tsv"

    # Delete temporary data
    rm -rf "$filtered_input_dir" 7_gunc_temp

    # Compress the output directory
    itens_to_compress=(7_gunc 7_gunc.tsv)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete the output directory
    rm -r 7_gunc
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 7.3) QUAST

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) QUAST"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="7_quast.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (7_quast.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ ! -f "7_checkm2_passed.txt" ]; then
        echo "${workflow_step}: 7_checkm2_passed.txt not found — attempting to extract it from 7_checkm2.tar.gz" | tee -a 0_workflow_progress.txt
        if [ -f "7_checkm2.tar.gz" ]; then
            tar -xzf 7_checkm2.tar.gz 7_checkm2_passed.txt
        fi
        if [ ! -f "7_checkm2_passed.txt" ]; then
            echo "✗  ERROR: ${workflow_step}: 7_checkm2_passed.txt could not be found or recovered. Run CheckM2 first." | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi

    # Build the sample list from genomes that passed CheckM2
    filtered_genomes=()
    for file in 6_assemblies/*.fasta; do
        [ -f "$file" ] || continue
        filename=$(basename "$file")
        sample="${filename%.fasta}"
        if grep -qxF "$sample" 7_checkm2_passed.txt; then
           filtered_genomes+=("$file")
        fi
    done

    # Check if there are input files to analyze
    if [ ${#filtered_genomes[@]} -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No genomes passed the CheckM2 quality filter." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create output directory
    mkdir -p 7_quast

    # Activate conda environment
    conda activate quast
    quast.py -t $(nproc --ignore=1) -m 0 -o \
        7_quast \
        "${filtered_genomes[@]}"
    # Deactivate Conda environment
    conda deactivate

    # Verify quast actually produced its report before trusting it.
    if [ ! -s "7_quast/transposed_report.tsv" ]; then
        echo "✗  ERROR: ${workflow_step}: QUAST did not produce 7_quast/transposed_report.tsv." | tee -a 0_workflow_progress.txt
        exit 1
    fi
    cp "7_quast/transposed_report.tsv" "7_quast.tsv"

    # Compress the output directory
    itens_to_compress=(7_quast 7_quast.tsv)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Only delete raw output after the archive is verified valid.
    rm -r 7_quast
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 7.4) Pybarrnap

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) Pybarrnap"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform the output file
compressed_file="7_pybarrnap.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (7_pybarrnap.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ ! -f "7_checkm2_passed.txt" ]; then
        echo "${workflow_step}: 7_checkm2_passed.txt not found — attempting to extract it from 7_checkm2.tar.gz" | tee -a 0_workflow_progress.txt
        if [ -f "7_checkm2.tar.gz" ]; then
            tar -xzf 7_checkm2.tar.gz 7_checkm2_passed.txt
        fi
        if [ ! -f "7_checkm2_passed.txt" ]; then
            echo "✗  ERROR: ${workflow_step}: 7_checkm2_passed.txt could not be found or recovered. Run CheckM2 first." | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi

    # Build the sample list from genomes that passed CheckM2
    filtered_genomes=()
    for file in 6_assemblies/*.fasta; do
        [ -f "$file" ] || continue
        filename=$(basename "$file")
        sample="${filename%.fasta}"
        if grep -qxF "$sample" 7_checkm2_passed.txt; then
           filtered_genomes+=("$file")
        fi
    done

    # Count and verify sample files
    i=1
    sample_count=${#filtered_genomes[@]}
    # Check if there are input files to analyze
    if [ ${#filtered_genomes[@]} -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No genomes passed the CheckM2 quality filter." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create output directory
    mkdir -p 7_pybarrnap

    # Activate conda environment
    conda activate pybarrnap
    for file in ${filtered_genomes[@]}; do
        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%.fasta}

        # Inform current sample
        echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Skip sample if bac output exists and is non-empty AND arc output at least exists
        # (arc output is legitimately empty for bacterial genomes with no archaeal rRNA hits,
        # so it must be checked with -f, not -s, or the sample is reprocessed forever)
        bac_output="7_pybarrnap/${sample}_bac_pybarrnap.fasta"
        arc_output="7_pybarrnap/${sample}_arc_pybarrnap.fasta"
        if [ -s "$bac_output" ] && [ -f "$arc_output" ]; then
            echo "${workflow_step} output files already exist and are valid for sample: $sample. Skipping sample."
            i=$((i + 1))
            continue
        fi

        # Run barrnap for archea
        cat "$file" | pybarrnap \
            --threads $(nproc --ignore=1) \
            --quiet \
            --kingdom arc \
            --outseq "$arc_output" \
            > "7_pybarrnap/${sample}_arc_pybarrnap.gff" \
            2> "7_pybarrnap/${sample}_arc_pybarrnap.log"

        # Run barrnap for bacteria
        cat "$file" | pybarrnap \
            --threads $(nproc --ignore=1) \
            --quiet \
            --kingdom bac \
            --outseq "$bac_output" \
            > "7_pybarrnap/${sample}_bac_pybarrnap.gff" \
            2> "7_pybarrnap/${sample}_bac_pybarrnap.log"

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time} " | tee -a 0_workflow_progress.txt

        # Increate sample count
        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Merge gff files
    echo "${workflow_step}: Rebuilding aggregated GFF files" | tee -a 0_workflow_progress.txt
    > 7_pybarrnap_bac_pybarrnap.gff
    > 7_pybarrnap_arc_pybarrnap.gff
    for file in "${filtered_genomes[@]}"; do
        filename=${file##*/}
        sample=${filename%.fasta}

        echo "${sample}" >> 7_pybarrnap_bac_pybarrnap.gff
        grep -v '^#' "7_pybarrnap/${sample}_bac_pybarrnap.gff" >> 7_pybarrnap_bac_pybarrnap.gff

        echo "${sample}" >> 7_pybarrnap_arc_pybarrnap.gff
        grep -v '^#' "7_pybarrnap/${sample}_arc_pybarrnap.gff" >> 7_pybarrnap_arc_pybarrnap.gff
    done

    # Compress the output directory
    itens_to_compress=(7_pybarrnap)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Only delete raw output after the archive is verified valid.
    rm -r 7_pybarrnap
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 7.6) Vertical coverage

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="7) Vertical coverage"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Ensure 7_checkm2_passed.txt is present
if [ ! -f "7_checkm2_passed.txt" ]; then
    echo "${workflow_step}: 7_checkm2_passed.txt not found — attempting to extract it from 7_checkm2.tar.gz" | tee -a 0_workflow_progress.txt
    if [ -f "7_checkm2.tar.gz" ]; then
        tar -xzf 7_checkm2.tar.gz 7_checkm2_passed.txt
    fi
    if [ ! -f "7_checkm2_passed.txt" ]; then
        echo "✗  ERROR: ${workflow_step}: 7_checkm2_passed.txt could not be found or recovered. Run CheckM2 first." | tee -a 0_workflow_progress.txt
        exit 1
    fi
fi

# Build the sample list from genomes that passed CheckM2
i=1
files=()
for file in 6_assemblies/*.fasta; do
    [ -f "$file" ] || continue
    filename=$(basename "$file")
    sample="${filename%.fasta}"
    if grep -qxF "$sample" 7_checkm2_passed.txt; then
        files+=("$file")
    fi
done

# Check if there are input files to analyze
sample_count=${#files[@]}
if [ "$sample_count" -eq 0 ]; then
    echo "✗  ERROR: ${workflow_step}: No genomes passed the CheckM2 quality filter." | tee -a 0_workflow_progress.txt
    exit 1
fi

# Start the analyis
if [ -f "7_coverage.tsv" ] && [ -f "7_coverage.tsv.md5" ] && md5sum -c "7_coverage.tsv.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (7_coverage.tsv verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create output file
    echo -e Sample"\t"Coverage > 7_coverage.tsv

    for file in "${files[@]}"; do
        # Extract file name
        filename=${file##*/}
        # Extract sample name
        sample=${filename%.fasta}

        # Inform current sample
        echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        r1=(3_fastp_downsampling/${sample}_trimmed_ds_1.fq.gz)
        r2=(3_fastp_downsampling/${sample}_trimmed_ds_2.fq.gz)
        echo -e Assembly file: ${file}
        echo -e R1 file: $r1
        echo -e R2 file: $r2

        # Check the input files
        if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
            echo "✗  ERROR: The input files of sample ${sample} are empty. Aborting." | tee -a 0_workflow_progress.txt
            exit 1
        fi

        # Calculate coverage
        bases_in_assembly=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$file")
        echo -e Bases in assembly: $bases_in_assembly
        bases_in_reads=$(zcat "$r1" "$r2" | awk 'NR%4==2 {sum += length($0)} END {print sum}')
        echo -e Bases in reads: $bases_in_reads
        if [ "$bases_in_assembly" -gt 0 ]; then
            coverage=$(echo "scale=2; $bases_in_reads / $bases_in_assembly" | bc)
            echo -e Coverage: "${coverage}\n"
            echo -e "${sample}\t${coverage}" >> 7_coverage.tsv
        else
            echo -e "${sample}\t0" >> 7_coverage.tsv
        fi

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time} " | tee -a 0_workflow_progress.txt

        # Increate sample count
        i=$((i + 1))

    done

    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    output_file="7_coverage.tsv"
    md5sum "${output_file}" > "${output_file}".md5
    if ! md5sum -c "${output_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${output_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 8) Taxonomic assignment
############################################################

############################################################
## 8.1) GTDB-Tk

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="8) GTDB-Tk"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Inform output file
compressed_file="8_gtdbtk.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (8_gtdbtk.tar.gz verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    if [ -d "8_gtdbtk" ]; then
        echo "${workflow_step}: Found incomplete 8_gtdbtk directory from a previous interrupted run. Removing it to start fresh." | tee -a 0_workflow_progress.txt
        rm -rf "8_gtdbtk"
    fi
    # Delete incomplete results
    rm -f "8_gtdbtk.tar.gz" "8_gtdbtk.tar.gz.md5"

    # Create output directory
    mkdir -p 8_gtdbtk

    # Check the input list
    if [ ! -f "7_checkm2_passed.txt" ]; then
        echo "${workflow_step}: 7_checkm2_passed.txt not found — attempting to extract it from 7_checkm2.tar.gz" | tee -a 0_workflow_progress.txt
        if [ -f "7_checkm2.tar.gz" ]; then
            tar -xzf 7_checkm2.tar.gz 7_checkm2_passed.txt
        fi
        if [ ! -f "7_checkm2_passed.txt" ]; then
            echo "✗  ERROR: ${workflow_step}: 7_checkm2_passed.txt could not be found or recovered. Run CheckM2 first." | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi

    # Inform the number of samples 
    n_passed=$(wc -l < 7_checkm2_passed.txt)
    echo "${workflow_step}: ${n_passed} genome(s) passed quality control filtering." | tee -a 0_workflow_progress.txt

    # Check the input directory
    sampledir="6_assemblies"
    if [ ! -d "$sampledir" ]; then
        echo "✗  ERROR: ${workflow_step}: No input directory found matching '6_assemblies'. Check that the previous steps completed successfully." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Create batch file
    echo "${workflow_step}: Building batchfile from 6_assemblies"
    batchfile="8_gtdbtk_batchfile.tsv"
    > "$batchfile"

    # Build the sample list from genomes that passed CheckM2
    for file in 6_assemblies/*.fasta; do
        [ -f "$file" ] || continue
        filename=$(basename "$file")
        samplename="${filename%.fasta}"
        if grep -qxF "$samplename" 7_checkm2_passed.txt; then
            printf '%s\t%s\n' "$(readlink -f "$file")" "$samplename" >> "$batchfile"
        fi
    done

    n_genomes=$(wc -l < "$batchfile")
    echo "${workflow_step} batchfile has ${n_genomes} quality-filtered genome(s) across $(ls -d 6_assemblies/*.fasta | wc -l) sample(s)"

    # Guard against an empty batchfile
    if [ "$n_genomes" -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No genomes passed the CheckM2 quality filter — batchfile is empty." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Inform the begining of the pipeline
    echo "${workflow_step} started the workflow @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt

    # Create scratch directory
    scratch_dir="/tmp/gtdbtk_scratch_${USER}_$$"
    mkdir -p "$scratch_dir"

    # Activate conda environment
    conda activate gtdbtk
    gtdbtk classify_wf \
        --cpus $(nproc --ignore=1) \
        --pplacer_cpus 12 \
        --scratch_dir "$scratch_dir" \
        --batchfile "$batchfile" \
        --out_dir "8_gtdbtk"
    # Deactivate Conda environment
    conda deactivate

    # Delete the temporary data
    rm -rf "$scratch_dir"

    # Check the output file
    if [ ! -f "8_gtdbtk/classify/gtdbtk.bac120.summary.tsv" ] && [ ! -f "8_gtdbtk/classify/gtdbtk.ar53.summary.tsv" ]; then
        echo "✗  ERROR: ${workflow_step}: GTDB-Tk did not produce any summary file (classify/gtdbtk.bac120.summary.tsv or .ar53.summary.tsv). Aborting without deleting 8_gtdbtk." | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Rename the output files
    for map in "classify/gtdbtk.bac120.summary.tsv:8_gtdbtk_bacteria.tsv" \
               "classify/gtdbtk.ar53.summary.tsv:8_gtdbtk_archaea.tsv" \
               "gtdbtk.log:8_gtdbtk.log" ; do
        src="8_gtdbtk/${map%%:*}"
        dst="8_gtdbtk/${map#*:}"
        [ -f "$src" ] && cp "$src" "$dst"
    done

    # Merge the output files
    first_file=""
    for f in "8_gtdbtk/8_gtdbtk_bacteria.tsv" "8_gtdbtk/8_gtdbtk_archaea.tsv"; do
        [ -f "$f" ] && first_file="$f" && break
    done
    head -n 1 "$first_file" > "8_gtdbtk/8_gtdbtk_archaea_bacteria.tsv"
    for f in "8_gtdbtk/8_gtdbtk_bacteria.tsv" "8_gtdbtk/8_gtdbtk_archaea.tsv"; do
        [ -f "$f" ] || continue
        tail -n +2 "$f" >> "8_gtdbtk/8_gtdbtk_archaea_bacteria.tsv"
    done

    # Verify the if combined file was actually built before copying it out
    if [ ! -s "8_gtdbtk/8_gtdbtk_archaea_bacteria.tsv" ]; then
        echo "✗  ERROR: ${workflow_step}: Failed to build 8_gtdbtk_archaea_bacteria.tsv. Aborting without deleting 8_gtdbtk." | tee -a 0_workflow_progress.txt
        exit 1
    fi
    cp 8_gtdbtk/8_gtdbtk_archaea_bacteria.tsv 8_gtdbtk.tsv

    # Compress the output directory
    itens_to_compress=(8_gtdbtk 8_gtdbtk.tsv)
    echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
    if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
        echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum "${compressed_file}" > "${compressed_file}".md5
    echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
    if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Delete the output directory
    rm -r 8_gtdbtk
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

############################################################
## 8.2) TYGS

# https://tygs.dsmz.de/user_requests/new
# Send the files from 6_assemblies


############################################################
## 9) Plasmid identification
############################################################

############################################################
## 9.1) MOB-suite

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="9) MOB-suite"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Function to merge per-sample TSV files into one.
merge_tsv_files() {
    local pattern="$1"
    local search_dir="$2"
    local output_file="$3"

    [ -f "$output_file" ] && rm -f "$output_file"

    local header_printed=0
    local file
    while IFS= read -r file; do
        if [ "$header_printed" -eq 0 ]; then
            cat "$file" >> "$output_file"
            header_printed=1
        else
            tail -n +2 "$file" >> "$output_file"
        fi
    done < <(find "$search_dir" -maxdepth 2 -type f -name "$pattern")
}

# Inform the output file
compressed_file="9_mobsuite.tar.gz"

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (${compressed_file} verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Create output directory
    mkdir -p 9_mobsuite

    # Calculate sample size
    i=1
    files=(6_assemblies/*.fasta)
    sample_count=${#files[@]}
    if [ "$sample_count" -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No input found matching 6_assemblies/*.fasta" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Track whether any sample was actually processed
    work_done=false

    # Activate conda environment
    conda activate mob_suite
    for file in "${files[@]}"; do
        # Extract file name
        filename=${file#*/}
        # Extract sample name
        sample=${filename%.fasta}
        outdir="9_mobsuite/${sample}_mobsuite"

        if [ -s "${outdir}/contig_report.txt" ]; then
            echo "${workflow_step} output already exists and is valid for sample: $sample. Skipping genome."
            # Increate sample count
            i=$((i + 1))
            continue
        elif [ -d "$outdir" ]; then
            echo "${workflow_step} found incomplete output for sample: $sample. Removing partial directory and reprocessing."
            rm -rf "$outdir"
        fi

        # The sample will actually be (re)processed
        work_done=true

        # Inform current sample
        echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Run the software
        mob_recon -n $(nproc --ignore=1) \
            --infile "${file}" \
            --outdir "$outdir"

        # Check the output file
        if [ ! -s "${outdir}/contig_report.txt" ]; then
            echo "✗  ERROR: ${workflow_step}: mob_recon did not produce ${outdir}/contig_report.txt for sample ${sample}. Aborting." | tee -a 0_workflow_progress.txt
            exit 1
        fi

        # Rename chromosome and plasmid files 
        (
            cd "$outdir" || exit
            [ -f chromosome.fasta ] && mv chromosome.fasta "${sample}"_chromosome.fasta
            for p in plasmid*; do
                [ -e "$p" ] || continue
                mv "$p" "${sample}_$p"
            done
        )

        # Stop counting the running time
        loop_elapsed_time=$((SECONDS - $loop_start_time))
        # Calculate the running time
        loop_hours=$((loop_elapsed_time / 3600))
        loop_minutes=$(((loop_elapsed_time % 3600) / 60))
        loop_seconds=$((loop_elapsed_time % 60))
        loop_running_time=$(printf "%02d:%02d:%02d" "$loop_hours" "$loop_minutes" "$loop_seconds")
        # Show the running time
        echo "✔  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${loop_running_time} " | tee -a 0_workflow_progress.txt
        
        # Increate sample count
        i=$((i + 1))

    done
    # Deactivate Conda environment
    conda deactivate

    # Merge per-sample MOB-suite outputs into combined files
    merge_tsv_files "contig_report.txt" "9_mobsuite" "9_mobsuite/contig_report_all.tsv"
    merge_tsv_files "mobtyper_results.txt" "9_mobsuite" "9_mobsuite/mobtyper_results_all.tsv"
    merge_tsv_files "mge.report.txt" "9_mobsuite" "9_mobsuite/mge.report_all.tsv"

    # Copy merged mobtyper result file to main directory
    if [ -f "9_mobsuite/mobtyper_results_all.tsv" ]; then
        cp 9_mobsuite/mobtyper_results_all.tsv 9_mobsuite_mobtyper_results_all.tsv
    fi

    # Compress the output directory
    itens_to_compress=(9_mobsuite)
    [ -f "9_mobsuite_mobtyper_results_all.tsv" ] && itens_to_compress+=("9_mobsuite_mobtyper_results_all.tsv")
    if [ "$work_done" = false ] && [ -f "$compressed_file" ] && [ -f "${compressed_file}.md5" ] \
        && md5sum -c "${compressed_file}.md5" >/dev/null 2>&1; then
        echo "✔  ${workflow_step} already completed successfully (${compressed_file} verified). Skipping compression." | tee -a 0_workflow_progress.txt
    else
        echo "${workflow_step}: Compressing output directory" | tee -a 0_workflow_progress.txt
        if ! tar -c --use-compress-program=pigz -f "${compressed_file}" "${itens_to_compress[@]}"; then
            echo "✗  ERROR: ${compressed_file}: FAILED to create archive" | tee -a 0_workflow_progress.txt
            exit 1
        fi
        md5sum "${compressed_file}" > "${compressed_file}".md5
        echo "${workflow_step}: Checking file integrity" | tee -a 0_workflow_progress.txt
        if ! md5sum -c "${compressed_file}".md5 | tee -a 0_workflow_progress.txt; then
            echo "✗  ERROR: ${compressed_file}: integrity check failed" | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt


#########################################################################
## 10) Assignment of contigs to molecules
############################################################

############################################################
## Add molecule attribution to contigs. Required for batch genome submission.

# Avoid literal glob pattern
shopt -s nullglob

# Software name for tracking progress in 0_workflow_progress.txt
workflow_step="10) Molecule attribution"
# Update the file 0_workflow_progress.txt
echo "▶▶▶  ${workflow_step} started @ $(date +'%Y-%m-%d %H:%M:%S') ▶▶▶" | tee -a 0_workflow_progress.txt
# Start counting the running time
start_time=$SECONDS

# Skip if this step already completed successfully (final archive present and valid)
if [ -f "10_assemblies_for_analysis.zip" ] && [ -f "10_assemblies_for_analysis.zip.md5" ] \
    && [ -f "9_mobsuite.zip" ] && [ -f "9_mobsuite.zip.md5" ] \
    && md5sum -c "10_assemblies_for_analysis.zip.md5" >/dev/null 2>&1 \
    && md5sum -c "9_mobsuite.zip.md5" >/dev/null 2>&1; then
    echo "✔  ${workflow_step} already completed successfully (both archives verified). Skipping step." | tee -a 0_workflow_progress.txt
else
    # Delete incomplete output
    rm -rf 10_assemblies_for_analysis
    rm -f 10_assemblies_for_analysis.zip 10_assemblies_for_analysis.zip.md5
    rm -f 9_mobsuite.zip 9_mobsuite.zip.md5

    # Check the MOB-suite directory
    if [ ! -d "9_mobsuite" ]; then
        if [ -f "9_mobsuite.tar.gz" ]; then
            echo "${workflow_step}: 9_mobsuite/ not found — extracting from 9_mobsuite.tar.gz" | tee -a 0_workflow_progress.txt
            tar -xzf 9_mobsuite.tar.gz
        fi
        if [ ! -d "9_mobsuite" ]; then
            echo "✗  ERROR: ${workflow_step}: '9_mobsuite' directory not found and could not be recovered. Check that step 9.1 completed successfully." | tee -a 0_workflow_progress.txt
            exit 1
        fi
    fi

    # Copy and rename fasta files
    cp -r 6_assemblies 10_assemblies_for_analysis
    rename "s/.fasta$/.fsa/" 10_assemblies_for_analysis/*.fasta

    # Calculate sample size
    i=1
    files=(10_assemblies_for_analysis/*.fsa)
    sample_count=${#files[@]}
    if [ "$sample_count" -eq 0 ]; then
        echo "✗  ERROR: ${workflow_step}: No input found matching 10_assemblies_for_analysis/*.fsa" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    for assembly in 10_assemblies_for_analysis/*.fsa; do
        # Extract file name
        filename=${assembly##*/}
        # Extract sample name
        sample=${filename%%.*}

        # Check the input file
        report="9_mobsuite/${sample}_mobsuite/contig_report.txt"
        # Check the report exists before parsing it
        if [ ! -s "$report" ]; then
            echo "✗  ERROR: ${workflow_step}: Missing or empty ${report} for sample ${sample}." | tee -a 0_workflow_progress.txt
            exit 1
        fi

        # Inform current sample
        echo "▶  ${workflow_step} — ${sample} (${i}/${sample_count}) @ $(date +'%Y-%m-%d %H:%M:%S')" | tee -a 0_workflow_progress.txt
        # Start counting the running time
        loop_start_time=$SECONDS

        # Assignment of contigs to molecules
        awk '1' "$report" | while IFS=$'\t' read -r sample_id molecule_type primary_cluster_id secondary_cluster_id contig_id others; do
            if [ "$sample_id" == "sample_id" ]; then
                continue
            fi
            # Escape contig_id before using it as a sed regex pattern
            escaped_contig_id=$(printf '%s' "$contig_id" | sed 's/[.[\*^$/]/\\&/g')
            if [ "$molecule_type" == "plasmid" ]; then
                new_contig_id=$(echo "$contig_id" "[plasmid-name="p"$primary_cluster_id""]")
                echo Sample: "$sample" - Contig: "$new_contig_id"
                sed -i "s/^>${escaped_contig_id}\$/>${new_contig_id}/" "$assembly"
            fi
            if [ "$molecule_type" == "chromosome" ]; then
                new_contig_id=$(echo "$contig_id" "[chromosome=1]")
                echo Sample: "$sample" - Contig: "$new_contig_id"
                sed -i "s/^>${escaped_contig_id}\$/>${new_contig_id}/" "$assembly"
            fi
        done

        # Increate sample count
        i=$((i + 1))
    done

    # Check the output
    if ! zip -r 10_assemblies_for_analysis.zip 10_assemblies_for_analysis; then
        echo "✗  ERROR: ${workflow_step}: FAILED to create 10_assemblies_for_analysis.zip" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum 10_assemblies_for_analysis.zip > 10_assemblies_for_analysis.zip.md5
    if ! md5sum -c 10_assemblies_for_analysis.zip.md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: 10_assemblies_for_analysis.zip: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    if ! zip -r 9_mobsuite.zip 9_mobsuite; then
        echo "✗  ERROR: ${workflow_step}: FAILED to create 9_mobsuite.zip" | tee -a 0_workflow_progress.txt
        exit 1
    fi
    md5sum 9_mobsuite.zip > 9_mobsuite.zip.md5
    if ! md5sum -c 9_mobsuite.zip.md5 | tee -a 0_workflow_progress.txt; then
        echo "✗  ERROR: 9_mobsuite.zip: integrity check failed" | tee -a 0_workflow_progress.txt
        exit 1
    fi

    # Only delete temporary data
    rm -r 10_assemblies_for_analysis
    rm -r 9_mobsuite
fi

# Stop counting the running time
elapsed_time=$((SECONDS - $start_time))
# Calculate the running time
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
running_time=$(printf "%02d:%02d:%02d" "$hours" "$minutes" "$seconds")
# Update the file 0_workflow_progress.txt
echo -e "■■■  ${workflow_step} finished @ $(date +'%Y-%m-%d %H:%M:%S') — Total: ${running_time} ■■■\n" | tee -a 0_workflow_progress.txt

# In case of a novel plasmid, you will have to change its temporary name given by MOB-suite to an appropriate and shorter name.

############################################################
## Genome submission to GenBank

# Submission portal
# https://submit.ncbi.nlm.nih.gov/subs/

# Create BioProject
# https://submit.ncbi.nlm.nih.gov/subs/bioproject/

# Create Biosample
# https://submit.ncbi.nlm.nih.gov/subs/biosample/
# BioSample template (Required form)
# https://submit.ncbi.nlm.nih.gov/biosample/template/

# Submit the genome
# https://submit.ncbi.nlm.nih.gov/subs/genome/
# Choose the "Batch submission", even for a single submission: “New submission” escolher a opção “Batch/multiple genomes (maximum 400 per submission)”.
# Choose the genome annotation with PGAP: "Annotate this prokaryotic genome in the NCBI Prokaryotic Annotation Pipeline (PGAP) before its release"
# Upload a table containing the genomes metadata. A template (Batch genomes: Genome Info file template) is available in https://submit.ncbi.nlm.nih.gov/templates/.
# Upload the .fsa files in the directory 10_assemblies_for_analysis

# Submit sequencing reads
# https://submit.ncbi.nlm.nih.gov/subs/sra/
# Upload a table containing the sequencing metadata. A template (SRA: Metadata spreadsheet with sample names) is available in https://submit.ncbi.nlm.nih.gov/templates/.
# Upload the .fq.gz files in directory 1_reads
