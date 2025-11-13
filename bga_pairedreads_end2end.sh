#!/bin/bash
# Bash script for bacterial genome assembly from short-read sequencing data (end-to-end worflow)
#
# Author: Marcus Vinicius Canário Viana
# Date: 13/11/2025
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
# 1. Create a **tab-separated file** named **"1_reads_accessions.tsv"**.
# 2. This file **must contain** the NCBI SRA **accession number** in the first column and the **sample name** in the second column. Other columns will be ignored.
# 3. **Do not use** special characters in the sample names.
# 4. Place the **"1_reads_accessions.tsv"** file in the working directory.
#
# **C. Execution**
#
# Place this script (**bga_pairedreads_end2end.sh**) in the working directory and execute it **using the following commands**:
# chmod +x bga_shortreads_end2end.sh
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
    echo "ERROR: The 'conda' command was not found." >&2
    echo "Ensure Conda or Miniconda is installed and configured in your PATH." >&2
    # Exit the script
    exit 1
fi


############################################################
## 1) Sequencing reads directory and files
############################################################

############################################################
## Reads from NCBI SRA

# Delete previous file of not used reads
rm -f 1_reads_single-end.tsv

# Verify the presence of the file 1_reads_accessions.tsv with a list of accessions
if [ -f 1_reads_accessions.tsv ]; then
    echo "############################################################"
    echo "# Downloading reads from NCBI SRA"
    echo "############################################################"
    echo "The file 1_reads_accessions.tsv was found. The sequencing reads will be downloaded."

    # Create output directory 
    mkdir -p 1_reads

    # SRA Tools
    # Activate Conda environment
    conda activate sra-tools

    # Loop through file lines
    awk '1' 1_reads_accessions.tsv | while IFS=$'\t' read -r accession sample others; do
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
                # Go back to main directory
                cd ..
            # Check if the reads are not paired-end
            elif [ -f "${accession}.fastq.gz" ]; then
                echo "Sample ${sample} has single-end reads."
                # Renaming file
                echo "Renaming file."
                mv "${accession}.fastq.gz" "${sample}.fq.gz"
                # Warning about the requirement of paired-reads
                echo "Warning: this script only use paired-end reads. This file will not be used." 
                # Create list of not used read files
                echo -e "${accession}\t${sample}.fq.gz" >> ../1_reads_single-end.tsv
                # Go back to main directory
                cd ..
            fi
        fi
    done
    echo "Download process complete. Deactivating the environment."
    # Deactivate Conda environment
    conda activate base
else
    echo "The file 1_reads_accessions.tsv was not found. Proceeding using local files."
fi

############################################################
## Reads stored as local files

echo "############################################################"
echo "# Checking local read files"
echo "############################################################"

# Checking wether the directory 1_reads exists
if [ ! -d 1_reads ]; then
    echo "The reads directory '1_reads' was not found" >&2
    echo "Please create it and put the files in it" >&2
    # Exit the script
    exit 1
fi

# The paired-end file names of each sample should be standardized to samplename_1.fq.gz and samplename_2.fq.gz

# Checking the presence of files in the format 1_reads/*_1.fastq.gz and 1_reads/*_2.fastq.gz
if ls 1_reads/*_1.fastq.gz > /dev/null 2>&1; then
    echo "Found files in the format *_1.fastq.gz and *_2.fastq.gz" >&2
    echo "Renaming them and their pairs to the format *_1.fq.gz and *_2.fq.gz" >&2
    # Rename the .fastq.gz files
    rename 's/\.fastq\.gz$/.fq.gz/g' 1_reads/*.fastq.gz

    # Checking the presence of .md5 files
    if ls 1_reads/*_1.fastq.gz.md5 > /dev/null 2>&1; then
        # Rename the file names inside the .md5 files
        sed -i 's/\.fastq\.gz$/.fq.gz/' 1_reads/*.md5
        # Rename the .md5 files
        rename 's/\.fastq\.gz.md5$/.fq.gz.md5/' 1_reads/*.fastq.gz.md5
    fi

else
    :
fi

# Checking the presence of files in the format 1_reads/*_R1_001.fastq.gz and 1_reads/*_2.fastq.gz
if ls 1_reads/*_R1_001.fastq.gz > /dev/null 2>&1; then
    echo "Found files in the format *_R1_001.fastq.gz and *_2.fastq.gz" >&2
    echo "Renaming them and their pairs to the format *_1.fq.gz and *_2.fq.gz" >&2
    # Rename the .fastq.gz files
    rename 's/_R1_001\.fastq\.gz/_1.fq.gz/; s/_R2_001\.fastq\.gz/_2.fq.gz/' 1_reads/*.fastq.gz

    # Checking the presence of .md5 files
    if ls 1_reads/*_R1_001.fastq.gz.md5 > /dev/null 2>&1; then
        # Rename the file names inside the .md5 files
        sed -i 's/_R1_001\.fastq\.gz/_1.fq.gz/; s/_R2_001\.fastq\.gz/_2.fq.gz/' 1_reads/*.md5
        # Rename the .md5 files
        rename 's/_R1_001\.fastq\.gz\.md5/_1.fq.gz\.md5/; s/_R2_001\.fastq\.gz\.md5/_2.fq.gz\.md5/' 1_reads/*.fastq.gz.md5
    fi

else
    :
fi

## Checking the presence, names and pairing of the input files
files_found=0
missing_pairs=0
#Loop
for r1 in 1_reads/*_1.fq.gz; do
    #Check the presence of r1 file
    if [[ "$r1" == "1_reads/*_1.fq.gz" ]]; then
        # Exit the script
        break 
    fi

    # If at least one *_1.fq.gz file was found
    files_found=1

    # Extract r2 file name
    r2="${r1/_1.fq.gz/_2.fq.gz}"

    #Check the presence of r2 file
    if [ ! -f "$r2" ]; then
        echo "Expected file "${r2}" not found" >&2
        missing_pairs=1
    fi
done

# Final reports
if [ "$files_found" -eq 0 ]; then
    echo "No read files in the format 'samplename_1.fq.gz' was found in the diretory '1_reads'" >&2
    exit 1
fi
if [ "$missing_pairs" -eq 1 ]; then
    echo "The script requires read pairs for all samples in the format *_1.fq.gz and *_2.fq.gz" >&2
    echo "Rename them or put them in the diretory 1_reads" >&2
    # Exit the script
    exit 1
fi


############################################################
## THE SCRIPT BELOW ARE LINES SELECTED FROM THE FILE bga_shortreads_script.sh
############################################################


############################################################
## 2) Raw reads quality assessment
############################################################

############################################################
## FastQC

echo "############################################################"
echo "# Running FastQC for raw reads"
echo "############################################################"

# Create an output directory
mkdir 2_fastqc
# Activate Conda environment
conda activate fastqc
# Run FastQC
fastqc -t $(nproc --ignore=1) 1_reads/*.fq.gz -o 2_fastqc
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 2_fastqc.zip 2_fastqc

############################################################
## FastQC -> MultiQC

echo "############################################################"
echo "# Running MultiQC for raw reads"
echo "############################################################"

# Activate Conda environment
conda activate multiqc
# Run MultiQC
multiqc 2_fastqc/*_fastqc.zip -o 2_fastqc_multiqc
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 2_fastqc_multiqc.zip 2_fastqc_multiqc
# Delete the output directory
rm -r 2_fastqc 2_fastqc_multiqc


############################################################
## 3) Raw reads trimming, estimation of genome size and downsampling 
############################################################

############################################################
## Fastp

echo "############################################################"
echo "# Running Fastp"
echo "############################################################"

# Create an output directory
mkdir 3_fastp
# Activate Conda environment
conda activate fastp
# Loop through a list of files
for r1 in 1_reads/*_1.fq.gz; do
    # Extract r2 path
    r2="${r1/_1.fq.gz/_2.fq.gz}"
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%%_*}

    # Verify if the input files are empty
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 3_fastp_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    # Inform the current sample being processed
    echo "Processing sample: ${sample}"
    echo "R1: ${r1}"
    echo "R2: ${r2}"
    # Run Fastp
    fastp \
        --thread $(nproc --ignore=1) \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --cut_front \
        --cut_right \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        --overrepresentation_analysis \
        --in1 "$r1" \
        --in2 "$r2" \
        --out1 "3_fastp/${sample}_trimmed_1.fq.gz" \
        --out2 "3_fastp/${sample}_trimmed_2.fq.gz" \
        --html "3_fastp/${sample}_trimmed_fastp.html" \
        --json "3_fastp/${sample}_trimmed_fastp.json"
done
# Deactivate Conda environment
conda activate base
# Compress report files
zip -r 3_fastp.zip 3_fastp/*.json 3_fastp/*.html
# Delete report files
rm 3_fastp/*.json 3_fastp/*.html

############################################################
## Estimation of genome size (KMC, GenomeScope) and downsampling  (Rasusa)

echo "############################################################"
echo "# Running estimation of genome size and downsampling"
echo "############################################################"

# Required coverage
coverage=100

# Loop through a list of files
for r1 in 3_fastp/*_trimmed_1.fq.gz; do
    # Extract r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    r1filename=${r1##*/}
    # Extract sample name
    sample=${r1filename%%_*}

    # Verify if the input files are empty
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 3_genomesize_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    # Create directory for the genome size estimation results
    genomesizedir="3_fastp/${sample}_genomesize"
    mkdir ${genomesizedir}
    
    # Count k-mers using KMC
    # Activate Conda environment
    conda activate kmc
    # Count k-mers
    echo "Counting k-mers from the sequencing reads of sample ${sample}"
    # Create kmc temporary directory
    mkdir kmc_tmp
    # Create kmc list of input read files
    ls -1 ${r1} ${r2} > kmc_input_reads.txt
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
    echo "Generating k-mers histogram from the sequencing reads of sample ${sample}"
    kmc_tools transform kmc_count histogram kmc_histogram.tsv -cx10000
    # Deactivate Conda environment
    conda activate base

    # Estimate genome size using GenomeScope
    # Activate Conda environment
    conda activate genomescope
    #Run software
    echo "Estimating genome size of sample ${sample}"
    genomescope2 \
    -k 21 \
    -i kmc_histogram.tsv \
    -o "${genomesizedir}/genomescope"
    # Deactivate Conda environment
    conda activate base

    # Verify that summary.txt exists
    summary_file="${genomesizedir}/genomescope/summary.txt"
    if [[ ! -f "$summary_file" ]]; then
        echo "Warning: GenomeScope failed for sample ${sample} (no summary.txt). Skipping sample..."
        echo -e "${sample}\tNA" >> 3_genomesize.tsv
        # Clean up temp files before skipping
        mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
        rm -r kmc_tmp
        # Skip to the next iteration
        continue
    fi

    # Extract genome size estimate
    genomesize_bp=$(grep "Genome Haploid Length" "$summary_file" | awk '{print $(NF-1)}' | tr -d ',')
    # Validate the value
    if [[ "$genomesize_bp" =~ ^[0-9]+$ ]]; then
        genomesize_mb=$(echo "scale=2; $genomesize_bp / 1000000" | bc)
        echo "Estimated genome size of sample ${sample}: ${genomesize_mb} Mb"
        echo -e "${sample}\t${genomesize_bp}" >> 3_genomesize.tsv
    else
        echo "Warning: Invalid genome size estimate for sample ${sample} (${genomesize_bp}). Skipping sample..."
        echo -e "${sample}\tNA" >> 3_genomesize.tsv
        # Move kmc temporary files to the genome size directory
        mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
        # Delete the directory kmc_tmp
        rm -r kmc_tmp
        # Skip to the next iteration
        continue
    fi

    # Move kmc temporary files to the genome size directory
    mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
    # Delete the directory kmc_tmp
    rm -r kmc_tmp

    # Downsampling reads using Rasusa
    # Activate Conda environment
    conda activate rasusa
    # Generate output files names with sufix "ds"
    output_r1="${r1filename/_1.fq.gz/_ds_1.fq.gz}"
    output_r2="${r1filename/_1.fq.gz/_ds_2.fq.gz}"
    # Run software
    echo "Downsampling the sequencing reads of sample ${sample}"
    rasusa reads \
    --coverage ${coverage} \
    --genome-size ${genomesize_mb}mb \
    -s 100 \
    -o "3_fastp/${output_r1}" \
    -o "3_fastp/${output_r2}" \
    "${r1}" "${r2}"
    # Deactivate Conda environment
    conda activate base

    # Delete the original trimmed reads files
    rm "${r1}" "${r2}"
done
# Compress output files
zip -r 3_fastp_genomesize.zip 3_fastp/*_genomesize
# Delete output files
rm -r 3_fastp/*_genomesize


############################################################
## 4) Trimmed reads quality assessment
############################################################

############################################################
## FastQC

echo "############################################################"
echo "# Running FastQC for trimmed reads"
echo "############################################################"

# Create an output directory
mkdir 4_fastqc
# Activate Conda environment
conda activate fastqc
# Run FastQC
fastqc -t $(nproc --ignore=1) 3_fastp/*.fq.gz -o 4_fastqc
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 4_fastqc.zip 4_fastqc

############################################################
## Fastp -> FastQC -> MultiQC

echo "############################################################"
echo "# Running for MultiQC for trimmed reads"
echo "############################################################"

# Activate Conda environment
conda activate multiqc
# Run MultiQC
multiqc 4_fastqc/*_fastqc.zip -o 4_fastqc_multiqc
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 4_fastqc_multiqc.zip 4_fastqc_multiqc
# Delete the output directory
rm -r 4_fastqc 4_fastqc_multiqc 


############################################################
## 5) De novo assembly
############################################################

echo "############################################################"
echo "# 5) Running de novo assembly with Unicycler"
echo "############################################################"

############################################################
## Unicycler (only for prokaryotic genomes)

# Create an output directory
mkdir 5_unicycler
# Activate Conda environment
conda activate unicycler
# Loop through a list of files
for r1 in 3_fastp/*_1.fq.gz; do
    # Extract r2 path
    r2=${r1/_1.fq.gz/_2.fq.gz}
    # Extract r1 file name
    filename=${r1##*/}
    # Extract sample name
    sample=${filename%%_*}
    # Verify if the assembly files are empty
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 5_unicycler_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi
    # Run Unicycler
    unicycler \
    -t $(nproc --ignore=1) \
    --spades_options "--cov-cutoff auto" \
    --min_fasta_length 200 \
    -1 "${r1}" \
    -2 "${r2}" \
    -o "5_unicycler/${sample}_unicycler"
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 5_unicycler.zip 5_unicycler


############################################################
## 6) Organization of de novo assembly files
############################################################

############################################################
## Create assemblies directory

echo "############################################################"
echo "# Creating assemblies directory"
echo "############################################################"

mkdir 6_assemblies

############################################################
## Assemblies from Unicycler

# Loop through a list of directories
for dir in 5_unicycler/*/; do
    # Extract sample name from directory name
    # sample=$(basename "$directory" _unicycler)
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_unicycler*}
    # Copy and rename the assembly file
    cp "${dir}assembly.fasta" "6_assemblies/${sample}.fasta"
done

############################################################
## Compress the assemblies directory

zip -r 6_assemblies.zip 6_assemblies

############################################################
## Delete the assembler output directory
rm -r 5_unicycler


############################################################
## 7) Assembly quality assessment
############################################################

############################################################
## CheckM2

echo "############################################################"
echo "# Running CheckM2"
echo "############################################################"

# Activate Conda environment
conda activate checkm2
# Run the program
checkm2 predict \
--threads $(nproc --ignore=1) \
-x fasta \
--input 6_assemblies \
--output-directory 7_checkm
# Deactivate Conda environment
conda activate base
# Copy and rename the output file
cp 7_checkm/quality_report.tsv 7_checkm2.tsv
# Delete the output directory
rm -r 7_checkm

############################################################
## GUNC

echo "############################################################"
echo "# Running GUNC"
echo "############################################################"

# Create an output directory
mkdir 7_gunc 7_gunc_temp
# Activate Conda environment
conda activate gunc
# Run the program
gunc run \
--threads $(nproc --ignore=1) \
--contig_taxonomy_output \
--file_suffix .fasta \
--input_dir 6_assemblies \
--temp_dir 7_gunc_temp \
--out_dir 7_gunc
# Copy and rename the output file
cp 7_gunc/*maxCSS_level.tsv 7_gunc.tsv
# Plotting the data
# Create an output directory
mkdir 7_gunc/plot
# Loop through a list of files
for file in 7_gunc/diamond_output/*.out; do\
   gunc plot \
   --contig_display_num 0\
   --diamond_file $file \
   --out_dir 7_gunc/plot;\
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 7_gunc.zip 7_gunc
# Delete the output directory
rm -r 7_gunc 7_gunc_temp 

############################################################
## QUAST

echo "############################################################"
echo "# Running QUAST"
echo "############################################################"

# Activate Conda environment
conda activate quast
# Run the program
quast.py -m 0 -o 7_quast 6_assemblies/*.fasta
# Deactivate Conda environment
conda activate base
cp 7_quast/transposed_report.tsv 7_quast.tsv
# Compress the output directory
zip -r 7_quast.zip 7_quast
# Delete the output directory
rm -r 7_quast

############################################################
## Barrnap

echo "############################################################"
echo "# Running Barrnap"
echo "############################################################"

# Create an output directory
mkdir 7_barrnap
# Activate Conda environment
conda activate barrnap
# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    # Extract sample name
    prefix=$(basename ${file} .fasta)
    # Run barrnap for bacteria
    barrnap \
    --threads $(nproc --ignore=1) \
    --kingdom bac \
    -o "7_barrnap/${prefix}_bac_barrnap.fasta" \
    < $file \
    > "7_barrnap/${prefix}_bac_barrnap.gff"
    # Run the program for archaea
    barrnap \
    --threads $(nproc --ignore=1) \
    --kingdom arc \
    -o "7_barrnap/${prefix}_arc_barrnap.fasta" \
    < $file \
    > "7_barrnap/${prefix}_arc_barrnap.gff"
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 7_barrnap.zip 7_barrnap
# Delete the output directory
rm -r 7_barrnap

############################################################
## Vertical sequencing coverage

echo "############################################################"
echo "# Calculating vertical sequencing coverage"
echo "############################################################"

# Create output file
echo -e Sample"\t"Coverage > 7_coverage.tsv
# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    # Extract assembly name
    assembly=$(basename $file .fasta)
    # Extract assembly name
    sample=${assembly%%_*}
    # Inform the sample
    echo -e Calculating sequencing coverage for assembly: $assembly
    echo -e Assembly file: ${file}
    echo -e R1 file: 3_fastp/${sample}_*_1.fq.gz
    echo -e R2 file: 3_fastp/${sample}_*_2.fq.gz
    # Count bases in assembly
    bases_in_assembly=$(grep -v '^>' "$file" | tr -d '\n' | wc -c)
    echo -e Bases in assembly: $bases_in_assembly
    # Count bases in sequencing files
    # bases_in_reads=$(zcat 3_fastp/"$sample"*.gz | awk 'NR%4==2 {print $0}' | tr -d '\n' | wc -c)
    bases_in_reads=$(zcat 3_fastp/${sample}_*_1.fq.gz 3_fastp/${sample}_*_2.fq.gz | 
                     awk 'NR%4==2 {print length}' | paste -sd+ | bc)
    echo -e Bases in reads: $bases_in_reads
    # Calculate vertical coverage
    if [ "$bases_in_assembly" -gt 0 ]; then
        # coverage=$(($bases_in_reads / $bases_in_assembly))
        coverage=$(echo "scale=2; $bases_in_reads / $bases_in_assembly" | bc)
        echo -e Coverage: "${coverage}\n"
        echo -e "${assembly}\t${coverage}" >> 7_coverage.tsv
    else
        echo -e "${assembly}\t0" >> 7_coverage.tsv
    fi
done


############################################################
## 8) Taxonomic assignment
############################################################

############################################################
## GTDB-Tk

echo "############################################################"
echo "# Running GTDB-Tk"
echo "############################################################"

# Requires 64GB of RAM if the species is not identified by the ANI screening step
# Activate Conda environment
conda activate gtdbtk
# Run the program
gtdbtk classify_wf \
--cpus $(nproc --ignore=1) \
--extension .fasta \
--mash_db /db/gtdbtk \
--genome_dir 6_assemblies \
--out_dir 8_gtdbtk
# Copy and rename output files
if [ -f "8_gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
    cp 8_gtdbtk/gtdbtk.bac120.summary.tsv 8_gtdbtk_bacteria.tsv
fi
if [ -f "8_gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
    cp 8_gtdbtk/gtdbtk.ar53.summary.tsv 8_gtdbtk_archaea.tsv
fi
# Compress the output directory
zip -r 8_gtdbtk.zip 8_gtdbtk
# Delete the output directory
rm -r 8_gtdbtk
# Deactivate Conda environment
conda activate base

############################################################
## TYGS

# https://tygs.dsmz.de/user_requests/new
# Send the files from 6_assemblies


############################################################
## 9) Plasmid identification
############################################################

############################################################
## MOB-suite

echo "############################################################"
echo "# Running MOB-suite"
echo "############################################################"

# Create an output directory
mkdir 9_mobsuite
# Activate Conda environment
conda activate mob_suite
# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    # Extract sample name
    sample=$(basename $file .fasta)
    # Run the program
    mob_recon -n $(nproc --ignore=1) \
    --infile $file \
    --outdir 9_mobsuite/${sample}_mobsuite
done
# Deactivate Conda environment
conda activate base

# Merge contig_report.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/contig_report_all.tsv ]; then
    rm 9_mobsuite/contig_report_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any contig_report.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "contig_report.txt" | grep -q .; then
    # Iterate over all found contig_report.txt files
    for file in 9_mobsuite/*/contig_report.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/contig_report_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/contig_report_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/contig_report_all.tsv
    done
fi

# Merge mobtyper_results.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/mobtyper_results.txt ]; then
    rm 9_mobsuite/mobtyper_results_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any mobtyper_results.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "mobtyper_results.txt" | grep -q .; then
    # Iterate over all found mobtyper_results.txt files
    for file in 9_mobsuite/*/mobtyper_results.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/mobtyper_results_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/mobtyper_results_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/mobtyper_results_all.tsv
    done
fi

# Merge mge.report.txt files
# Delete merged file if it exists
if [ -f 9_mobsuite/mge.report_all.tsv ]; then
    rm 9_mobsuite/mge.report_all.tsv
fi
# Initialize control variable to check in the header was printed
header_printed=0
# Check if any mge.report.txt files exist
if find 9_mobsuite/ -maxdepth 2 -type f -name "mge.report.txt" | grep -q .; then
    # Iterate over all found mge.report.txt files
    for file in 9_mobsuite/*/mge.report.txt; do
        # Test if the header was not printed yet
        if [ $header_printed -eq 0 ]; then
            # If not, contatenate the entire file
            cat "$file" >> 9_mobsuite/mge.report_all.tsv
            # Mark that the header was printed
            header_printed=1
        else
        # The header was printed, so contanate file and ignore its header
        tail -n +2 "$file" >> 9_mobsuite/mge.report_all.tsv  
        fi
        # Add a new line to separate the results of each sample
        # echo >> 9_mobsuite/mge.report_all.tsv
    done
fi

# Copy merged mobtyper result file to main directory
if ls 9_mobsuite/mobtyper_results_all.tsv > /dev/null 2>&1; then
    cp 9_mobsuite/mobtyper_results_all.tsv 9_mobsuite_mobtyper_results_all.tsv
fi

# Rename MOB-suite fasta files
# Loop through a list of directories
for dir in 9_mobsuite/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_mobsuite*}
    # Go to input directory
    # Add sample name to the chromosome file name 
    mv ${dir}chromosome.fasta ${dir}"$sample"_chromosome.fasta
    # Add sample name to the plasmid file name 
    rename s/plasmid_/"$sample"_plasmid_/ ${dir}plasmid*
    # Leave input diretory
done


#########################################################################
## 10) Assignment of contigs to molecules
############################################################

############################################################
## Add molecule attribution to contigs. Required for batch genome submission.

echo "############################################################"
echo "# Adding molecule attribution to contigs"
echo "############################################################"

# Copy assemblies directory
cp -r 6_assemblies 10_assemblies_for_analysis
# Change extension to .fsa
rename "s/.fasta$/.fsa/" 10_assemblies_for_analysis/*.fasta
# Loop through a list of files (assemblies)
for assembly in 10_assemblies_for_analysis/*.fsa; do
    # Extract file name
    filename=${assembly##*/}
    # Extract sample name
    sample=${filename%%.*}
    # Inform sample
    echo "Analyzing $sample"
    # Loop through a list of files (mob-suite result files)
    awk '1' 9_mobsuite/"$sample"_mobsuite/contig_report.txt | while IFS=$'\t' read -r sample_id molecule_type primary_cluster_id secondary_cluster_id contig_id others; do
        # Ignore column names
        if [ "$sample_id" == "sample_id" ]; then
            # Skip to the next iteration
            continue
        fi
        # Conditional: plasmid contig
        if [ "$molecule_type" == "plasmid" ]; then
            # Create new contig header
            new_contig_id=$(echo "$contig_id" "[plasmid-name="p"$primary_cluster_id""]")
            # Inform sample and contig
            echo Sample: "$sample" - Contig: "$new_contig_id"
            # Replace old header with the new one
            sed -i "s/^>${contig_id}/>${new_contig_id}/" $assembly
        fi
        # Conditional: chromosome contig
        if [ "$molecule_type" == "chromosome" ]; then
            # Create new contig header
            new_contig_id=$(echo "$contig_id" "[chromosome=1]")
            # Inform sample and contig
            echo Sample: "$sample" - Contig: "$new_contig_id"
            # Replace old header with the new one
            sed -i "s/^>${contig_id}/>${new_contig_id}/" $assembly
        fi
    done
done
# Compress the output directory
zip -r 10_assemblies_for_analysis.zip 10_assemblies_for_analysis
# Compress mob-suite directory
zip -r 9_mobsuite.zip 9_mobsuite
# Delete output file
rm -r 9_mobsuite

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
