#!/bin/bash
# Bash script for bacterial genome assembly from long-read sequencing data
#
# Author: Marcus Vinicius Canário Viana
# Date: 20/10/2025
# More info: see README.md in the repository
#
# Instructions:
#
# **A. Using Local Read Files**
#
# 1. Standardize the paired-end file names of each sample to samplename.fq.gz or samplename_*.fq.gz
# 2. In the working directory, create the directory **1_reads** and place the read files **inside it**.
# 
# **B. Execution**
#
# Place this script (**bga_longreads_end2end.sh**) in the working directory and execute it **using the following commands**:
# chmod +x bga_longreads_end2end.sh 
# ./bga_longreads_end2end.sh <sequencing_type>
# Valid options for sequencing type: --pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw | --nano-corr | --nano-hq


############################################################
## SUMMARY OF END-TO-END GENOME ASSEMBLY WORKFLOW FROM LONG-READS
############################################################

## 0) Sequencing type, error handling and checking Conda installation
## 1) Sequencing reads directory and files
  # Reads stored as local files
## 2) Raw reads quality assessment
    # NanoPlot
## 3) Raw reads trimming and estimation of genome size
    # Fastplong
    # Estimation of genome size (KMC and GenomeScope)
## 4) Trimmed reads quality assessment
    ## NanoPlot
## 5) De novo assembly
    # Flye
    # metaMDBG (Not used in the end-to-end pipeline. Unsatisfactory using a test dataset.)
    # myloasm (Not used in the end-to-end pipeline. Unsatisfactory using a test dataset.)
    # NextDenovo (Not used in the end-to-end pipeline. Unsatisfactory using a test dataset.)
    # Raven
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
    # MOB-suite and an inhouse script


############################################################
## 0) Sequencing type, error handling and checking Conda installation
############################################################

############################################################
## Sequencing type

# Access the sequencing type from the command line
seq_type=$1

# Valid options (based on Flye assembler)
valid_options="--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw | --nano-corr | --nano-hq"

# Check if the argument was provided
if [ -z "$seq_type" ]; then
    echo "Error: You must provide the sequencing type."
    echo "Options: --pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw | --nano-corr | --nano-hq"
    echo "Usage: $0 <sequencing_type>"
    exit 1
fi

# Check if the sequencing type argument was provided
case "$seq_type" in
    # Match against all valid options
    --pacbio-raw|--pacbio-corr|--pacbio-hifi|--nano-raw|--nano-corr|--nano-hq)
        # Argument is valid, now proceed with the script logic
        echo "Valid sequencing type selected: $seq_type"
        ;;
    # Handle any other input
    *)
        echo "Error: '$seq_type' is not a valid sequencing type."
        echo "Options: $valid_options"
        exit 1
        ;;
esac

# Adjust the sequencing type argument from Flye to metaMDBG
# --in-hifi (Will be used in place of --pacbio-raw | --pacbio-corr | --pacbio-hifi)
# --in-ont (Will be used in place of --nano-raw | --nano-corr | --nano-hq)
case "$seq_type" in
    # PacBio options
    --pacbio-raw|--pacbio-corr|--pacbio-hifi)
    seq_type_metamdbg="--in-hifi"
        ;;

    # Nanopore options
    --nano-raw|--nano-corr|--nano-hq)
        seq_type_metamdbg="--in-ont"
        ;;
esac

# Adjust the sequencing type argument from Flye to myloasm
# --nano-r10 (Will be used in place of --nano-hq)
# --nano-r9 (Will be used in place of --nano-raw | --nano-corr )
# --hifi (Will be used in place of --pacbio-corr | --pacbio-hifi)
case "$seq_type" in
    # PacBio options
    --pacbio-raw|--pacbio-corr|--pacbio-hifi)
    seq_type_myloasm="--hifi"
        ;;
        
    # Nanopore options
    --nano-raw|--nano-corr)
        seq_type_myloasm="--nano-r9"
        ;;
    --nano-hq)
        seq_type_myloasm="--nano-r10"
        ;;
esac

# Adjust the sequencing type argument from Flye to NextDenovo
# clr (Will be used in place of --pacbio-raw | --pacbio-corr)
# hifi (Will be used in place of --pacbio-hifi)
# ont (Will be used in place of --nano-raw | --nano-corr | --nano-hq)
case "$seq_type" in
    # PacBio options
    --pacbio-raw|--pacbio-corr)
    seq_type_nextdenovo="clr"
        ;;
    --pacbio-hifi)
        seq_type_nextdenovo="hifi"
        ;;

    # Nanopore options
    --nano-raw | --nano-corr | --nano-hq)
        seq_type_nextdenovo="ont"
        ;;
esac

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
## Reads stored as local files

############################################################
## Checking reads directory

# Checking wether the directory 1_reads exists
if [ ! -d 1_reads ]; then
    echo "The reads directory '1_reads' was not found" >&2
    echo "Please create it and put the files in it" >&2
    exit 1
fi

############################################################
## Standardize the file names of each sample from *.fastq.gz to *.fq.gz

# Checking the presence of files in the format 1_reads/*.fastq.gz
if ls 1_reads/*.fastq.gz > /dev/null 2>&1; then
    echo "Found files in the format 1_reads/*.fastq.gz" >&2
    echo "Renaming them to the format *.fq.gz" >&2
    # Rename the .fq.gz files
    rename 's/\.fastq\.gz/.fq.gz/' 1_reads/*.fastq.gz

    # Checking the presence of .md5 files
    if ls 1_reads/*.fastq.gz.md5 > /dev/null 2>&1; then
        echo "Found files in the format 1_reads/*.fastq.gz.md5" >&2
        echo "Renaming them to the format *.fq.gz.md5" >&2
        # Rename the file names inside the .md5 files
        sed -i 's/.fastq\.gz/.fq.gz/' 1_reads/*.md5
        # Rename the .md5 files
        rename 's/\.fastq\.gz\.md5/.fq.gz.md5/' 1_reads/*.fastq.gz.md5
    fi

else
    :
fi

############################################################
## Checking the presence and names of the input files

files_found=0
#Loop
for reads in 1_reads/*.fq.gz; do
    #Check the presence of reads file
    if [[ "$reads" == "1_reads/*.fq.gz" ]]; then
        break 
    fi
    # If at least one *.fq.gz file was found
    files_found=1
done

# Final reports
if [ "$files_found" -eq 0 ]; then
    echo "No read files in the format '*.fq.gz' was found in the diretory '1_reads'" >&2
    exit 1
fi


############################################################
## 2) Raw reads quality assessment
############################################################

############################################################
## NanoPlot

# Create an output directory
mkdir 2_nanoplot
# Activate Conda environment
conda activate nanoplot
# Loop through a list of files
for reads in 1_reads/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz
    echo "NanoPlot is processing sample: ${sample}"
    echo "A warning is expected if the Conda environment was installed by the root."
    mkdir "2_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "2_nanoplot/${sample}"
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 2_nanoplot.zip 2_nanoplot
# Delete the output directory
rm -r 2_nanoplot


############################################################
## 3) Raw reads trimming and estimation of genome size
############################################################

############################################################
## Fastplong

# Create an output directory
mkdir 3_fastplong
# Activate Conda environment
conda activate fastplong
# Loop through a list of files
for reads in 1_reads/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz
    # Inform the current sample being processed
    echo "Processing sample: ${sample}"
    echo "Reads file ${reads}"
    # Run Fastplong
    fastplong \
    --thread $(nproc --ignore=1) \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --length_required 1000 \
    -i "$reads" \
    -o 3_fastplong/"${sample}_trimmed.fq.gz" \
    --html 3_fastplong/"$sample"_fastp.html \
    --json 3_fastplong/"$sample"_fastp.json
done
# Deactivate Conda environment
conda activate base
# Compress report files
zip -r 3_fastplong.zip 3_fastplong/*.json 3_fastplong/*.html
# Delete report files
rm 3_fastplong/*.json 3_fastplong/*.html

############################################################
## Estimation of genome size

# Create output directory
mkdir 3_genomesize
# Create output table file
> 3_genomesize.tsv

# Loop through a list of files
for reads in 3_fastplong/*.fq.gz; do
    # Extract file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*}

    # Create directory for the genome size estimation results
    genomesizedir="3_genomesize/${sample}_genomesize"
    mkdir ${genomesizedir}
    
    # Count k-mers using KMC
    # Activate Conda environment
    conda activate kmc
    # Count k-mers
    echo "Counting k-mers from the sequencing reads of sample ${sample}"
    # Create kmc temporary directory
    mkdir kmc_tmp
    # Create kmc list of input read files
    ls -1 ${reads} > kmc_input_reads.txt
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

    # Get estimated genome size
    genomesize_bp=$(grep "Genome Haploid Length" "${genomesizedir}/genomescope/summary.txt" | awk '{print $(NF-1)}' | tr -d ',')
    genomesize_mb=$(echo "scale=2; $genomesize_bp / 1000000" | bc)
    # Inform estimated genome size
    echo "Estimated genome size of sample ${sample}: ${genomesize_mb}Mb"
    echo -e "${sample}\t${genomesize_bp}" >> 3_genomesize.tsv 

    # Move kmc temporary files to the genome size directory
    mv kmc_count* kmc_histogram.tsv kmc_input_reads.txt ${genomesizedir}
    # Delete the directory kmc_tmp
    rm -r kmc_tmp

    # Compress and delete output directory
    zip -r 3_genomesize.zip 3_genomesize
    rm -r 3_genomesize
done


############################################################
## 4) Trimmed reads quality assessment
############################################################

############################################################
## NanoPlot

# Create an output directory
mkdir 4_nanoplot
# Activate Conda environment
conda activate nanoplot
# Loop through a list of files
for reads in 3_fastplong/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz
    echo "NanoPlot is processing sample: ${sample}"
    echo "A warning is expected if the Conda environment was installed by the root."
    mkdir "4_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "4_nanoplot/${sample}"
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 4_nanoplot.zip 4_nanoplot
# Delete the output directory
rm -r 4_nanoplot


############################################################
## 5) De novo assembly
############################################################

############################################################
## Flye

# Create an output directory
mkdir 5_flye

# Activate Conda environment
conda activate flye
# Loop through a list of files
for reads in 3_fastplong/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*}
    # Run Flye
    echo "Assembling the genome of sample: $sample"
    flye -t $(nproc --ignore=1) \
    "${seq_type}" \
    "$reads" \
    -o 5_flye/"$sample"_flye
done
# Deactivate Conda environment
conda activate base
# Compress output files
zip -r 5_flye.zip 5_flye

# ############################################################
# ## metaMDBG (for metagenomes)

# # Create an output directory
# mkdir 5_metamdbg
# # Activate Conda environment
# conda activate metamdbg
# # Loop through a list of files
# for reads in 3_fastplong/*.fq.gz; do
#     # Extract reads file name
#     readsfilename=${reads##*/}
#     # Extract sample name
#     sample=${readsfilename%%_*}
#     # Run metaMDBG
#     echo "metaMDBG is assembling the genome of sample: $sample"
#     metaMDBG asm \
#     --threads $(nproc --ignore=1) \
#     "$seq_type_metamdbg" "$reads" \
#     --out-dir "5_metamdbg/${sample}_metamdbg"
# done
# # Deactivate Conda environment
# conda activate base
# # Compress the output directory
# zip -r 5_metamdbg.zip 5_metamdbg

# ############################################################
# ## myloasm

# # Create an output directory
# mkdir 5_myloasm
# # Activate Conda environment
# conda activate myloasm
# # Loop through a list of files
# for reads in 3_fastplong/*.fq.gz; do
#     # Extract reads file name
#     readsfilename=${reads##*/}
#     # Extract sample name
#     sample=${readsfilename%%_*}
#     # Run Myloasm
#     echo "Myloasm is assembling the genome of sample: $sample"
#     myloasm \
#     -t $(nproc --ignore=1) \
#     -o "5_myloasm/${sample}_myloasm" \
#     "$reads" "$seq_type_myloasm"
# done
# # Deactivate Conda environment
# conda activate base
# # Compress the output directory
# zip -r 5_myloasm.zip 5_myloasm

# ############################################################
# ## NextDenovo

# # Create an output directory
# mkdir 5_nextdenovo

# # Activate Conda environment
# conda activate nextdenovo
# # Loop through a list of files
# for reads in 3_fastplong/*.fq.gz; do
#     # Extract reads file name
#     readsfilename=${reads##*/}
#     # Extract sample name
#     sample=${readsfilename%%_*}

#     # Create output directory
#     working_dir="5_nextdenovo/${sample}_nextdenovo"
#     mkdir "${working_dir}"

#     #Get estimated genome size from file 5_genomesize.tsv
#     genomesize_pb=$(awk -F '\t' -v target="${sample}" '$1 == target {print $2; exit}' "3_genomesize.tsv")
#     genomesize_mb=$(echo "scale=2; $genomesize_bp / 1000000" | bc)
    
#     # Create list of input files
#     #reads_list="$PWD/5_nextdenovo/${sample}_nextdenovo/${sample}.fofn"
#     #ls "$PWD/$reads" > ${reads_list}
#     reads_list="${sample}.fofn"
#     ls "$reads" > ${reads_list}
#     # Create run.cfg
#     #config_file="5_nextdenovo/${sample}_nextdenovo/${sample}.cfg"
#     config_file="${sample}.cfg"
#     find $CONDA_PREFIX -path "*/doc/run.cfg" -type f \
#     -exec cp {} ./"${config_file}" \; -quit

#     # Edit prefix
#     sed -i "s|^job_prefix =.*|job_prefix = ${sample}|" "${config_file}"
#     # Edit number of parallel jobs
#     sed -i "s|^parallel_jobs =.*|parallel_jobs = $(nproc --ignore=1)|" "${config_file}"
#     # Edit input_type: raw, corrected
#     sed -i "s|^input_type =.*|input_type = raw|" "${config_file}"
#     # Edit read type: clr, ont, hifi
#     sed -i "s|^read_type =.*|read_type = $seq_type_nextdenovo|" "${config_file}"
#     # Edit list of input files
#     sed -i "s|^input_fofn =.*|input_fofn = ${reads_list}|" "${config_file}"
#     # Edit working directory
#     sed -i "s|^workdir =.*|workdir = ${working_dir}|" "${config_file}"
#     # Edit estimated genome size
#     sed -i "s|^genome_size =.*|genome_size = ${genomesize_mb}m|" "${config_file}"
#     # Run NextDenovo
#     echo "NextDenovo is assembling the genome of sample: $sample"
#     nextDenovo "${config_file}"
#     # Move log file
#     mv pid* ${sample}.fofn ${sample}.cfg "${working_dir}"
# done
# # Deactivate Conda environment
# conda activate base
# # Compress the output directory
# zip -r 5_nextdenovo.zip 5_nextdenovo

############################################################
## Raven

# Create an output directory
mkdir 5_raven
# Activate Conda environment
conda activate raven
# Loop through a list of files
for reads in 3_fastplong/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*}
    # Create output directory
    mkdir "5_raven/${sample}_raven"
    # Run Raven
    echo "Raven is assembling the genome of sample: $sample"
    raven \
    -t $(nproc --ignore=1) \
    "$reads" \
    > "5_raven/${sample}_raven/${sample}_raven.fasta"
done
# Deactivate Conda environment
conda activate base
# Delete temporary file
rm raven.cereal
# Compress the output directory
zip -r 5_raven.zip 5_raven


############################################################
## 6) Organization of de novo assembly files
############################################################

############################################################
## Create assemblies directory
mkdir 6_assemblies

############################################################
## Assemblies from Flye

# Loop through a list of directories
for dir in 5_flye/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_flye*}
    # Copy and rename the assembly file
    cp "${dir}assembly.fasta" "6_assemblies/${sample}_flye.fasta"
done

# ############################################################
# ## Assemblies from metaMDBG

# # Loop through a list of directories
# for dir in 5_metamdbg/*/; do
#     # Extract directory name
#     dirname=${dir#*/}
#     # Extract sample name
#     sample=${dirname%%_metamdbg*}
#     # Copy and rename the assembly file
#     gunzip -c "${dir}contigs.fasta.gz" > "6_assemblies/${sample}_metamdbg.fasta"
# done

# ############################################################
# ## Assemblies from myloasm 

# # Loop through a list of directories
# for dir in 5_myloasm/*/; do
#     # Extract directory name
#     dirname=${dir#*/}
#     # Extract sample name
#     sample=${dirname%%_myloasm*}
#     # Copy and rename the assembly file
#     cp "${dir}assembly_primary.fa" "6_assemblies/${sample}_myloasm.fasta"
# done

# ############################################################
# ## Assemblies from NextDenovo

# # Loop through a list of directories
# for dir in 5_nextdenovo/*/; do
#     # Extract directory name
#     dirname=${dir#*/}
#     # Extract sample name
#     sample=${dirname%%_nextdenovo*}
#     # Copy and rename the assembly file
#     cp "${dir}/03.ctg_graph/nd.asm.fasta" "6_assemblies/${sample}_nextdenovo.fasta"
# done

############################################################
## Assemblies from Raven

cp 5_raven/*/*.fasta 6_assemblies

############################################################
## Compress the assemblies directory

zip -r 6_assemblies.zip 6_assemblies

############################################################
## Delete the assemblers output directories

rm -r 5_flye
# rm -r 5_metamdbg
# rm -r 5_myloasm
# rm -r 5_nextdenovo
rm -r 5_raven


############################################################
## 7) Assembly quality assessment
############################################################

############################################################
## CheckM2

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
    # Run barrnap for archaea
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
    echo -e Read file: 3_fastplong/${sample}_trimmed.fq.gz
    # Count bases in assembly
    bases_in_assembly=$(grep -v '^>' "$file" | tr -d '\n' | wc -c)
    echo -e Bases in assembly: $bases_in_assembly
    # Count bases in sequencing files
    # bases_in_reads=$(zcat 3_fastp/"$sample"*.gz | awk 'NR%4==2 {print $0}' | tr -d '\n' | wc -c)
    bases_in_reads=$(zcat 3_fastplong/${sample}_trimmed.fq.gz | 
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
    while IFS=$'\t' read -r sample_id molecule_type primary_cluster_id secondary_cluster_id contig_id others; do
        # Ignore column names
        if [ "$sample_id" == "sample_id" ]; then
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
    done < 9_mobsuite/"$sample"_mobsuite/contig_report.txt
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

