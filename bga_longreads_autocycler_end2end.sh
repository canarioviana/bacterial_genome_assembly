#!/bin/bash
# Bash script for bacterial genome assembly from long-read sequencing data
#
# Author: Marcus Vinicius Canário Viana
# Date: 24/03/2026
# Repository: https://github.com/canarioviana/bacterial_genome_assembly
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
# Place this script (**bga_longreads_autocycler_end2end.sh**) in the working directory and execute it **using the following commands**:
# chmod +x bga_longreads_autocycler_end2end.sh 
# ./bga_longreads_autocycler_end2end.sh <sequencing_type> <medaka_model>
# Valid options for sequencing type: ont_r9 | ont_r10 | pacbio_clr | pacbio_hifi
# For valid options for medaka models execute:
# conda activate medaka
# medaka tools list_models
# The information to choose the model is in the fastq headers


############################################################
## SUMMARY OF END-TO-END GENOME ASSEMBLY WORKFLOW FROM LONG-READS
############################################################

## 0) Sequencing type, error handling and checking Conda installation
## 1) Sequencing reads directory and files
    # Reads stored as local files
## 2) Raw reads quality assessment
    # NanoPlot
## 3) Raw reads trimming and estimation of genome size
    # Chopper
## 4) Trimmed reads quality assessment
    # NanoPlot
## 5) De novo assembly
    # Autocycler -> dnaapler 
    # Medaka (optional)
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
## 0) Sequencing type, Medaka model, error handling and checking Conda installation
############################################################

############################################################
## Sequencing type

# Access the sequencing type from the command line
seq_type=$1

# Valid options (based on Flye assembler)
valid_options="ont_r9 | ont_r10 | pacbio_clr | pacbio_hifi"

# Check if the argument was provided
if [ -z "$seq_type" ]; then
    echo "Error: You must provide the sequencing type."
    echo "Options: $valid_options"
    echo "Usage: $0 <sequencing_type>"
    exit 1
fi

# Check if the sequencing type argument was provided
case "$seq_type" in
    # Match against all valid options
    ont_r9|ont_r10|pacbio_clr|pacbio_hifi)
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

############################################################
## Medaka model

# Access Medaka model from the command line
model=$2

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

echo "############################################################"
echo "# Checking local read files"
echo "############################################################"

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

echo "############################################################"
echo "# Running NanoPlot for raw reads"
echo "############################################################"

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

    # Verify if the input files are empty
    if [ ! -s "$reads" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 2_nanoplot_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    echo "NanoPlot is processing sample: ${sample}"
    echo "A warning is expected if the Conda environment was installed by the root."
    mkdir "2_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "2_nanoplot/${sample}"
done
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -r 2_nanoplot.zip 2_nanoplot
# Delete the output directory
rm -r 2_nanoplot


############################################################
## 3) Raw reads trimming and estimation of genome size
############################################################

###########################################################
## chopper

echo "############################################################"
echo "# Running chopper"
echo "############################################################"

# Create an output directory
mkdir 3_chopper
# Activate Conda environment
conda activate chopper
# Loop through a list of files
for reads in 1_reads/*.fq.gz; do

    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Verify if the input files are empty
    if [ ! -s "$reads" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 3_chopper_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    # Inform the current sample being processed
    echo "Processing sample: ${sample}"
    echo "Reads file ${reads}"

    # Run chopper
    chopper \
    -t $(nproc --ignore=1) \
    --minlength 500 \
    --quality 10 \
    --trim-approach trim-by-quality \
    --cutoff 10 \
    -i "$reads" \
    | gzip > 3_chopper/"${sample}_trimmed.fq.gz"
done
# Deactivate Conda environment
conda deactivate


############################################################
## 4) Trimmed reads quality assessment
############################################################

############################################################
## NanoPlot

echo "############################################################"
echo "# Running NanoPlot for trimmed reads"
echo "############################################################"

# Create an output directory
mkdir 4_nanoplot
# Activate Conda environment
conda activate nanoplot
# Loop through a list of files
for reads in 3_chopper/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Verify if the input files are empty
    if [ ! -s "$reads" ]; then
        echo "Warning: The input files of sample ${sample} are empty. Skiping..."
        echo -e "${sample}" >> 4_nanoplot_skiped_samples.tsv
        # Skip to the next iteration
        continue
    fi

    echo "NanoPlot is processing sample: ${sample}"
    echo "A warning is expected if the Conda environment was installed by the root."
    mkdir "4_nanoplot/${sample}"
    NanoPlot -t $(nproc --ignore=1) --fastq ${reads} -p "${sample}_" -o "4_nanoplot/${sample}"
done
# Deactivate Conda environment
conda deactivate
# Compress the output directory
zip -r 4_nanoplot.zip 4_nanoplot
# Delete the output directory
rm -r 4_nanoplot


############################################################
## 5) De novo assembly
############################################################

############################################################
## Autocycler
# https://github.com/rrwick/Autocycler/wiki/Illustrated-pipeline-overview

echo "############################################################"
echo "# Running de novo assembly with Autocycler"
echo "############################################################"

# Create an output directory
mkdir 5_autocycler

# Activate Conda environment
conda activate autocycler
# Loop through a list of files
for reads in 3_chopper/*.fq.gz; do
    # Extract reads file name
    readsfilename=${reads##*/}
    # Extract sample name
    sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
    sample=${sample%%.*} # In case the file name is samplename.fq.gz

    # Create sample output diretory
    mkdir -p 5_autocycler/${sample}_autocycler
    # Go to sample output diretory
    cd 5_autocycler/${sample}_autocycler

    # Estimate genome size
    genome_size=$(autocycler helper genome_size --reads ../../"$reads" --threads $(nproc --ignore=1))

    # Step 1: Autocyler subsample
    autocycler subsample --reads ../../"$reads" --out_dir 1_subsampled_reads --genome_size "$genome_size"

    # Step 2: Generate input assemblies
    mkdir 2_assemblies
    for assembler in flye metamdbg miniasm necat nextdenovo plassembler raven; do
        for i in 01 02 03 04; do
            autocycler helper \
                "$assembler" \
                --reads 1_subsampled_reads/sample_"$i".fastq \
                --out_prefix 2_assemblies/"$assembler"_"$i" \
                --threads $(nproc --ignore=1) \
                --genome_size "$genome_size" \
                --read_type "$seq_type"
        done
    done
    # Delete the subsampled reads files
    rm 1_subsampled_reads/*.fastq

    # Step 3: Autocycler compress
    autocycler compress -i 2_assemblies -a 3_autocycler

    # Step 4: Autocycler cluster
    autocycler cluster -a 3_autocycler

    # Steps 5 and 6: Autocycler trim and Autocycler resolve
    for file in 3_autocycler/clustering/qc_pass/cluster_*; do
        autocycler trim -c "$file"
        autocycler resolve -c "$file"
    done

    # Step 7: Autocycler combine
    autocycler combine -a 3_autocycler -i 3_autocycler/clustering/qc_pass/cluster_*/5_final.gfa
   
    # Reorient circular sequences with Dnaapler
    # conda activate dnaapler # Uncomment if dnaapler does not work in the autocycler Conda environment
    dnaapler all -i 3_autocycler/consensus_assembly.gfa -o 4_dnaapler -t $(nproc --ignore=1)
    # conda activate autocycler # Uncomment if dnaapler does not work in the autocycler Conda environment
    autocycler gfa2fasta -i 4_dnaapler/dnaapler_reoriented.gfa -o 4_dnaapler/dnaapler_reoriented.fasta

    # Back to main directory
    cd ../../
done
# Deactivate Conda environment
conda deactivate

# Create metrics file with a header
echo -e "name\tinput_read_count\tinput_read_bases\tinput_read_n50\tpass_cluster_count\tfail_cluster_count\toverall_clustering_score\tuntrimmed_cluster_size\tuntrimmed_cluster_distance\ttrimmed_cluster_size\ttrimmed_cluster_median\ttrimmed_cluster_mad\tconsensus_assembly_bases\tconsensus_assembly_unitigs\tconsensus_assembly_fully_resolved" > 5_autocycler_metrics.tsv
# Activate Conda environment
conda activate autocycler
# Extract metrics from each assembly
for dir in 5_autocycler/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_autocycler*}
    # Add metrics from each assembly
    autocycler table -a "${dir}" -n "$sample" >> 5_autocycler_metrics.tsv
done
# Deactivate Conda environment
conda deactivate

# Compress output files
# zip -r 5_autocycler.zip 5_autocycler


############################################################
## Autocycler -> Medaka

case "$seq_type" in
    # Check if the sequencing type is NanoPore
    ont_r9|ont_r10)

    # Check if argument $model is not empty
    if [ -n "$model" ]; then
        # Create an output directory
        mkdir 5_autocycler_medaka

        # Avoid using the GPU
        export CUDA_VISIBLE_DEVICES=""

        # Activate Conda environment
        conda activate medaka

        # Loop through a list of files
        for reads in 3_chopper/*.fq.gz; do
            # Extract reads file name
            readsfilename=${reads##*/}
            # Extract sample name
            sample=${readsfilename%%_*} # In case the file name is samplename_*.fq.gz
            sample=${sample%%.*} # In case the file name is samplename.fq.gz

            # Run medaka
            medaka_consensus \
                -i ${reads} \
                -d 5_autocycler/"${sample}"_autocycler/4_dnaapler/dnaapler_reoriented.fasta \
                -o 5_autocycler_medaka/"${sample}"_medaka \
                -t $(nproc --ignore=1) \
                -m "${model}"
        done
    else
        echo "Skipping Medaka: seq_type is ont_r9 or ont_r10, but the model was not informed." | tee 5_autocycler_medaka_skipped.log
    fi
    ;;
    *)
        echo "Skipping Medaka: seq_type '$seq_type' is not ont_r9 or ont_r10." | tee 5_autocycler_medaka_skipped.log
        ;;
esac

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
## Assemblies from Autocycler

# Loop through a list of directories
for dir in 5_autocycler/*/; do
    # Extract directory name
    dirname=${dir#*/}
    # Extract sample name
    sample=${dirname%%_autocycler*}
    # Copy and rename the assembly file
    cp "${dir}4_dnaapler/dnaapler_reoriented.fasta" "6_assemblies/${sample}.fasta"
done

############################################################
## Assemblies from Autocycler -> Medaka

# Loop through a list of directories
for dir in 5_autocycler_medaka/*/; do
    # Check if the directory actually exists (prevents errors if no directories are found)
    [ -d "$dir" ] || continue

    # Execute only if the input file exists
    if [ -f "${dir}consensus.fasta" ]; then
        # Extract directory name
        dirname=${dir#*/}
        # Extract sample name
        sample=${dirname%%_medaka*}
        # Copy and rename the assembly file
        cp "${dir}consensus.fasta" "6_assemblies/${sample}_medaka.fasta"
    fi
done


############################################################
## Compress the assemblies directory

# zip -r 6_assemblies.zip 6_assemblies

############################################################
## Delete the assemblers output directories

# rm -r 5_autocycler


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
conda deactivate
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
conda deactivate
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
conda deactivate
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
conda deactivate
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
    echo -e Read file: 3_chopper/${sample}_trimmed.fq.gz
    # Count bases in assembly
    bases_in_assembly=$(grep -v '^>' "$file" | tr -d '\n' | wc -c)
    echo -e Bases in assembly: $bases_in_assembly
    # Count bases in sequencing files
    # bases_in_reads=$(zcat 3_fastp/"$sample"*.gz | awk 'NR%4==2 {print $0}' | tr -d '\n' | wc -c)
    bases_in_reads=$(zcat 3_chopper/${sample}_trimmed.fq.gz | 
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
conda deactivate

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
conda deactivate

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

