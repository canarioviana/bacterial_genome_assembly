# Bash script for bacterial genome annotation, including software installation
#
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.
#
# Author: Marcus Vinicius Canário Viana
# Date: 30/09/2025
# More info: see README.md in the repository



##########################################################################
## Software installation
##########################################################################

##########################################################################
# Prokka (Genome annotation)
conda create -n prokka -c bioconda prokka -y

##########################################################################
# PanViTa (Virulence and antimicrobial resistance prediction)
conda create -n panvita -y
conda activate panvita
conda install -c anaconda wget basemap -y
conda install seaborn pandas -y
conda install -c conda-forge matplotlib -y
cd
git clone https://github.com/dlnrodrigues/panvita.git
cd panvita
sed -i 's/VFDB_setA_pro/VFDB_setB_pro/' panvita.py #For full dataset
sed -i 's/vfdb_core/vfdb_full/g' panvita.py
python3 panvita.py -u
python3 panvita.py -h
conda deactivate

##########################################################################
# RGI (Antimicrobial resistance prediction)
conda create -n rgi -c bioconda rgi -y
sudo chmod 777 -R /usr/local/anaconda3/envs/rgi/lib/python3.8/site-packages/app/_db/.ncbirc #so all the users can run the program

##########################################################################
# CRISPRcasFinder local (CRISPR/Cas system prediction)
git clone https://github.com/dcouvin/CRISPRCasFinder.git
cd CRISPRCasFinder
conda env create -f ccf.environment.yml -n crisprcasfinder
conda activate crisprcasfinder
mamba install -c bioconda macsyfinder=2.1.2 -y
macsydata install -u CASFinder==3.1.0
conda deactivate

##########################################################################
# FastMLST (Multi-Locus Sequence Typing)
conda create -n fastmlst -c bioconda fastmlst -y
conda activate fastmlst
fastmlst -t 1 --update-mlst
conda deactivate

##########################################################################
# EggNOG-mapper (Gene functional annotation)
conda create -n eggnog-mapper -c bioconda eggnog-mapper=2.1.12 -y
conda activate eggnog-mapper
conda env config vars set EGGNOG_DATA_DIR="/db/eggnog/"
conda deactivate
conda activate eggnog-mapper
echo $EGGNOG_DATA_DIR
download_eggnog_data.py --data_dir /db/eggnog
conda deactivate

##########################################################################
# COG classifier (Gene functional annotation)
conda create -n cogclassifier -c bioconda -c conda-forge cogclassifier -y


#########################################################################
## 10) Genome annotation

# Getting the genome annotation files from GenBank
# During the genome submission process, choose the genome annotation with PGAP: "Annotate this prokaryotic genome in the NCBI Prokaryotic Annotation Pipeline (PGAP) before its release"
# In the genome submission portal, check wether the genome submission was processed
# https://submit.ncbi.nlm.nih.gov/subs/genome/
# Before the chosen release date GenBank will give you only the .gb annotation file.
# If you do not want to release the genome before your analysis, you can used the .gb file to generate the chromosone (.fsa), gene (.ffn) and protein (.faa) files 
# Download the .gb files and run the script gbk_to_fasta_faa_ffn.py 
# rename s/_.*/.gbk/ *.gb
# python gbk_to_fasta_faa_ffn.py *.gbk
# for file in *.gbk; do name=$(basename $file .gbk); mkdir $name; done
# for file in *.gbk; do name=$(basename $file .gbk); mv "$name".* $name; done 


# Prokka (Local genomic annotation)
# Activate Conda environment
conda activate prokka
# Loop through a list of files
for file in 6_assemblies/*.fasta; do
    #Extract file name
    filename=${file##*/}
    #Extract sample name
    prefix=${filename%%.*}
    prokka \
    --cpus $(nproc) \
    --addgenes \
    --centre "" \
    --outdir 10_genome_annotation/${prefix} \
    --prefix ${prefix} \
    $file
done
# Deactivate Conda environment
conda activate base

# Generate annotation summary files
echo -e "Sample\tGene\tCDS\trRNA\ttRNA\ttmRNA" > 10_genome_annotation.tsv
find 10_genome_annotation -type f -name "*.txt" | while read -r file; do
    # Extract file name
    filename=${file##*/}
    # Extract sample name
    sample=${filename%%.*}
    # Convert lines to tab separated table lines
    cds=$(grep -i "^CDS:" "$file" | cut -d':' -f2 | tr -d ' ')
    gene=$(grep -i "^gene:" "$file" | cut -d':' -f2 | tr -d ' ')
    rrna=$(grep -i "^rRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    trna=$(grep -i "^tRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    tmrna=$(grep -i "^tmRNA:" "$file" | cut -d':' -f2 | tr -d ' ')
    # Add line to the output file
    echo -e "${sample}\t${gene}\t${cds}\t${rrna}\t${trna}\t${tmrna}" >> 10_genome_annotation.tsv
done
# Compress the output directory
zip -r 10_genome_annotation.zip 10_genome_annotation
# Delete the output directory
# rm -r 10_genome_annotation


#########################################################################
## 10) Genome annotation - Special features

# PanViTa
# Create an output directory
mkdir 10_panvita
# Go to the output directory
cd 10_panvita
# Create symbolic links to input files
ln -s ../10_genome_annotation/*/*.gbk .
# Activate Conda environment
conda activate panvita
# Run PanViTa
python3 ~/panvita/panvita.py -card -vfdb -bacmet *.gbk
# Leave output directory
cd ..
# Compress the output directory
zip -r 10_panvita.zip 10_panvita
# Delete the output directory
rm -r 10_panvita


# RGI (Antimicrobial resistance prediction)
conda activate rgi
mkdir 10_rgi
# Create file to save the results table
echo -e "Genome\tORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_type\tSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug Class\tResistance Mechanism\tAMR Gene Family\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage Length of Reference Sequence\tID\tModel_ID\tNudged\tNote\tHit_Start\tHit_End\tAntibiotic" > 10_rgi/rgi_all.txt
# Loop through a list of files
for file in 10_genome_annotation/*/*.faa; do
    sample=$(basename ${file} .faa)
    #Run RGI
    rgi main -n $(nproc) --clean -t protein -i $file -o 10_rgi/"$sample"
    sed -i "s/^/$sample\t/" 10_rgi/"$sample".txt #Sample name
    tail -n +2 10_rgi/"$sample".txt >> 10_rgi/rgi_all.txt #Results
    echo >> 10_rgi/rgi_all.txt #Empty line
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 10_rgi.zip 10_rgi
# Delete the output directory
# rm -r 10_rgi


# ResFinder (Antimicrobial resistance prediction)
# http://genepi.food.dtu.dk/resfinder


# CRISPRcasFinder local (CRISPR/Cas system prediction)
mkdir 10_crisprcasfinder
cd 10_crisprcasfinder
# Copy binaries to the current directory
cp -r /mnt/4tb_1/CRISPRCasFinder/* .
# Activate Conda environment
conda activate crisprcasfinder
# Install dependencies
macsydata install -u CASFinder==3.1.0 # You only have to run this command the first time you use this software
# Loop through a list of files
for file in ../6_assemblies/*.fasta; do
    prefix=$(basename ${file} .fasta)
    # Run CRISPR-CasFinder
    perl CRISPRCasFinder.pl -cas -html -in $file -out "${prefix}"
done
# Deactivate Conda environment
conda activate base

# Remove CRISPRCasFinder files
rm -r CasFinder-2.0.3 ccf.environment.yml COPYING COPYRIGHT CRISPRCasFinder.patch CRISPRCasFinder.pl CRISPRCasFinder_Viewer_manual.pdf crispr-cas_logo.png *.sh *.txt install_test README.md sel392v2.so singularity supplementary_files supplementary_tools
# Concatenate CRISPR-CasFinder summaries
> CRISPR-Cas_summary_all.tsv
for directory in */; do
    sample=$(basename $directory /)
    echo -e $sample >> CRISPR-Cas_summary_all.tsv
    cat "$sample"/TSV/CRISPR-Cas_summary.tsv >> CRISPR-Cas_summary_all.tsv
    echo -e "\n" >> CRISPR-Cas_summary_all.tsv
done

# Edit the header
sed -i 's/Evidence-levels/Nb_arrays_evidence-level_1\tNb_arrays_evidence-level_2\tNb_arrays_evidence-level_3\tNb_arrays_evidence-level_4/g' CRISPR-Cas_summary_all.tsv
# Separate into columns
sed -i 's/,Nb/\tNb/g' CRISPR-Cas_summary_all.tsv
# Remove the prefix
sed -i 's/Nb_arrays_evidence-level_.=//g' CRISPR-Cas_summary_all.tsv
# Remove commas
sed -i 's/,\t/\t/g' CRISPR-Cas_summary_all.tsv

# Leave the directory
cd ..

# Copy the summary file
cp 10_crisprcasfinder/CRISPR-Cas_summary_all.tsv 10_crisprcasfinder_summary.tsv
# Compress the output directory
zip -r 10_crisprcasfinder.zip 10_crisprcasfinder
# Delete the output directory
rm -r 10_crisprcasfinder


# CRISPRcasFinder online (CRISPR/Cas system prediction)
# https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index
# Send the chromosome file (*.fasta)


# PHASTEST (Prophage prediction)
# https://phastest.ca/
# Send the GenBank file (*.gb, *.gbk or *.gbff) or the chromosome file (*.fasta)


# PHASTER (Prophage prediction)
# https://phaster.ca/
# Send the GenBank file (*.gb, *.gbk or *.gbff) or chromosome file (*.fasta)


# FastMLST (Multi-Locus Sequence Typing)
# https://github.com/EnzoAndree/FastMLST
conda activate fastmlst
mkdir 10_fastmlst
fastmlst \
-t $(nproc) \
-s '\t' \
-to 10_fastmlst/fastmlst.tsv \
-fo 10_fastmlst/mlst_concat_alleles.fasta \
--splited-output 10_fastmlst/mlst_alleles \
10_genome_annotation/*/*.fna
sed -i s/.fna// 10_fastmlst/fastmlst.tsv
# Compress the output directory
zip -r 10_fastmlst.zip 10_fastmlst
# Delete the output directory
rm -r 10_fastmlst


# PubMLST
# https://pubmlst.org/organisms
# Enter the organism page, enter Typing page
# "Query a sequence" -> "Single sequence"
# Send the chromosome file (*.fasta) and click on "Submit"


# EggNOG-mapper (Gene functional annotation)
# Activate Conda environment
conda activate eggnog-mapper
mkdir -p 10_eggnog_mapper
# Loop through a list of files
for file in 10_genome_annotation/*/*.faa; do
    prefix=$(basename ${file} .faa)
    # Run EggNOG-mapper
    emapper.py --cpu $(nproc) -i $file --output "${prefix}" --output_dir 10_eggnog_mapper
    tail +5 10_eggnog_mapper/"${prefix}".emapper.annotations > 10_eggnog_mapper/"${prefix}".emapper.annotations.tsv
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 10_eggnog_mapper.zip 10_eggnog_mapper
# Delete the output directory
rm -r 10_eggnog_mapper emappertmp*


# COG classifier (Gene functional annotation)
# Create an output directory
mkdir 10_cogclassifier
# Activate Conda environment
conda activate cogclassifier
# Loop through a list of files
for file in 10_genome_annotation/*/*.faa; do
    prefix=$(basename ${file} .faa)
    # Run COG classifier
    COGclassifier -t $(nproc) -i $file -o 10_cogclassifier/"${prefix}"
done
# Deactivate Conda environment
conda activate base
# Compress the output directory
zip -r 10_cogclassifier.zip 10_cogclassifier
# Delete the output directory
rm -r 10_cogclassifier


# 10) Blastkoala (Gene functional annotation)
# https://www.kegg.jp/blastkoala/