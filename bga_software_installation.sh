# Bash script for installation of software used for bacterial genome assembly from short-read sequencing data and genome annotation
#
# ⚠️ DO NOT execute this script entirely at once!
# Copy and paste individual command lines into the Linux terminal as needed.
# This file uses the .sh extension only to enable Bash syntax highlighting in text editors.
#
# Author: Marcus Vinicius Canário Viana
# Date: 17/10/2025
# More info: see README.md in the repository


# Summary
# A) System requirements
# B) Software installation
# B.1) Software installation - Miniconda
# B.2) Software installation - Genome assembly
# B.3) Software installation - Genome annotation
# C) Connecting to a server and using Screen
# D) Sending all results to another computer


##########################################################################
## A) System requirements
##########################################################################
# RAM: 16GB (64GB if you run taxonomic analysis with GTDB-Tk)
# Storage: 25GB for Miniconda and its environments, 170GB for databases (142GB of which is from GTDB-Tk)
# Static IP address (For use on a computer server)


##########################################################################
## B) Software installation
##########################################################################

##########################################################################
## B.1) Software installation - Miniconda
##########################################################################

##########################################################################
# Miniconda (For a single user. Local computer.)
# Download the installer
cd
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x87_64.sh
sh Miniconda3-latest-Linux-x87_64.sh
# Follow the instructions to install Miniconda in /home/$USER/miniconda3

# Do you accept the license terms? [yes|no]
# >>> yes

# Miniconda3 will now be installed into this location:
# /root/miniconda3
#
#   - Press ENTER to confirm the location
#   - Press CTRL-C to abort the installation
#   - Or specify a different location below
#
# [/home/user_name/miniconda3] >>>

# Do you wish to update your shell profile to automatically initialize conda?
# This will activate conda on startup and change the command prompt when activated.
# If you'd prefer that conda's base environment not be activated on startup,
#    run the following command when conda is activated:

# conda config --set auto_activate_base false

# You can undo this by running `conda init --reverse $SHELL`? [yes|no]
# [no] >>> yes

# Apply changes in current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before the user name
# Add Conda channels
conda config --show channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
conda config --add channels anaconda
conda config --add channels r
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# Set libmamba as enviroment solver. It is faster than the classic solver
conda config --set solver libmamba
# Update Miniconda
conda update conda -y # If the message "Terms of Service have not been accepted" appears, execute the command lines shown in the message and try this command again.
# Show Conda configuration
conda config --show-sources

# Directory for databases
# Create a directory
sudo mkdir /db
# Change the directory owner and group to the current user. This is required to download databases later.
sudo chown $USER /db
sudo chgrp $USER /db

##########################################################################
# Miniconda (For all users. Server computer.)
# Log in as root
sudo su
cd
# Download the installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x87_64.sh
sh Miniconda3-latest-Linux-x87_64.sh
# Follow the instructions and install Miniconda in /usr/local/miniconda3

# Do you accept the license terms? [yes|no]
# >>> yes

# Miniconda3 will now be installed into this location:
# /root/miniconda3
#
#   - Press ENTER to confirm the location
#   - Press CTRL-C to abort the installation
#   - Or specify a different location below
#
# [/root/miniconda3] >>> /usr/local/miniconda3

# Do you wish to update your shell profile to automatically initialize conda?
# This will activate conda on startup and change the command prompt when activated.
# If you'd prefer that conda's base environment not be activated on startup,
#    run the following command when conda is activated:

# conda config --set auto_activate_base false

# You can undo this by running `conda init --reverse $SHELL`? [yes|no]
# [no] >>> yes

# Apply changes in the current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before "root¨
# Add Conda channels
conda config --show channels
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels defaults
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
# Set libmamba as enviroment solver. It is faster than the classic solver
conda config --set solver libmamba
# Update Miniconda
conda update conda -y
# Show Conda configuration
conda config --show-sources
# Create a simbolic link for all users
ln -s /usr/local/miniconda3/bin/conda /usr/bin/conda
# Log out as root
exit

# For each user, execute this command only once to add Conda.
conda init # To write Conda information to ~/.bashrc
# Apply changes in the current session
source ~/.bashrc
# Now, "(base)" should appear in the terminal before the username.

# Directory for databases
# Create a directory
sudo mkdir /db

# IMPORTANT! Log in as root before creating a Conda environment to make it available to all users.


##########################################################################
## B.2) Software installation - Genome assembly
##########################################################################

##########################################################################
# barrnap (Evaluate the completeness of rRNA genes)
conda create -n barrnap -c bioconda barrnap -y

##########################################################################
# CheckM2 (Evaluate the assembly completeness and contamination)
mkdir /db/checkm2
cd
git clone --recursive https://github.com/chklovski/checkm2.git
cd checkm2
conda env create -n checkm2 -f checkm2.yml
conda activate checkm2
python setup.py install
checkm2 -h
checkm2 database --download --path /db/checkm2
conda env config vars set CHECKM2DB="$(find /db/checkm2 -type f -name "*.dmnd" -print -quit)"
conda deactivate
conda activate checkm2
echo $CHECKM2DB
conda deactivate

##########################################################################
# fastq-dl (Download reads from ENA or GenBank SRA)
conda create -n fastq-dl -c bioconda fastq-dl -y

##########################################################################
# FastQC (Evaluate read quality)
conda create -n fastqc -c bioconda fastqc -y

##########################################################################
# Fastp (Trim sequencing reads)
conda create -n fastp -c bioconda fastp -y

############################################################
# Fastplong
conda create -n fastplong -c bioconda fastplong -y

############################################################
# Flye
conda create -n flye -c bioconda flye -y

############################################################
# GenomeScope (Estimation of genome size)
conda create -n genomescope -c bioconda genomescope2 -y

##########################################################################
# Git (Clone software directories from github.com)
sudo apt-get update -y
sudo apt-get install git -y

##########################################################################
# Gzip (Compress and decompress directories and files)
sudo apt-get install gzip -y

##########################################################################
# GTDB-Tk (Taxonomyc analysis)
mkdir -p /db/gtdbtk/
cd /db/
wget -c --progress=bar https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
# wget -c --progress=bar https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
# echo '24b476ea5a4ef30519d461e56cc4a27f gtdbtk_data.tar.gz' | md5sum -c #md5sum for version r226
tar -xvzf gtdbtk_data.tar.gz -C "/db/gtdbtk" --strip 1 > /dev/null
conda create -n gtdbtk -c bioconda gtdbtk=2.4.1 -y
conda activate gtdbtk
conda env config vars set GTDBTK_DATA_PATH="/db/gtdbtk"
conda deactivate
conda activate gtdbtk
echo $GTDBTK_DATA_PATH
gtdbtk check_install
conda deactivate

##########################################################################
# GUNC (Evaluate the assembly contamination)
mkdir /db/gunc
conda create -n gunc -c bioconda gunc -y
conda activate gunc
conda install setuptools -y
gunc download_db /db/gunc/
conda env config vars set GUNC_DB=$(find /db/gunc/ -type f -name "*.dmnd" -print -quit)
conda deactivate
conda activate gunc
echo $GUNC_DB
conda deactivate

############################################################
## KMC (K-mer counting of sequencing reads)
conda create -n kmc -c bioconda kmc -y

##########################################################################
# MultiQC (Unify read quality report)
conda create -n multiqc -c bioconda multiqc -y

############################################################
## NanoPlot
conda create -n nanoplot -c bioconda nanoplot -y

##########################################################################
# QUAST (Evaluate assembly fragmentation)
conda create -n quast -c bioconda quast -y

##########################################################################
# MOB-suite (Identify genetic mobile elements)
conda create -n mob_suite -c bioconda mob_suite -y

##########################################################################
# Openssh (Access computer remotely)
sudo apt install openssh-server openssh-client
sudo ufw allow ssh
sudo ufw enable
sudo ufw status
sudo echo "ClientAliveInterval 60" >> /etc/ssh/sshd_config
sudo systemctl restart ssh.service

############################################################
# Rasusa
conda create -n rasusa -c bioconda rasusa -y

##########################################################################
# Rename (Rename files using regular expressions)
sudo apt-get install rename -y

##########################################################################
# Screen (Run the program in different screens)
sudo apt-get install screen -y

##########################################################################
# seqkit (Filter sequences by size)
conda create -n seqkit -c bioconda seqkit -y

##########################################################################
# Shovill (De novo assembly)
conda create -n spades -c bioconda shovill -y

##########################################################################
# SPAdes (De novo assembly)
conda create -n spades -c bioconda spades -y

##########################################################################
# Unicycler (De novo assembly)
conda create -n unicycler -c bioconda unicycler -y

##########################################################################
# Zip (Compress and decompress directories and files)
sudo apt-get install zip -y


##########################################################################
## B.3) Software installation - Genome annotation
##########################################################################

##########################################################################
# COG classifier (Gene functional annotation)
conda create -n cogclassifier -c bioconda -c conda-forge cogclassifier -y

##########################################################################
# CRISPRcasFinder local (CRISPR/Cas system prediction)
cd /db
git clone https://github.com/dcouvin/CRISPRCasFinder.git
cd CRISPRCasFinder
conda env create -f ccf.environment.yml -n crisprcasfinder
conda activate crisprcasfinder
mamba install -c bioconda macsyfinder=2.1.2 -y
macsydata install -u CASFinder==3.1.0
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
# FastMLST (Multi-Locus Sequence Typing)
conda create -n fastmlst -c bioconda fastmlst -y
conda activate fastmlst
fastmlst -t 1 --update-mlst
conda deactivate

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
# Prokka (Genome annotation)
conda create -n prokka -c bioconda prokka -y

##########################################################################
# RGI (Antimicrobial resistance prediction)
conda create -n rgi -c conda-forge -c bioconda -c defaults rgi -y

##########################################################################
## C) Connecting to a server and using Screen
##########################################################################

# Connection to server
ssh user_name@server_ip

# Create screen
screen -S screen_name

# Leave screen without closing it (Detach)
CTRL+a+d

# List screens
screen -ls

# Retrive a screen (Attach)
screen -rd screen_name

# Close a screen (When attached to it)
exit

# Kill a screen (When a process inside the screen freezes)
screen -ls #To list screens
screen -XS screen_name kill


##########################################################################
## D) Sending all results to another computer
##########################################################################

# Create a checksum file
md5sum 3_fastp/*.gz *.zip *.tsv > assembly_pipeline.md5
# Send trimmed reads, results, and checksum files to the destination computer
scp 3_fastp/*.gz *.zip *.tsv assembly_pipeline.md5 user@computerip:/destination/directory
# Check the files on the destination computer
md5sum -c assembly_pipeline.md5
