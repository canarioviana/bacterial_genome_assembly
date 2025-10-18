# Bash Scripts for Bacterial Genome Assembly

This repository contains a collection of command-lines and **Bash scripts** essential for bacterial genome **assembly from short or long reads** and **genome annotation**.

It also includes detailed instructions for the installation of all necessary software.

---

## 📋 Usage Guidelines

* The scripts provide comprehensive **software installation instructions**, **a worfkflow for genome assembly with guidelines for genome and sequencing reads submission** and **a worfkflow for genome annotation**

* **_end2end.sh Files (Workflow):** These scripts encapsulate the complete workflow and **can be executed at once** (end-to-end).

    To execute them:
    1.  **Grant execution permission:**
        ```bash
        chmod +x gba_XXX_end2end.sh
        ```
    2.  **Run the script:**
        ```bash
        ./gba_XXX_end2end.sh
        ```

* **_script.sh Files (Modular Commands):** These files are collections of commands grouped by function (e.g., QC only, assembly only).

    **⚠️ IMPORTANT: These scripts SHOULD NOT be executed in their entirety.**

    Instead, you should copy (or modify) and paste the relevant command lines directly into your Linux terminal as needed for modular use.

---

## 🛠️ The Short-Reads Genome Assembly Workflow

### 1) Sequencing reads directory and files
* Reads stored as local files
* Reads from ENA or GenBank
### 2) Raw reads quality assessment
* FastQC
* MultiQC
### 3) Raw reads trimming and downsampling 
* Fastp
* Downsampling (KMC, GenomeScope and Rasusa)
### 4) Trimmed reads quality assessment
* FastQC
* MultiQC
### 5) De novo assembly
* Unicycler
* Shovill
* SPAdes
### 6) Organizing de novo assembly files
### 7) Assembly quality assessment
* CheckM2
* GUNC
* QUAST
* Barrnap
* Calculation of vertical sequencing coverage
### 8) Taxonomic assignment
* GTDB-Tk
* TYGS (online)
### 9) Plasmids identification
* MOB-suite 
### 10) Assignment of contigs to molecules
* MOB-suite and an inhouse script

## 🛠️ The Long-Reads Genome Assembly Workflow

### 1) Sequencing reads directory and files
* Reads stored as local files
### 2) Raw reads quality assessment
* NanoPlot
* FastQC
* MultiQC
### 3) Raw reads trimming and downsampling 
* Fastplong
### 4) Trimmed reads quality assessment
* NanoPlot
* FastQC
* MultiQC
### 5) De novo assembly
* Flye
* Unicycler
* Raven
### 6) Organizing de novo assembly files
### 7) Assembly quality assessment
* CheckM2
* GUNC
* QUAST
* Barrnap
* Calculation of vertical sequencing coverage
### 8) Taxonomic assignment
* GTDB-Tk
* TYGS (online)
### 9) Plasmids identification
* MOB-suite 
### 10) Assignment of contigs to molecules
* MOB-suite and an inhouse script
