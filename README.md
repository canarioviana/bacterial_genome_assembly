# Bash Scripts for Bacterial Genome Assembly

This repository contains a collection of command-lines and **Bash scripts** essential for bacterial genome **assembly from short paired-end reads (Illumina) or single-end long-reads (PacBio/ONT)** and **genome annotation**.

It also includes detailed instructions for the installation of all necessary software.

---

## Usage Guidelines

* The scripts provide comprehensive **software installation instructions**, **a worfkflow for genome assembly with guidelines for genome and sequencing reads submission** and **a worfkflow for genome annotation**

* **_end2end.sh Files (Workflow):** These scripts encapsulate the complete workflow and **can be executed at once** (end-to-end). To execute them:
    
    * **Grant execution permission:**
        ```bash
        chmod +x gba_XXX_end2end.sh
        ```
    * **Run the script:**
        ```bash
        ./gba_XXX_end2end.sh
        ```

* **_script.sh Files (Modular Commands):** These files are collections of commands grouped by function (e.g., QC only, assembly only).

    * **⚠️ IMPORTANT: These scripts SHOULD NOT be executed in their entirety.**

    * Instead, you should copy (or modify) and paste the relevant command lines directly into your Linux terminal as needed for modular use.

---

## The Short Paired-End Reads Genome Assembly Workflow

1) Sequencing reads directory and files
    * Download reads from NCBI SRA (SRA Tools)
    * Check local read files
2) Raw reads quality assessment
    * FastQC
    * MultiQC
3) Raw reads trimming, estimation of genome size and downsampling
    * Fastp
    * Estimation of genome size (KMC and GenomeScope) and downsampling (Rasusa)
4) Trimmed reads quality assessment
    * FastQC
    * MultiQC
5) De novo assembly
    * Unicycler (end-to-end workflow)
    * Shovill
    * SPAdes
6) Organizing de novo assembly files
7) Assembly quality assessment
    * CheckM2
    * GUNC
    * QUAST
    * Barrnap
    * Calculation of vertical sequencing coverage
8) Taxonomic assignment
    * GTDB-Tk
    * TYGS (online)
9) Plasmids identification
    * MOB-suite 
10) Assignment of contigs to molecules
    * MOB-suite and an in-house script

---

## The Long-Reads Genome Assembly Workflow

1) Sequencing reads directory and files
    * Download reads from NCBI SRA (SRA Tools)
    * Check local read files
2) Raw reads quality assessment
    * NanoPlot
3) Raw reads trimming and estimation of genome size
    * Chopper
    * Estimation of genome size (KMC and GenomeScope)
4) Trimmed reads quality assessment
    * NanoPlot
5) De novo assembly
    * Flye
6) Organization of de novo assembly files
7) Assembly quality assessment
    * CheckM2
    * GUNC
    * QUAST
    * Barrnap
    * Calculation of vertical sequencing coverage
8) Taxonomic assignment
    * GTDB-Tk
    * TYGS (online)
9) Plasmids identification
    * MOB-suite 
10) Assignment of contigs to molecules
    * MOB-suite and an in-house script







