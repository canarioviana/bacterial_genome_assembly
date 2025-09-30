# Python script to convert GenBank (.gbk) files into FASTA (.fasta), 
# Protein (.faa), and Gene Nucleotide (.ffn) sequence files.
#
# Author: Marcus Vinicius Canário Viana
# Date: 30/09/2025
# More info: see README.md in the repository

import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Function to format sequences with 60 characters per line
def format_sequence(sequence, line_length=60):
    return "\n".join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

# Directory where the script is located
gbk_directory = os.path.dirname(os.path.abspath(__file__))

# Iterates through all files in the specified directory
for filename in os.listdir(gbk_directory):
    if filename.endswith(".gbk"):
        gbk_filename = os.path.join(gbk_directory, filename)
        root_name = os.path.splitext(gbk_filename)[0]
        faa_filename = root_name + ".faa"
        fasta_filename = root_name + ".fasta"
        ffn_filename = root_name + ".ffn"

        # Open output files
        print(f"Processando o arquivo {filename}")
        with open(faa_filename, "w") as faa_output, open(fasta_filename, "w") as fasta_output, open(ffn_filename, "w") as ffn_output:
            for seq_record in SeqIO.parse(gbk_filename, "genbank"):
                
                # Save contigs in FASTA format
                fasta_output.write(f">{seq_record.id}\n{format_sequence(str(seq_record.seq))}\n")
                
                # Extract proteins and nucleotide sequences from CDS
                for seq_feature in seq_record.features:
                    if seq_feature.type == "CDS" and "translation" in seq_feature.qualifiers:
                        if "pseudo" not in seq_feature.qualifiers:
                            locus_tag = seq_feature.qualifiers.get("locus_tag", ["unknown"])[0]
                            product = seq_feature.qualifiers.get("product", ["unknown"])[0]
                            translation = seq_feature.qualifiers["translation"][0]
                            nucleotide_seq = seq_feature.location.extract(seq_record.seq)  # Obtendo a sequência de nucleotídeos
                            
                            faa_output.write(f">{locus_tag} {product}\n{format_sequence(translation)}\n")
                            ffn_output.write(f">{locus_tag} {product}\n{format_sequence(str(nucleotide_seq))}\n")
