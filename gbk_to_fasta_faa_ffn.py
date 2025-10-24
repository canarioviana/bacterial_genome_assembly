# Python script to convert GenBank (.gbk or *.gbff) files into Nucleotide (.fasta), 
# Protein (.faa), and Gene Nucleotide (.ffn) sequence files.
#
# Date: 24/10/2025
# More info: see README.md in the repository
#
# Instructions:
#
# Put this script and the .gbk files in the same directory
# Execute the script using the following command
# python gbk_to_fasta_faa_ffn.py

import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Function to format sequences with 60 characters per line
def format_sequence(sequence, line_length=60):
    return "\n".join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))

# Directory where the script is located
gbk_directory = os.path.dirname(os.path.abspath(__file__))
processed_files = []

# Veryfing and renaming *.gb and *.gbff files to *.gbk
extensions_to_rename = [".gb", ".gbff"]
print("Verifying and renaming *.gb and *.gbff files to *.gbk...") # CORREÇÃO: Typo ("Veryfing" -> "Verifying")
for filename in os.listdir(gbk_directory):
    for ext in extensions_to_rename:
        if filename.lower().endswith(ext):
            old_path = os.path.join(gbk_directory, filename)
            
            root = filename[:-len(ext)]
            new_filename = root + ".gbk"
            new_path = os.path.join(gbk_directory, new_filename)
            
            if old_path != new_path:
                try:
                    os.rename(old_path, new_path)
                    print(f"Renamed: {filename} -> {new_filename}")
                except Exception as e:
                    print(f"Error during renaming {filename}: {e}")
            break

# Iterates through all files in the specified directory
for filename in os.listdir(gbk_directory):
    if filename.endswith(".gbk"):
        gbk_filename = os.path.join(gbk_directory, filename)
        
        # CORREÇÃO: Aplica os.path.splitext() a 'filename' (que é só o nome)
        root_name = os.path.splitext(filename)[0] 
        
        # CORREÇÃO: Adiciona o nome raiz à lista de arquivos processados
        processed_files.append(root_name)
        
        # AQUI os paths dos arquivos de saída devem ser corrigidos para usar apenas o nome
        # pois o os.path.splitext(gbk_filename)[0] original incluía o diretório.
        # Agora, os arquivos de saída serão criados no diretório principal.
        faa_filename = os.path.join(gbk_directory, root_name + ".faa")
        fasta_filename = os.path.join(gbk_directory, root_name + ".fasta")
        ffn_filename = os.path.join(gbk_directory, root_name + ".ffn")

        # Open output files
        print(f"Processing the file {filename}") # CORREÇÃO: Typo ("Processind" -> "Processing")
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
                            nucleotide_seq = seq_feature.location.extract(seq_record.seq)
                            
                            faa_output.write(f">{locus_tag} {product}\n{format_sequence(translation)}\n")
                            ffn_output.write(f">{locus_tag} {product}\n{format_sequence(str(nucleotide_seq))}\n")


print("\nOrganizing files into directories...")
# Use set to process each unique root_name only once
for root_name in set(processed_files): 
    
    # Create the target directory
    target_dir = os.path.join(gbk_directory, root_name)
    try:
        # Create the directory if it doesn't exist
        os.makedirs(target_dir, exist_ok=True) 
        print(f"Directory created: {root_name}")
    except Exception as e:
        print(f"Error creating directory {root_name}: {e}")
        continue # Skip moving if the directory could not be created

    # Move the .gbk file and generated files (.faa, .fasta, .ffn)
    files_to_move = [f"{root_name}.gbk", f"{root_name}.faa", f"{root_name}.fasta", f"{root_name}.ffn"]
    
    for filename in files_to_move:
        source_path = os.path.join(gbk_directory, filename)
        destination_path = os.path.join(target_dir, filename)
        
        # Check if the file exists before attempting to move
        if os.path.exists(source_path):
            try:
                os.rename(source_path, destination_path)
            except Exception as e:
                print(f"Error moving file {filename} to {root_name}: {e}")


print("\nProcessing and organization completed successfully.")