import os
import csv
from Bio import SeqIO

def extract_conserved_fasta(csv_path, source_fasta_path, output_fasta_path):
    # Step 1: Read the conserved ortholog IDs from the CSV
    conserved_ids = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            conserved_ids.append(row["LSDV_Ortholog"])
            
    print(f"Found {len(conserved_ids)} conserved protein IDs in CSV.")
    
    # Step 2: Parse the original LSDV FASTA and filter
    extracted_count = 0
    with open(output_fasta_path, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(source_fasta_path, "fasta"):
            if record.id in conserved_ids:
                SeqIO.write(record, out_handle, "fasta")
                extracted_count += 1
                
    print(f"Successfully extracted {extracted_count} sequences to {output_fasta_path}")
    
    if extracted_count != len(conserved_ids):
        print("WARNING: Mismatch in counts. Some IDs from CSV were not found in FASTA.")

if __name__ == "__main__":
    csv_file = "data/processed/conserved_pan_capriproteins.csv"
    source_fasta = "data/processed/LSDV_proteins.fasta"
    output_fasta = "data/processed/conserved_targets.fasta"
    
    if not os.path.exists(csv_file) or not os.path.exists(source_fasta):
        print("ERROR: Required input files not found. Run previous scripts first.")
    else:
        extract_conserved_fasta(csv_file, source_fasta, output_fasta)
        print("Preparation for epitope prediction complete.")