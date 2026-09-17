import os
from Bio import SeqIO

def split_fasta(input_fasta, output_dir, batch_size=20):
    os.makedirs(output_dir, exist_ok=True)
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    total = len(records)
    
    print(f"Splitting {total} sequences into batches of {batch_size}...")
    
    for i in range(0, total, batch_size):
        batch_records = records[i:i + batch_size]
        batch_num = (i // batch_size) + 1
        output_file = os.path.join(output_dir, f"batch_{batch_num:02d}.fasta")
        
        with open(output_file, "w", encoding="utf-8") as out_handle:
            SeqIO.write(batch_records, out_handle, "fasta")
            
        print(f"  -> Created {output_file} with {len(batch_records)} sequences")
        
    print("Splitting complete. Ready for web server submission.")

if __name__ == "__main__":
    split_fasta("data/processed/conserved_targets.fasta", "data/processed/batches", batch_size=20)