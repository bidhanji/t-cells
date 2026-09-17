import os
import re

def ultimate_fasta_cleaner(input_file, output_dir, batch_size=50):
    if not os.path.exists(input_file):
        print(f"ERROR: File not found: {input_file}")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    clean_sequences = []
    
    current_id = ""
    current_seq = ""
    
    print("Cleaning sequences and removing invalid characters...")
    
    for line in lines:
        line = line.strip() # Removes \r, \n, and spaces
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence if valid
            if current_id and current_seq:
                # Check if sequence contains ONLY valid amino acids
                if all(aa in valid_amino_acids for aa in current_seq):
                    if 8 <= len(current_seq) <= 50: # Standard epitope length
                        clean_sequences.append((current_id, current_seq))
            
            # Create ultra-simple header: >Seq_0001 (No spaces, no special chars)
            # Extract only alphanumeric characters from the old header
            clean_header = re.sub(r'[^a-zA-Z0-9_]', '', line[1:].split()[0])
            if not clean_header:
                clean_header = f"Seq_{len(clean_sequences)+1:04d}"
            current_id = clean_header
            current_seq = ""
        else:
            # Append to sequence, ensure uppercase
            current_seq += line.upper()
            
    # Don't forget the last sequence
    if current_id and current_seq:
        if all(aa in valid_amino_acids for aa in current_seq):
            if 8 <= len(current_seq) <= 50:
                clean_sequences.append((current_id, current_seq))
                
    print(f"Total valid sequences after strict cleaning: {len(clean_sequences)}")
    
    # Split into web-safe batches
    batch_count = 1
    for i in range(0, len(clean_sequences), batch_size):
        batch_seqs = clean_sequences[i:i+batch_size]
        out_file = os.path.join(output_dir, f"web_safe_batch_{batch_count:02d}.fasta")
        
        # newline='\n' forces Unix line endings, preventing \r errors
        with open(out_file, 'w', encoding='utf-8', newline='\n') as f:
            for seq_id, seq in batch_seqs:
                f.write(f">{seq_id}\n")
                f.write(f"{seq}\n")
                
        print(f"  Saved: {out_file} ({len(batch_seqs)} sequences)")
        batch_count += 1
        
    print(f"\nSuccess! {batch_count - 1} 100% web-safe batches created.")
    print(f"Location: {output_dir}")
    print("These files will NOT trigger 'invalid character' errors.")

if __name__ == "__main__":
    # Assuming the previous promiscuous file is the source
    input_path = "data/processed/promiscuous_epitopes.fasta"
    output_directory = "data/processed/ultimate_web_safe_batches"
    
    ultimate_fasta_cleaner(input_path, output_directory, batch_size=50)