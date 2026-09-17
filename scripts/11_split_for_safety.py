import os

def split_fasta_safely(input_file, output_dir, batch_size=100):
    if not os.path.exists(input_file):
        print(f"ERROR: File not found: {input_file}")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    batch_count = 1
    current_batch_lines = []
    sequences_in_batch = 0
    
    print(f"Splitting {len(lines)//2} sequences into batches of {batch_size}...")
    
    for line in lines:
        if line.startswith('>'):
            if sequences_in_batch >= batch_size and current_batch_lines:
                out_file = os.path.join(output_dir, f"promiscuous_batch_{batch_count:02d}.fasta")
                with open(out_file, 'w', encoding='utf-8') as out_f:
                    out_f.writelines(current_batch_lines)
                print(f"  Saved: {out_file} ({sequences_in_batch} sequences)")
                
                batch_count += 1
                current_batch_lines = []
                sequences_in_batch = 0
        
        current_batch_lines.append(line)
        if line.startswith('>'):
            sequences_in_batch += 1
            
    if current_batch_lines:
        out_file = os.path.join(output_dir, f"promiscuous_batch_{batch_count:02d}.fasta")
        with open(out_file, 'w', encoding='utf-8') as out_f:
            out_f.writelines(current_batch_lines)
        print(f"  Saved: {out_file} ({sequences_in_batch} sequences)")
        
    print(f"\nSuccess! Total batches created: {batch_count}")
    print(f"Location: {output_dir}")

if __name__ == "__main__":
    input_path = "data/processed/promiscuous_epitopes.fasta"
    output_directory = "data/processed/safety_batches_final"
    
    split_fasta_safely(input_path, output_directory, batch_size=100)