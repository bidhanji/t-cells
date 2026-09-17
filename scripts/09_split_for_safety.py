import os

def split_fasta_for_safety(input_file, output_dir, batch_size=100):
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return

    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    batch_count = 1
    current_batch_lines = []
    sequences_in_batch = 0
    
    print(f"Splitting {input_file} into batches of {batch_size} sequences...")
    
    for line in lines:
        if line.startswith('>'):
            # यदि पुरानो batch भरिसकेको छ भने, त्यसलाई save गर र नयाँ सुरु गर
            if sequences_in_batch >= batch_size and current_batch_lines:
                out_file = os.path.join(output_dir, f"safety_batch_{batch_count:03d}.fasta")
                with open(out_file, 'w', encoding='utf-8') as out_f:
                    out_f.writelines(current_batch_lines)
                print(f"  Saved: {out_file} ({sequences_in_batch} sequences)")
                
                batch_count += 1
                current_batch_lines = []
                sequences_in_batch = 0
        
        current_batch_lines.append(line)
        if line.startswith('>'):
            sequences_in_batch += 1
            
    # अन्तिम बाँकी रहेको batch लाई save गर्ने
    if current_batch_lines:
        out_file = os.path.join(output_dir, f"safety_batch_{batch_count:03d}.fasta")
        with open(out_file, 'w', encoding='utf-8') as out_f:
            out_f.writelines(current_batch_lines)
        print(f"  Saved: {out_file} ({sequences_in_batch} sequences)")
        
    print(f"\nSuccess! Total batches created: {batch_count}")
    print(f"Location: {output_dir}")
    print("NEXT ACTION: Upload these batches one by one to AllerTOP and ToxinPred.")

if __name__ == "__main__":
    input_path = "data/processed/unique_strong_binders.fasta"
    output_directory = "data/processed/safety_batches"
    
    split_fasta_for_safety(input_path, output_directory, batch_size=100)