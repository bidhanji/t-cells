import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def parse_mhci_results(input_dir, output_file):
    # अब .csv फाइलहरू खोज्नेछौं
    batch_files = sorted([f for f in os.listdir(input_dir) 
                          if f.startswith('iedb_mhci_batch_') and f.endswith('.csv')])
    
    if not batch_files:
        print(f"ERROR: No IEDB MHC-I batch CSV files found in {input_dir}")
        return
    
    print(f"Found {len(batch_files)} MHC-I batch files to parse...")
    
    all_epitopes = []
    
    for batch_file in batch_files:
        file_path = os.path.join(input_dir, batch_file)
        print(f"  Processing {batch_file}...")
        
        try:
            # CSV फाइल पढ्ने
            df = pd.read_csv(file_path)
            
            # Percentile column पत्ता लगाउने
            percentile_col = None
            for col in df.columns:
                if 'percentile' in str(col).lower() and 'el' in str(col).lower():
                    percentile_col = col
                    break
            if not percentile_col:
                for col in df.columns:
                    if 'percentile' in str(col).lower():
                        percentile_col = col
                        break
            
            if not percentile_col:
                print(f"    WARNING: No percentile column found in {batch_file}.")
                continue
            
            # Strong binders (<= 0.5%) मात्र छुट्याउने
            strong_binders = df[df[percentile_col] <= 0.5].copy()
            strong_binders['source_batch'] = batch_file
            
            all_epitopes.append(strong_binders)
            print(f"    -> Found {len(strong_binders)} strong binders (≤0.5%) out of {len(df)} total")
            
        except Exception as e:
            print(f"    ERROR processing {batch_file}: {e}")
            
    if not all_epitopes:
        print("ERROR: No strong binders found.")
        return
        
    combined_df = pd.concat(all_epitopes, ignore_index=True)
    
    # Duplicate हटाउने
    if 'peptide' in combined_df.columns and 'allele' in combined_df.columns:
        combined_df = combined_df.drop_duplicates(subset=['peptide', 'allele'], keep='first')
        
    print(f"\nTotal unique strong-binding MHC-I epitopes: {len(combined_df)}")
    combined_df.to_csv(output_file, index=False, encoding='utf-8')
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    parse_mhci_results("data/processed/batches", "data/processed/final_mhci_epitopes.csv")