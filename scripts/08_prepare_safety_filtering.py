import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def prepare_unique_peptides(mhci_file, mhcii_file, output_fasta):
    all_peptides = set()
    
    # 1. MHC-I बाट unique peptides निकाल्ने
    if os.path.exists(mhci_file):
        df_mhci = pd.read_csv(mhci_file)
        # Case-insensitive search for 'peptide' column
        peptide_cols = [col for col in df_mhci.columns if 'peptide' in col.lower()]
        if peptide_cols:
            actual_col = peptide_cols[0]
            print(f"MHC-I file found. Using column: '{actual_col}'")
            unique_mhci = df_mhci[actual_col].dropna().astype(str).str.strip().unique()
            all_peptides.update(unique_mhci)
            print(f"Extracted {len(unique_mhci)} unique MHC-I peptides.")
        else:
            print(f"ERROR: 'peptide' column not found in MHC-I file. Available columns: {list(df_mhci.columns)}")
            return
    else:
        print(f"ERROR: MHC-I file not found at {mhci_file}")
        return

    # 2. MHC-II बाट unique peptides निकाल्ने
    if os.path.exists(mhcii_file):
        df_mhcii = pd.read_csv(mhcii_file)
        peptide_cols = [col for col in df_mhcii.columns if 'peptide' in col.lower()]
        if peptide_cols:
            actual_col = peptide_cols[0]
            print(f"MHC-II file found. Using column: '{actual_col}'")
            unique_mhcii = df_mhcii[actual_col].dropna().astype(str).str.strip().unique()
            all_peptides.update(unique_mhcii)
            print(f"Extracted {len(unique_mhcii)} unique MHC-II peptides.")
        else:
            print(f"ERROR: 'peptide' column not found in MHC-II file. Available columns: {list(df_mhcii.columns)}")
            return
    else:
        print(f"ERROR: MHC-II file not found at {mhcii_file}")
        return

    if not all_peptides:
        print("ERROR: No peptides found to process.")
        return

    # 3. Duplicate हटाएर FASTA file बनाउने
    unique_list = sorted(list(all_peptides))
    print(f"\nTotal unique strong-binding peptides (MHC-I + MHC-II combined): {len(unique_list)}")
    print("Writing to FASTA file...")
    
    with open(output_fasta, 'w', encoding='utf-8') as f:
        for i, peptide in enumerate(unique_list, 1):
            # FASTA header मा ID र length राख्ने
            f.write(f">Peptide_{i}_Len{len(peptide)}\n")
            f.write(f"{peptide}\n")
            
    print(f"Success! Unique peptides saved to: {output_fasta}")
    print("\nFirst 5 peptides in the file:")
    for p in unique_list[:5]:
        print(f"  - {p}")

if __name__ == "__main__":
    mhci_path = "data/processed/final_mhci_epitopes.csv"
    mhcii_path = "data/processed/mhcii_strong_binders.csv"
    output_path = "data/processed/unique_strong_binders.fasta"
    
    prepare_unique_peptides(mhci_path, mhcii_path, output_path)