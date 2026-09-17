import pandas as pd
import os

def main():
    df = pd.read_csv('data/processed/final_vaccine_candidates.csv')
    
    # Print actual columns for transparency
    print("Columns in final_vaccine_candidates.csv:", list(df.columns))
    
    # Create a unified 'Final_Peptide' column by checking all possible names
    possible_cols = ['Peptide', 'Epitope', 'Sequence', 'peptide', 'epitope', 'sequence']
    df['Final_Peptide'] = None
    for col in possible_cols:
        if col in df.columns:
            df['Final_Peptide'] = df['Final_Peptide'].fillna(df[col])
            
    # Separate
    mhci_df = df[df['T_Cell_Class'] == 'CD8+ (MHC-I)'].copy()
    mhcii_df = df[df['T_Cell_Class'] == 'CD4+ (MHC-II)'].copy()
    
    # Drop NaNs and string 'nan'
    mhci_df = mhci_df[mhci_df['Final_Peptide'].notna() & (mhci_df['Final_Peptide'] != 'nan') & (mhci_df['Final_Peptide'] != 'None')]
    mhcii_df = mhcii_df[mhcii_df['Final_Peptide'].notna() & (mhcii_df['Final_Peptide'] != 'nan') & (mhcii_df['Final_Peptide'] != 'None')]
    
    print(f"Valid MHC-I epitopes found: {len(mhci_df)}")
    print(f"Valid MHC-II epitopes found: {len(mhcii_df)}")
    
    if len(mhci_df) == 0:
        print("ERROR: Still 0 MHC-I epitopes found. Check the column names printed above.")
        return

    # Select top 3 unique MHC-I
    mhci_peps = mhci_df['Final_Peptide'].astype(str).drop_duplicates().head(3).tolist()
    
    # Select top 3 unique MHC-II (with strict 9-mer overlap filtering)
    mhcii_all = mhcii_df['Final_Peptide'].astype(str).drop_duplicates().tolist()
    mhcii_peps = []
    for pep in mhcii_all:
        if len(mhcii_peps) >= 3:
            break
        overlap_found = False
        for selected in mhcii_peps:
            # Check if they share a 9-mer window
            for i in range(len(pep) - 8):
                if pep[i:i+9] in selected:
                    overlap_found = True
                    break
            if overlap_found:
                break
        if not overlap_found:
            mhcii_peps.append(pep)
            
    print(f"\nSelected {len(mhci_peps)} UNIQUE MHC-I epitopes: {mhci_peps}")
    print(f"Selected {len(mhcii_peps)} UNIQUE MHC-II epitopes: {mhcii_peps}")
    
    if len(mhci_peps) < 3 or len(mhcii_peps) < 3:
        print("WARNING: Could not find 3 unique epitopes for both classes. Proceeding with available.")

    # Build final sequence
    adj = "GQDPCVKNLGTCVQCR" # Bovine beta-defensin 3 fragment
    seq_parts = [adj, "EAAAK"] # Adjuvant + rigid linker
    
    # Add MHC-I with AAY linkers
    for i, pep in enumerate(mhci_peps):
        seq_parts.append(str(pep))
        if i < len(mhci_peps) - 1:
            seq_parts.append("AAY")
            
    # Add GPGPG linker between MHC-I and MHC-II blocks
    if mhcii_peps:
        seq_parts.append("GPGPG")
        
    # Add MHC-II with GPGPG linkers
    for i, pep in enumerate(mhcii_peps):
        seq_parts.append(str(pep))
        if i < len(mhcii_peps) - 1:
            seq_parts.append("GPGPG")
            
    # Add His-tag for purification
    seq_parts.append("HHHHHH")
    
    final_seq = "".join(seq_parts)
    
    print("\n" + "="*70)
    print("FINAL PERFECT SEQUENCE (Copy this exact string):")
    print("="*70)
    print(final_seq)
    print(f"\nTotal Length: {len(final_seq)} amino acids")
    
    # Save to file
    out_path = "data/processed/mev_construct_PERFECT.txt"
    with open(out_path, "w") as f:
        f.write(final_seq)
    print(f"Sequence saved to: {out_path}")
    print("\nNEXT ACTION: Copy the sequence above and submit it to VaxiJen v2.0 (Target: Virus, Threshold: 0.4).")

if __name__ == "__main__":
    main()