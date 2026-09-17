import pandas as pd

def main():
    df = pd.read_csv('data/processed/final_vaccine_candidates.csv')
    
    # Extract top 3 unique MHC I
    mhci = df[df['T_Cell_Class'] == 'CD8+ (MHC-I)']['Peptide'].drop_duplicates().head(3).tolist()
    
    # Extract top 3 unique MHC II
    mhcii = df[df['T_Cell_Class'] == 'CD4+ (MHC-II)']['Peptide'].drop_duplicates().head(3).tolist()

    print(f"Found {len(mhci)} MHC I and {len(mhcii)} MHC II epitopes.")

    if len(mhcii) == 0:
        print("ERROR: No MHC II epitopes found in the CSV.")
        return

    adj = "GQDPCVKNLGTCVQCR"
    seq_parts = [adj, "EAAAK"]

    # Add MHC I with AAY linkers
    for i, pep in enumerate(mhci):
        seq_parts.append(str(pep))
        if i < len(mhci) - 1:
            seq_parts.append("AAY")
            
    # Add GPGPG between MHC I block and MHC II block
    if mhcii:
        seq_parts.append("GPGPG")
        
    # Add MHC II with GPGPG linkers
    for i, pep in enumerate(mhcii):
        seq_parts.append(str(pep))
        if i < len(mhcii) - 1:
            seq_parts.append("GPGPG")
            
    # Add His tag at the end
    seq_parts.append("HHHHHH")

    final_seq = "".join(seq_parts)
    
    print("\nFINAL PERFECT SEQUENCE (Copy this):")
    print(final_seq)
    print(f"\nTotal Length: {len(final_seq)} amino acids")
    print("ACTION: Copy the sequence above and submit it to VaxiJen v2.0 (Target: Virus).")

if __name__ == "__main__":
    main()