import os
import pandas as pd

def select_top_and_build_construct(input_csv, output_csv, output_seq, top_n=8):
    if not os.path.exists(input_csv):
        print(f"ERROR: {input_csv} not found.")
        return

    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} total candidates.")
    
    # Separate MHC-I and MHC-II
    df_mhci = df[df['T_Cell_Class'] == 'CD8+ (MHC-I)'].copy()
    df_mhcii = df[df['T_Cell_Class'] == 'CD4+ (MHC-II)'].copy()
    
    print(f"MHC-I candidates: {len(df_mhci)} | MHC-II candidates: {len(df_mhcii)}")

    # --- Find Peptide Column Dynamically ---
    pep_col = [c for c in df.columns if 'peptide' in str(c).lower()]
    pep_col = pep_col[0] if pep_col else 'Peptide'

    # --- MHC-I Sorting (By Rank, Ascending) ---
    mhci_rank_col = [c for c in df_mhci.columns if 'rank' in str(c).lower()]
    if mhci_rank_col:
        df_mhci = df_mhci.sort_values(by=mhci_rank_col[0], ascending=True)
        print(f"Sorted MHC-I by '{mhci_rank_col[0]}' (Ascending).")
    else:
        print("WARNING: Could not find Rank column for MHC-I.")
        
    top_mhci = df_mhci.head(top_n)
    
    # --- MHC-II Sorting (By NB Descending, then Rank_EL Ascending) ---
    nb_col = [c for c in df_mhcii.columns if 'nb' in str(c).lower()]
    rank_col = [c for c in df_mhcii.columns if 'rank' in str(c).lower()]
    
    if nb_col:
        df_mhcii = df_mhcii.sort_values(by=nb_col[0], ascending=False)
    if rank_col:
        df_mhcii = df_mhcii.sort_values(by=rank_col[0], ascending=True)
        
    top_mhcii = df_mhcii.head(top_n)
    print(f"Sorted MHC-II by NB (Desc) and Rank (Asc).")
    
    # --- Save Shortlist for Manual Verification ---
    selected_df = pd.concat([top_mhci, top_mhcii], ignore_index=True)
    selected_df.to_csv(output_csv, index=False)
    print(f"\n✅ Saved Top {top_n} MHC-I and Top {top_n} MHC-II shortlist to: {output_csv}")
    print("ACTION: Open this CSV in Excel. If you want to swap any epitope, do it now.")
    
    # --- Build MEV Construct Sequence (Bug Fixed) ---
    # Filter out NaN or empty values and ensure everything is a string
    mhci_seqs = [str(s).strip() for s in top_mhci[pep_col].tolist() if pd.notna(s) and str(s).strip()]
    mhcii_seqs = [str(s).strip() for s in top_mhcii[pep_col].tolist() if pd.notna(s) and str(s).strip()]
    
    construct_parts = []
    
    # 1. Add MHC-I epitopes with AAY linker
    if mhci_seqs:
        construct_parts.append(mhci_seqs[0])
        for seq in mhci_seqs[1:]:
            construct_parts.append("AAY")
            construct_parts.append(seq)
            
    # 2. Add MHC-II epitopes with GPGPG linker
    if mhcii_seqs:
        if construct_parts:
            construct_parts.append("GPGPG") # Linker between MHC-I block and MHC-II block
        construct_parts.append(mhcii_seqs[0])
        for seq in mhcii_seqs[1:]:
            construct_parts.append("GPGPG")
            construct_parts.append(seq)
            
    if not construct_parts:
        print("ERROR: No valid sequences found to build the construct. Check your CSV.")
        return
        
    final_sequence = "".join(construct_parts)
    
    # Add His-tag at the end for purification (Standard practice)
    final_sequence += "HHHHHH"
    
    with open(output_seq, 'w', encoding='utf-8') as f:
        f.write(final_sequence)
        
    print(f"\n✅ MEV Construct Sequence saved to: {output_seq}")
    print(f"Total Construct Length: {len(final_sequence)} amino acids (including His-tag)")
    print("NEXT ACTION: Submit this sequence to ProtParam (Expasy) for physicochemical properties.")

if __name__ == "__main__":
    select_top_and_build_construct(
        'data/processed/final_vaccine_candidates.csv',
        'data/processed/top_selected_epitopes.csv',
        'data/processed/mev_construct_sequence.txt',
        top_n=8  # Change this number if you want more or fewer epitopes
    )