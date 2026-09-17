import os
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def build_final_construct(mhci_peps, mhcii_peps, adjuvant="GQDPCVKNLGTCVQCR"):
    parts = []
    
    # 1. Adjuvant at N-terminus with EAAAK linker
    if adjuvant:
        parts.append(adjuvant)
        parts.append("EAAAK")
        
    # 2. MHC-I epitopes with AAY linker
    if mhci_peps:
        parts.append(mhci_peps[0])
        for p in mhci_peps[1:]:
            parts.extend(["AAY", p])
            
    # 3. MHC-II epitopes with GPGPG linker
    if mhcii_peps:
        if parts:
            parts.append("GPGPG")
        parts.append(mhcii_peps[0])
        for p in mhcii_peps[1:]:
            parts.extend(["GPGPG", p])
            
    # 4. His-tag at C-terminus
    parts.append("HHHHHH")
    
    return "".join(parts)

def get_stats(seq):
    X = ProteinAnalysis(seq)
    return {
        'II': X.instability_index(),
        'GRAVY': X.gravy(),
        'pI': X.isoelectric_point(),
        'MW': X.molecular_weight()
    }

def main():
    print("="*70)
    print("FINAL CONSTRUCT BUILDER (Adjuvant + Unique MHC-I/II + His-Tag)")
    print("="*70)
    
    input_csv = 'data/processed/final_vaccine_candidates.csv'
    if not os.path.exists(input_csv):
        print(f"ERROR: {input_csv} not found.")
        return
        
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} total candidates.")
    print(f"Columns: {list(df.columns)}")
    
    # Separate MHC-I and MHC-II
    df_mhci = df[df['T_Cell_Class'] == 'CD8+ (MHC-I)'].copy()
    df_mhcii = df[df['T_Cell_Class'] == 'CD4+ (MHC-II)'].copy()
    
    print(f"MHC-I candidates: {len(df_mhci)} | MHC-II candidates: {len(df_mhcii)}")
    
    # Find peptide column
    pep_col = [c for c in df.columns if 'peptide' in str(c).lower()]
    pep_col = pep_col[0] if pep_col else 'Peptide'
    
    # Clean NaNs
    df_mhci = df_mhci[df_mhci[pep_col].notna()]
    df_mhcii = df_mhcii[df_mhcii[pep_col].notna()]
    
    # Sort MHC-I by Rank
    mhci_rank_col = [c for c in df_mhci.columns if 'rank' in str(c).lower()]
    if mhci_rank_col:
        df_mhci = df_mhci.sort_values(by=mhci_rank_col[0], ascending=True)
        print(f"Sorted MHC-I by '{mhci_rank_col[0]}'")
        
    # Sort MHC-II by NB (if exists) or Rank
    nb_col = [c for c in df_mhcii.columns if 'nb' in str(c).lower()]
    rank_col = [c for c in df_mhcii.columns if 'rank' in str(c).lower()]
    
    if nb_col:
        df_mhcii = df_mhcii.sort_values(by=nb_col[0], ascending=False)
        print(f"Sorted MHC-II by '{nb_col[0]}' (Descending)")
    elif rank_col:
        df_mhcii = df_mhcii.sort_values(by=rank_col[0], ascending=True)
        print(f"Sorted MHC-II by '{rank_col[0]}' (Ascending)")
    else:
        print("WARNING: Could not find NB or Rank column for MHC-II. Using default order.")
        
    # Calculate GRAVY for top candidates
    df_mhci['GRAVY'] = df_mhci[pep_col].apply(lambda p: ProteinAnalysis(str(p)).gravy())
    df_mhcii['GRAVY'] = df_mhcii[pep_col].apply(lambda p: ProteinAnalysis(str(p)).gravy())
    
    # Sort by GRAVY (most hydrophilic first)
    df_mhci = df_mhci.sort_values(by='GRAVY', ascending=True)
    df_mhcii = df_mhcii.sort_values(by='GRAVY', ascending=True)
    
    # Select exactly 4 UNIQUE epitopes for each class
    top_n = 4
    selected_mhci = df_mhci.head(top_n * 2)[pep_col].drop_duplicates().head(top_n).tolist()
    selected_mhcii = df_mhcii.head(top_n * 2)[pep_col].drop_duplicates().head(top_n).tolist()
    
    print(f"\nSelected {len(selected_mhci)} UNIQUE MHC-I epitopes.")
    print(f"Selected {len(selected_mhcii)} UNIQUE MHC-II epitopes.")
    
    if len(selected_mhcii) == 0:
        print("\nWARNING: No MHC-II epitopes found. Building MHC-I only construct.")
    
    # Build construct
    final_seq = build_final_construct(selected_mhci, selected_mhcii)
    stats = get_stats(final_seq)
    
    # Save outputs
    out_seq = 'data/processed/mev_construct_FINAL.txt'
    with open(out_seq, 'w', encoding='utf-8') as f:
        f.write(final_seq)
        
    out_csv = 'data/processed/final_epitope_list.csv'
    final_df = pd.DataFrame({
        'Peptide': selected_mhci + selected_mhcii,
        'Class': ['MHC-I (CD8+)']*len(selected_mhci) + ['MHC-II (CD4+)']*len(selected_mhcii)
    })
    final_df.to_csv(out_csv, index=False)
    
    print("\n" + "="*70)
    print("FINAL CONSTRUCT GENERATED")
    print("="*70)
    print(f"Sequence Length: {len(final_seq)} aa (Includes Adjuvant + His-Tag)")
    print(f"Instability Index: {stats['II']:.2f} (Target < 40)")
    print(f"GRAVY: {stats['GRAVY']:.3f} (Target < 0)")
    print(f"pI: {stats['pI']:.2f}")
    print(f"Molecular Weight: {stats['MW']:.2f} Da")
    print(f"\nSequence saved to: {out_seq}")
    print(f"Epitope list saved to: {out_csv}")
    print("\nNEXT ACTION: Copy the sequence from mev_construct_FINAL.txt and submit to VaxiJen v2.0 (Target: Virus).")

if __name__ == "__main__":
    main()