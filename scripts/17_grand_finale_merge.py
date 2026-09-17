import os
import pandas as pd

def main():
    print("="*70)
    print("GRAND FINALE: MERGING MHC-I AND MHC-II FINAL CANDIDATES")
    print("="*70)
    
    # 1. Load existing MHC-I candidates (already filtered for safety)
    mhci_final_file = "data/processed/final_vaccine_candidates.csv"
    if os.path.exists(mhci_final_file):
        df_mhci = pd.read_csv(mhci_final_file)
        if 'T_Cell_Class' in df_mhci.columns:
            df_mhci = df_mhci[df_mhci['T_Cell_Class'] == 'CD8+ (MHC-I)']
        print(f"Loaded {len(df_mhci)} MHC-I candidates from previous safety filtering.")
    else:
        print(f"WARNING: {mhci_final_file} not found. MHC-I list will be empty.")
        df_mhci = pd.DataFrame()

    # 2. Load MHC-II Allergen-Safe peptides
    mhcii_allergen_safe_file = "data/processed/mhcii_allergen_safe.csv"
    if not os.path.exists(mhcii_allergen_safe_file):
        print(f"ERROR: {mhcii_allergen_safe_file} not found.")
        return
    df_allergen_safe = pd.read_csv(mhcii_allergen_safe_file)
    pep_col_allergen = [c for c in df_allergen_safe.columns if 'peptide' in str(c).lower() or 'clean_seq' in str(c).lower()]
    pep_col_allergen = pep_col_allergen[0] if pep_col_allergen else 'Peptide'
    allergen_safe_seqs = set(df_allergen_safe[pep_col_allergen].dropna().astype(str).str.strip().str.upper())
    print(f"Loaded {len(allergen_safe_seqs)} MHC-II peptides safe from allergenicity.")

    # 3. Load MHC-II ToxinPred results and extract Non-Toxin sequences
    toxinpred_dir = "data/processed/toxinpredmhcii"
    if not os.path.exists(toxinpred_dir):
        print(f"ERROR: {toxinpred_dir} not found.")
        return
    
    non_toxin_seqs = set()
    print(f"Scanning ToxinPred results in {toxinpred_dir}...")
    for file in os.listdir(toxinpred_dir):
        if file.endswith('.csv'):
            df_tox = pd.read_csv(os.path.join(toxinpred_dir, file))
            if 'Sequence' in df_tox.columns and 'Prediction' in df_tox.columns:
                non_tox = df_tox[df_tox['Prediction'] == 'Non-Toxin']['Sequence'].dropna().astype(str).str.strip().str.upper()
                non_toxin_seqs.update(non_tox)
    print(f"Loaded {len(non_toxin_seqs)} MHC-II peptides safe from toxicity.")

    # 4. Find intersection (Final MHC-II safe peptides)
    final_mhcii_seqs = allergen_safe_seqs.intersection(non_toxin_seqs)
    print(f"\nFinal MHC-II peptides passing BOTH safety filters: {len(final_mhcii_seqs)}")

    if not final_mhcii_seqs:
        print("ERROR: No MHC-II peptides passed both filters.")
        return

    # 5. Filter original MHC-II binding data to get full details (Alleles, Rank, etc.)
    mhcii_binding_file = "data/processed/mhcii_promiscuous_prefilt.csv"
    if not os.path.exists(mhcii_binding_file):
        print(f"WARNING: {mhcii_binding_file} not found. Trying mhcii_strong_binders.csv")
        mhcii_binding_file = "data/processed/mhcii_strong_binders.csv"
        
    df_mhcii_full = pd.read_csv(mhcii_binding_file)
    pep_col_full = [c for c in df_mhcii_full.columns if 'peptide' in str(c).lower() or 'clean_seq' in str(c).lower()]
    pep_col_full = pep_col_full[0] if pep_col_full else 'Peptide'
    
    df_mhcii_final = df_mhcii_full[df_mhcii_full[pep_col_full].astype(str).str.strip().str.upper().isin(final_mhcii_seqs)].copy()
    df_mhcii_final['T_Cell_Class'] = 'CD4+ (MHC-II)'
    print(f"Retained {len(df_mhcii_final)} MHC-II binding records for final candidates.")

    # 6. Combine and Save
    final_df = pd.concat([df_mhci, df_mhcii_final], ignore_index=True)
    
    output_file = "data/processed/final_vaccine_candidates.csv"
    final_df.to_csv(output_file, index=False, encoding='utf-8')
    
    print("\n" + "="*70)
    print("PIPELINE COMPLETE!")
    print("="*70)
    print(f"Total Final Vaccine Candidates: {len(final_df)}")
    print(f"  - MHC-I (CD8+): {len(df_mhci)}")
    print(f"  - MHC-II (CD4+): {len(df_mhcii_final)}")
    print(f"\nUltimate publication-ready file saved to: {output_file}")
    print("You can now proceed to Vaccine Construct Design (Linkers, 3D Modeling, Docking).")

if __name__ == "__main__":
    main()