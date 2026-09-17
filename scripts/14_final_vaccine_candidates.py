import os
import re
import pandas as pd

def main():
    fasta_path = 'data/processed/promiscuous_epitopes.fasta'
    if not os.path.exists(fasta_path):
        print("ERROR: promiscuous_epitopes.fasta not found.")
        return

    # 1. Process AllerCatPro results
    allergen_dir = 'data/processed/allertop'
    safe_headers = set()
    
    if os.path.exists(allergen_dir):
        for file in os.listdir(allergen_dir):
            if file.endswith('.csv'):
                df = pd.read_csv(os.path.join(allergen_dir, file))
                # AllerCatPro uses 'Protein' for headers and 'Result' for status
                if 'Result' in df.columns and 'Protein' in df.columns:
                    safe = df[df['Result'] == 'no evidence']['Protein'].tolist()
                    safe_headers.update(safe)
        print(f"Safe headers from AllerCatPro: {len(safe_headers)}")
    else:
        print("WARNING: allertop folder not found.")

    # Map safe headers back to actual sequences from the original FASTA
    safe_from_allergen = set()
    with open(fasta_path, 'r', encoding='utf-8') as f:
        current_clean_header = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                raw_header = line[1:].split()[0]
                # Clean header exactly as AllerCatPro does (remove colons, spaces, etc.)
                current_clean_header = re.sub(r'[^a-zA-Z0-9_]', '', raw_header)
            elif current_clean_header:
                seq = line.upper()
                if current_clean_header in safe_headers and 8 <= len(seq) <= 50:
                    safe_from_allergen.add(seq)
    print(f"Safe from AllerCatPro (Non-allergenic): {len(safe_from_allergen)}")

    # 2. Process ToxinPred results
    toxin_dir = 'data/processed/toxinpred'
    safe_from_toxin = set()
    
    if os.path.exists(toxin_dir):
        for file in os.listdir(toxin_dir):
            if file.endswith('.csv'):
                df = pd.read_csv(os.path.join(toxin_dir, file))
                if 'Sequence' in df.columns and 'Prediction' in df.columns:
                    seqs = df[df['Prediction'] == 'Non-Toxin']['Sequence'].tolist()
                    safe_from_toxin.update(seqs)
        print(f"Safe from ToxinPred (Non-Toxic): {len(safe_from_toxin)}")
    else:
        print("WARNING: toxinpred folder not found.")

    # 3. Intersection (Sequences that passed BOTH filters)
    final_safe_seqs = safe_from_allergen.intersection(safe_from_toxin)
    print(f"\nFinal safe sequences (Passed both filters): {len(final_safe_seqs)}")

    if not final_safe_seqs:
        print("ERROR: No sequences passed both filters.")
        return

    # 4. Filter original MHC-I and MHC-II data
    mhci_file = 'data/processed/final_mhci_epitopes.csv'
    mhcii_file = 'data/processed/mhcii_strong_binders.csv'
    final_dfs = []

    if os.path.exists(mhci_file):
        df_mhci = pd.read_csv(mhci_file)
        pep_col = [c for c in df_mhci.columns if 'peptide' in c.lower()][0]
        df_mhci_safe = df_mhci[df_mhci[pep_col].isin(final_safe_seqs)].copy()
        df_mhci_safe['T_Cell_Class'] = 'CD8+ (MHC-I)'
        final_dfs.append(df_mhci_safe)
        print(f"MHC-I candidates retained: {len(df_mhci_safe)}")

    if os.path.exists(mhcii_file):
        df_mhcii = pd.read_csv(mhcii_file)
        pep_col = [c for c in df_mhcii.columns if 'peptide' in c.lower()][0]
        df_mhcii_safe = df_mhcii[df_mhcii[pep_col].isin(final_safe_seqs)].copy()
        df_mhcii_safe['T_Cell_Class'] = 'CD4+ (MHC-II)'
        final_dfs.append(df_mhcii_safe)
        print(f"MHC-II candidates retained: {len(df_mhcii_safe)}")

    # 5. Combine and save
    if final_dfs:
        final_df = pd.concat(final_dfs, ignore_index=True)
        output_path = 'data/processed/final_vaccine_candidates.csv'
        final_df.to_csv(output_path, index=False, encoding='utf-8')
        print(f"\nSuccess! Final publication-ready file saved to: {output_path}")
        print("This file contains your ultimate, fully validated vaccine candidates.")
    else:
        print("ERROR: Could not find original MHC binding files to merge.")

if __name__ == "__main__":
    main()