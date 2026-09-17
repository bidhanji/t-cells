import os
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def filter_promiscuous_epitopes(mhci_file, mhcii_file, output_fasta, min_alleles=2):
    print(f"Filtering for peptides binding to >= {min_alleles} alleles...")
    
    # Dictionary to store peptide -> set of alleles
    peptide_alleles = {}
    
    # 1. Process MHC-I
    if os.path.exists(mhci_file):
        df_mhci = pd.read_csv(mhci_file)
        pep_col = [c for c in df_mhci.columns if 'peptide' in c.lower()][0]
        # MHC-I file मा 'allele' वा 'Allele' column हुन्छ
        allele_col = [c for c in df_mhci.columns if 'allele' in c.lower()]
        if allele_col:
            allele_col = allele_col[0]
            for _, row in df_mhci.iterrows():
                pep = str(row[pep_col]).strip()
                alle = str(row[allele_col]).strip()
                if pep and alle and pep != 'nan':
                    peptide_alleles.setdefault(pep, set()).add(alle)
        print(f"Processed MHC-I: {len(peptide_alleles)} unique peptides so far.")

    # 2. Process MHC-II
    if os.path.exists(mhcii_file):
        df_mhcii = pd.read_csv(mhcii_file)
        pep_col = [c for c in df_mhcii.columns if 'peptide' in c.lower()][0]
        # MHC-II file मा 'Allele' column हुन्छ
        allele_col = [c for c in df_mhcii.columns if 'allele' in c.lower()]
        if allele_col:
            allele_col = allele_col[0]
            for _, row in df_mhcii.iterrows():
                pep = str(row[pep_col]).strip()
                alle = str(row[allele_col]).strip()
                if pep and alle and pep != 'nan':
                    peptide_alleles.setdefault(pep, set()).add(alle)
        print(f"Processed MHC-II. Total unique peptides: {len(peptide_alleles)}")

    # 3. Filter by promiscuity
    promiscuous_peptides = {pep: alles for pep, alles in peptide_alleles.items() if len(alles) >= min_alleles}
    
    if not promiscuous_peptides:
        print(f"ERROR: No peptides found binding to >= {min_alleles} alleles.")
        return

    print(f"\nFound {len(promiscuous_peptides)} promiscuous epitopes (binding to >= {min_alleles} alleles).")
    
    # 4. Write to FASTA
    with open(output_fasta, 'w', encoding='utf-8') as f:
        for i, (pep, alles) in enumerate(sorted(promiscuous_peptides.items()), 1):
            f.write(f">Promiscuous_{i}_Alleles:{len(alles)}\n")
            f.write(f"{pep}\n")
            
    print(f"Saved to: {output_fasta}")
    print("NEXT ACTION: We will now automate AllerTOP and ToxinPred for this reduced list.")

if __name__ == "__main__":
    mhci_path = "data/processed/final_mhci_epitopes.csv"
    mhcii_path = "data/processed/mhcii_strong_binders.csv"
    output_path = "data/processed/promiscuous_epitopes.fasta"
    
    # कम्तीमा २ वटा allele मा bind हुने epitopes मात्र छुट्याउने
    filter_promiscuous_epitopes(mhci_path, mhcii_path, output_path, min_alleles=2)