import os
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def build_construct(mhci_peps, mhcii_peps, include_histag=False):
    parts = []
    if mhci_peps:
        parts.append(mhci_peps[0])
        for p in mhci_peps[1:]:
            parts.extend(["AAY", p])
    if mhcii_peps:
        if parts:
            parts.append("GPGPG")
        parts.append(mhcii_peps[0])
        for p in mhcii_peps[1:]:
            parts.extend(["GPGPG", p])
    
    seq = "".join(parts)
    if include_histag:
        seq += "HHHHHH"
    return seq

def get_stats(seq):
    X = ProteinAnalysis(seq)
    return {
        'II': X.instability_index(),
        'GRAVY': X.gravy(),
        'pI': X.isoelectric_point(),
        'MW': X.molecular_weight()
    }

def calculate_score(stats):
    """Multi-objective score: lower is better."""
    # Target: II < 40, GRAVY < 0, pI between 7-9
    ii_penalty = max(0, stats['II'] - 40) * 2  # Heavy penalty for II > 40
    gravy_penalty = max(0, stats['GRAVY']) * 10  # Heavy penalty for GRAVY > 0
    pi_penalty = abs(stats['pI'] - 8) * 0.5  # Prefer pI around 8
    return ii_penalty + gravy_penalty + pi_penalty

def main():
    print("="*70)
    print("IMPROVED CONSTRUCT OPTIMIZER (II, GRAVY, pI Multi-Objective)")
    print("="*70)
    
    input_csv = 'data/processed/final_vaccine_candidates.csv'
    if not os.path.exists(input_csv):
        print(f"ERROR: {input_csv} not found.")
        return
        
    df = pd.read_csv(input_csv)
    
    # Separate MHC-I and MHC-II
    df_mhci = df[df['T_Cell_Class'] == 'CD8+ (MHC-I)'].copy()
    df_mhcii = df[df['T_Cell_Class'] == 'CD4+ (MHC-II)'].copy()
    
    # Find peptide column
    pep_col = [c for c in df.columns if 'peptide' in str(c).lower()]
    pep_col = pep_col[0] if pep_col else 'Peptide'
    
    # Clean NaNs
    df_mhci = df_mhci[df_mhci[pep_col].notna()]
    df_mhcii = df_mhcii[df_mhcii[pep_col].notna()]
    
    # Sort pools
    mhci_rank_col = [c for c in df_mhci.columns if 'rank' in str(c).lower()]
    if mhci_rank_col:
        df_mhci = df_mhci.sort_values(by=mhci_rank_col[0], ascending=True)
        
    nb_col = [c for c in df_mhcii.columns if 'nb' in str(c).lower()]
    rank_col = [c for c in df_mhcii.columns if 'rank' in str(c).lower()]
    if nb_col:
        df_mhcii = df_mhcii.sort_values(by=nb_col[0], ascending=False)
    if rank_col:
        df_mhcii = df_mhcii.sort_values(by=rank_col[0], ascending=True)
        
    top_n = 8
    current_mhci = df_mhci.head(top_n)[pep_col].tolist()
    current_mhcii = df_mhcii.head(top_n)[pep_col].tolist()
    
    # Unused pools
    unused_mhci = df_mhci[~df_mhci[pep_col].isin(current_mhci)].copy()
    unused_mhcii = df_mhcii[~df_mhcii[pep_col].isin(current_mhcii)].copy()
    
    print("\nStarting multi-objective optimization (NO His-tag)...")
    max_iter = 50
    best_score = float('inf')
    best_mhci = current_mhci[:]
    best_mhcii = current_mhcii[:]
    
    for i in range(max_iter):
        seq = build_construct(current_mhci, current_mhcii, include_histag=False)
        stats = get_stats(seq)
        score = calculate_score(stats)
        
        if i % 5 == 0:
            print(f"Iteration {i+1:2d} | II: {stats['II']:6.2f} | GRAVY: {stats['GRAVY']:6.3f} | pI: {stats['pI']:5.2f} | Score: {score:.2f}")
        
        # Track best
        if score < best_score:
            best_score = score
            best_mhci = current_mhci[:]
            best_mhcii = current_mhcii[:]
        
        # Check if targets are met
        if stats['II'] < 40 and stats['GRAVY'] < 0 and 7 <= stats['pI'] <= 9:
            print("\n✅ TARGET MET: Construct is stable, hydrophilic, and neutral pI!")
            break
            
        # Find the most problematic epitope (highest individual GRAVY or most basic)
        current_mhci_scores = [(p, ProteinAnalysis(p).gravy(), ProteinAnalysis(p).isoelectric_point()) for p in current_mhci]
        current_mhcii_scores = [(p, ProteinAnalysis(p).gravy(), ProteinAnalysis(p).isoelectric_point()) for p in current_mhcii]
        
        # Prioritize swapping epitopes with high pI (basic) or high GRAVY (hydrophobic)
        max_mhci = max(current_mhci_scores, key=lambda x: x[2] + x[1]*5) if current_mhci_scores else ("", 0, 0)
        max_mhcii = max(current_mhcii_scores, key=lambda x: x[2] + x[1]*5) if current_mhcii_scores else ("", 0, 0)
        
        swapped = False
        if max_mhci[2] > max_mhcii[2] and not unused_mhci.empty:
            culprit = max_mhci[0]
            # Find replacement with lowest pI and negative GRAVY
            replacement_row = unused_mhci.iloc[0]
            replacement = replacement_row[pep_col]
            current_mhci.remove(culprit)
            current_mhci.append(replacement)
            unused_mhci = unused_mhci.iloc[1:]
            print(f"  -> Swapped MHC-I: {culprit} (pI {max_mhci[2]:.2f}) with {replacement}")
            swapped = True
        elif not unused_mhcii.empty:
            culprit = max_mhcii[0]
            replacement_row = unused_mhcii.iloc[0]
            replacement = replacement_row[pep_col]
            current_mhcii.remove(culprit)
            current_mhcii.append(replacement)
            unused_mhcii = unused_mhcii.iloc[1:]
            print(f"  -> Swapped MHC-II: {culprit} (pI {max_mhcii[2]:.2f}) with {replacement}")
            swapped = True
            
        if not swapped:
            print("No more replacements available.")
            break
    
    # Use best found
    current_mhci = best_mhci
    current_mhcii = best_mhcii
    
    # Final save (NO His-tag)
    final_seq = build_construct(current_mhci, current_mhcii, include_histag=False)
    final_stats = get_stats(final_seq)
    
    out_seq = 'data/processed/mev_construct_optimized.txt'
    with open(out_seq, 'w') as f:
        f.write(final_seq)
        
    out_csv = 'data/processed/top_selected_epitopes_optimized.csv'
    final_df = pd.DataFrame({
        'Peptide': current_mhci + current_mhcii,
        'Class': ['MHC-I']*len(current_mhci) + ['MHC-II']*len(current_mhcii)
    })
    final_df.to_csv(out_csv, index=False)
    
    print("\n" + "="*70)
    print("OPTIMIZATION COMPLETE")
    print("="*70)
    print(f"Final Construct Length: {len(final_seq)} aa (NO His-tag)")
    print(f"Final Instability Index: {final_stats['II']:.2f} (Target < 40)")
    print(f"Final GRAVY: {final_stats['GRAVY']:.3f} (Target < 0)")
    print(f"Final pI: {final_stats['pI']:.2f} (Target 7-9)")
    print(f"\nOptimized sequence saved to: {out_seq}")
    print(f"Optimized epitope list saved to: {out_csv}")

if __name__ == "__main__":
    main()