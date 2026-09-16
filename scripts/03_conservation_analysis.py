import os
import csv
from Bio import SeqIO
from Bio import pairwise2
import warnings

# Suppress Biopython deprecation warnings for pairwise2
warnings.filterwarnings("ignore", category=DeprecationWarning)

def find_orthologs(query_seq, target_dict, min_identity=90.0, min_coverage=80.0):
    best_hit = None
    best_identity = 0.0
    
    for target_id, target_seq in target_dict.items():
        # globalxx calculates alignment based purely on matches, ignoring gap penalties
        alignments = pairwise2.align.globalxx(query_seq, target_seq)
        if not alignments:
            continue
            
        best_aln = alignments[0]
        matches = best_aln.score
        align_len = len(best_aln.seqA)
        
        identity = (matches / align_len) * 100 if align_len > 0 else 0
        query_gaps = best_aln.seqA.count('-')
        coverage = ((align_len - query_gaps) / len(query_seq)) * 100
        
        if identity >= min_identity and coverage >= min_coverage:
            if identity > best_identity:
                best_identity = identity
                best_hit = {
                    "target_id": target_id,
                    "identity": round(identity, 2),
                    "coverage": round(coverage, 2)
                }
                
    return best_hit

def main():
    input_dir = "data/processed"
    output_file = os.path.join(input_dir, "conserved_pan_capriproteins.csv")
    
    print("Loading protein sequences...")
    lsdv_proteins = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(input_dir, "LSDV_proteins.fasta"), "fasta")}
    sppv_proteins = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(input_dir, "SPPV_proteins.fasta"), "fasta")}
    gtpv_proteins = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(input_dir, "GTPV_proteins.fasta"), "fasta")}
    
    print(f"Loaded: LSDV ({len(lsdv_proteins)}), SPPV ({len(sppv_proteins)}), GTPV ({len(gtpv_proteins)})")
    
    conserved_results = []
    total = len(lsdv_proteins)
    
    print("Starting pairwise conservation analysis (this may take 1-2 minutes)...")
    for i, (lsdv_id, lsdv_seq) in enumerate(lsdv_proteins.items()):
        if (i + 1) % 20 == 0:
            print(f"  -> Processed {i + 1}/{total} LSDV proteins...")
            
        sppv_hit = find_orthologs(lsdv_seq, sppv_proteins)
        if not sppv_hit:
            continue
            
        gtpv_hit = find_orthologs(lsdv_seq, gtpv_proteins)
        if not gtpv_hit:
            continue
            
        conserved_results.append({
            "LSDV_Ortholog": lsdv_id,
            "SPPV_Ortholog": sppv_hit["target_id"],
            "SPPV_Identity_%": sppv_hit["identity"],
            "GTPV_Ortholog": gtpv_hit["target_id"],
            "GTPV_Identity_%": gtpv_hit["identity"]
        })
        
    print(f"\nAnalysis complete. Found {len(conserved_results)} pan-capripoxvirus conserved proteins.")
    
    if conserved_results:
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=conserved_results[0].keys())
            writer.writeheader()
            writer.writerows(conserved_results)
        print(f"Results saved to {output_file}")
    else:
        print("No conserved proteins found. Check your thresholds.")

if __name__ == "__main__":
    main()