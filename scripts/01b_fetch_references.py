import os
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "bidhanji@gmail.com"

def fetch_specific_references(accessions, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    for acc in accessions:
        # Force clean filename format
        filename = os.path.join(output_dir, f"{acc.replace('.', '_')}.gb")
        
        print(f"Downloading verified RefSeq: {acc}...")
        try:
            fetch_handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(fetch_handle, "genbank")
            fetch_handle.close()
            
            # Strict length check for Capripoxvirus (~150kb)
            if len(record.seq) < 140000:
                print(f"  -> REJECTED: {acc} is only {len(record.seq)} bp. Not a complete genome.")
                if os.path.exists(filename):
                    os.remove(filename)
                continue
                
            with open(filename, "w", encoding="utf-8") as f:
                SeqIO.write(record, f, "genbank")
            print(f"  -> SUCCESS: Saved {acc} ({len(record.seq)} bp)")
            
        except Exception as e:
            print(f"  -> FAILED: {acc}. Error: {e}")

if __name__ == "__main__":
    # EXACT, VERIFIED REFSEQ ACCESSIONS
    targets = [
        "MN072619.1",  # LSDV
        "NC_004002.1", # SPPV (Sheep Pox)
        "NC_004003.1"  # GTPV (Goat Pox)
    ]
    fetch_specific_references(targets, "data/raw")
    print("RefSeq download complete.")