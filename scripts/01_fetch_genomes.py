import os
import time
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "bidhanji@gmail.com"

def fetch_and_save_genomes(search_query, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Executing search: {search_query}")
    search_handle = Entrez.esearch(db="nucleotide", term=search_query, usehistory="y", retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    print(f"Found {count} total records. Filtering for complete genomes (>= 140,000 bp) locally.")
    
    if count == 0:
        print("No records found.")
        return

    batch_size = 50
    saved_count = 0
    
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        
        attempt = 0
        max_attempts = 5
        
        while attempt < max_attempts:
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key
                )
                
                records = SeqIO.parse(fetch_handle, "genbank")
                for record in records:
                    # CRITICAL FIX: Filter by actual sequence length, not lazy metadata
                    if len(record.seq) >= 140000:
                        accession = record.id.replace(".", "_").replace("/", "_")
                        filename = os.path.join(output_dir, f"{accession}.gb")
                        with open(filename, "w") as f:
                            SeqIO.write(record, f, "genbank")
                        saved_count += 1
                        
                fetch_handle.close()
                print(f"Processed batch {start + 1} to {end}. Saved {saved_count} complete genomes so far.")
                break
                
            except Exception as e:
                attempt += 1
                print(f"Network error. Retrying attempt {attempt}. Error: {e}")
                time.sleep(5 * attempt)
                
        if attempt == max_attempts:
            print(f"Failed batch {start + 1} to {end}.")
            
        time.sleep(0.4)

    print(f"Download complete. Total complete genomes saved: {saved_count}")

if __name__ == "__main__":
    # Removed the invalid [Filter] tag. We filter by length locally.
    query = '"Lumpy skin disease virus"[Organism]'
    target_directory = "data/raw"
    
    fetch_and_save_genomes(query, target_directory)