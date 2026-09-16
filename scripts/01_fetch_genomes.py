import os
import time
from Bio import Entrez
from Bio import SeqIO

# Register your identity with NCBI. 
# Replace with your actual NCBI API key if you have one to increase rate limits.
Entrez.email = "bidhanji@gmail.com"
# Entrez.api_key = "YOUR_API_KEY_HERE"

def fetch_and_save_genomes(search_query, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Executing search: {search_query}")
    search_handle = Entrez.esearch(db="nucleotide", term=search_query, usehistory="y", retmax=10000)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    print(f"Found {count} matching records. Starting individual file downloads.")
    
    if count == 0:
        print("No records found. Check your search query.")
        return

    batch_size = 50
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
                    # Sanitize filename to prevent OS errors
                    accession = record.id.replace(".", "_").replace("/", "_")
                    filename = os.path.join(output_dir, f"{accession}.gb")
                    with open(filename, "w") as f:
                        SeqIO.write(record, f, "genbank")
                        
                fetch_handle.close()
                print(f"Successfully downloaded batch {start + 1} to {end}.")
                break
                
            except Exception as e:
                attempt += 1
                print(f"Network error at batch {start + 1} to {end}. Retrying attempt {attempt}. Error: {e}")
                time.sleep(5 * attempt)
                
        if attempt == max_attempts:
            print(f"Failed to download batch {start + 1} to {end} after {max_attempts} attempts. Skipping.")
            
        # Respect NCBI rate limits (3 requests per second without API key)
        time.sleep(0.4)

if __name__ == "__main__":
    # This query specifically targets complete genomes, avoiding partial fragments
    query = '"Lumpy skin disease virus"[Organism] AND "complete genome"[Filter]'
    target_directory = "data/raw"
    
    fetch_and_save_genomes(query, target_directory)
    print("Download process completed.")