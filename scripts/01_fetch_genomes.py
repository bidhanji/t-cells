import os
from Bio import Entrez

# Register your identity with NCBI Entrez API
Entrez.email = "bidhanji@gmail.com"
# Optional: Set Entrez.api_key = "YOUR_API_KEY" if you have an NCBI account key

def download_ncbi_records(
    search_query: str, 
    output_file: str = "sequences.gb", 
    ret_type: str = "gb", 
    batch_size: str = 100
):
    """
    Executes a search query against NCBI Nucleotide database and downloads records in batches.
    
    Parameters:
    search_query : Term to search (e.g., 'Lumpy skin disease virus[Organism] AND complete genome[Title]')
    output_file  : Local file path to save concatenated records
    ret_type     : 'gb' for full GenBank format, 'fasta' for FASTA format
    batch_size   : Number of records to pull per HTTP request
    """
    print(f"Executing search: {search_query}")
    
    # Search database and store results on NCBI server environment
    search_handle = Entrez.esearch(
        db="nucleotide", 
        term=search_query, 
        usehistory="y"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    print(f"Found {count} matching records.")
    if count == 0:
        return

    # Download in chunks to avoid HTTP 500 timeouts
    with open(output_file, "w") as out_stream:
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            print(f"Downloading records {start + 1} to {end} of {count}...")
            
            fetch_handle = Entrez.efetch(
                db="nucleotide",
                rettype=ret_type,
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key
            )
            
            data = fetch_handle.read()
            fetch_handle.close()
            out_stream.write(data)

    print(f"Finished. File saved as {output_file}")

if __name__ == "__main__":
    # Specify your target pathogen, locus, or host query string
    QUERY = '"Lumpy skin disease virus"[Organism] AND complete genome[Title]'
    
    # Save as full GenBank flat file format
    download_ncbi_records(
        search_query=QUERY, 
        output_file="lsdv_complete_genomes.gb", 
        ret_type="gb", 
        batch_size=50
    )