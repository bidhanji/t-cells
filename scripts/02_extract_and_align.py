import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_proteins_from_genbank(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # STRICT MAPPING: Only exact matches allowed
    target_accessions = {
        "MN072619.1": "LSDV",
        "NC_004002.1": "SPPV",
        "NC_004003.1": "GTPV"
    }
    
    for filename in os.listdir(input_dir):
        if not filename.endswith(".gb"):
            continue
            
        # CRITICAL FIX: Robust filename parsing
        # Filename is like "NC_004002_1.gb". We need "NC_004002.1"
        accession_base = filename.replace(".gb", "")
        
        # Find the LAST underscore (which separates the version number) and replace only that with a dot
        if "_" in accession_base:
            last_underscore_idx = accession_base.rfind("_")
            accession_base = accession_base[:last_underscore_idx] + "." + accession_base[last_underscore_idx+1:]
            
        if accession_base not in target_accessions:
            print(f"  -> SKIPPED: {filename} (Parsed as {accession_base}, not in target list)")
            continue
            
        matched_virus = target_accessions[accession_base]
        output_fasta = os.path.join(output_dir, f"{matched_virus}_proteins.fasta")
        print(f"Processing {filename} as {matched_virus}...")
        
        protein_records = []
        with open(os.path.join(input_dir, filename), "r", encoding="utf-8") as gb_handle:
            for record in SeqIO.parse(gb_handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        prot_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
                        locus_tag = feature.qualifiers.get("locus_tag", ["no_tag"])[0]
                        gene_name = feature.qualifiers.get("gene", ["no_gene"])[0]
                        
                        seq_id = f"{matched_virus}|{locus_tag}|{gene_name}|{prot_id}"
                        seq_desc = feature.qualifiers.get("product", ["hypothetical protein"])[0]
                        
                        if "translation" in feature.qualifiers:
                            prot_seq = feature.qualifiers["translation"][0]
                        else:
                            try:
                                nuc_seq = feature.location.extract(record.seq)
                                prot_seq = str(nuc_seq.translate())
                                if prot_seq.endswith("*"):
                                    prot_seq = prot_seq[:-1]
                            except Exception as e:
                                continue
                        
                        if prot_seq.count("*") > 1 or len(prot_seq) < 50:
                            continue
                            
                        protein_records.append(SeqRecord(Seq(prot_seq), id=seq_id, description=seq_desc))
            
        if protein_records:
            with open(output_fasta, "w", encoding="utf-8") as out_handle:
                SeqIO.write(protein_records, out_handle, "fasta")
            print(f"  -> Extracted {len(protein_records)} proteins to {matched_virus}_proteins.fasta")
        else:
            print(f"  -> WARNING: No valid CDS found in {filename}.")

if __name__ == "__main__":
    extract_proteins_from_genbank("data/raw", "data/processed")
    print("ORF extraction complete.")