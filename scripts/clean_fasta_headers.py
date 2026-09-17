import os
from Bio import SeqIO

input_dir = "data/processed/batches"

for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        filepath = os.path.join(input_dir, filename)
        records = []
        for rec in SeqIO.parse(filepath, "fasta"):
            # Space र special character हटाएर underscore राख्ने
            clean_id = rec.id.replace(" ", "_").replace("|", "_").replace(":", "_")
            rec.id = clean_id
            rec.description = "" # Description हटाउँदा parser crash हुँदैन
            records.append(rec)
        
        with open(filepath, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")
        print(f"Cleaned: {filename}")

print("All FASTA files are now IEDB MHC-II safe.")