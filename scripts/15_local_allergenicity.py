"""
Local allergenicity screening - FAO/WHO 6-mer sliding window method.

PIPELINE:
  Step 1: Extract allergen accession numbers from AllergenOnline PDF
  Step 2: Fetch actual protein sequences from NCBI (requires internet)
  Step 3: Build 6-mer lookup table from all allergen sequences
  Step 4: Screen MHC-II candidate peptides against the table
  Step 5: Output safe and flagged CSVs

REQUIREMENTS:
  pip install biopython pdfminer.six

INPUT FILES NEEDED:
  data/external/AllergenOnlineV24.pdf   -- the PDF you already downloaded
  data/processed/mhcii_promiscuous_prefilt.csv  -- from 15_prepare_mhcii_safety.py

CITATION FOR METHODS SECTION:
  "Allergenicity of candidate MHC-II epitopes was assessed using the FAO/WHO
  (2001) sliding window criterion. All 6-mer subsequences [contiguous windows
  of 6 amino acids] of each 15-mer candidate were compared against the full
  protein sequences of 2,373 allergens catalogued in the AllergenOnline
  database (v24, January 2026). Peptides containing any exact 6-mer match
  to a known allergen sequence were excluded from further analysis."
"""

import os
import re
import sys
import time
import subprocess

import pandas as pd

# ---- config ------------------------------------------------------------
ALLERGEN_PDF       = "data/external/AllergenOnlineV24.pdf"
ALLERGEN_FASTA_OUT = "data/external/allergen_sequences.fasta"   # cached after first run
INPUT_CSV          = "data/processed/mhcii_promiscuous_prefilt.csv"
OUTPUT_SAFE_CSV    = "data/processed/mhcii_allergen_safe.csv"
OUTPUT_FLAGGED_CSV = "data/processed/mhcii_allergen_flagged.csv"
PEPTIDE_COL        = "Peptide"
KMER_SIZE          = 6
NCBI_BATCH_SIZE    = 200   # max accessions per NCBI fetch call
NCBI_EMAIL         = "your_email@example.com"   # CHANGE THIS -- NCBI requires it
NCBI_SLEEP         = 0.4   # seconds between NCBI calls -- stay under their rate limit
# ------------------------------------------------------------------------


def extract_accessions_from_pdf(pdf_path):
    """Extract NCBI/UniProt accession numbers from AllergenOnline PDF."""
    try:
        result = subprocess.run(
            ["pdftotext", pdf_path, "-"],
            capture_output=True, text=True, check=True
        )
        text = result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        # fallback: try pdfminer
        try:
            from pdfminer.high_level import extract_text
            text = extract_text(pdf_path)
        except ImportError:
            print("ERROR: pdftotext not found and pdfminer.six not installed.")
            print("Install pdfminer: pip install pdfminer.six")
            sys.exit(1)

    # Match NCBI GenPept (ABL09307.1) and UniProt (P00785.4, A5HII1.1) accessions
    pattern = re.compile(r'\b([A-Z]{1,3}[0-9][A-Z0-9]{4,}\.[0-9]+)\b')
    accessions = list(dict.fromkeys(pattern.findall(text)))  # deduplicate, preserve order
    return accessions


def fetch_sequences_from_ncbi(accessions, email, batch_size, sleep_sec, out_fasta):
    """Fetch protein FASTA sequences from NCBI Entrez for all accessions."""
    from Bio import Entrez, SeqIO
    Entrez.email = email

    all_records = []
    total = len(accessions)
    print(f"  Fetching {total} sequences from NCBI in batches of {batch_size}...")
    print("  This takes 3-8 minutes depending on your connection. Do not interrupt.")

    for i in range(0, total, batch_size):
        batch = accessions[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size
        print(f"  Batch {batch_num}/{total_batches} ({len(batch)} accessions)...", end=" ")

        try:
            handle = Entrez.efetch(
                db="protein",
                id=",".join(batch),
                rettype="fasta",
                retmode="text"
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            all_records.extend(records)
            print(f"got {len(records)} sequences")
        except Exception as e:
            print(f"WARNING: batch {batch_num} failed ({e}) -- skipping")

        time.sleep(sleep_sec)

    # Cache to disk so re-runs don't refetch
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    with open(out_fasta, "w") as f:
        SeqIO.write(all_records, f, "fasta")
    print(f"  Sequences cached to {out_fasta} ({len(all_records)} total)")
    return all_records


def load_fasta_sequences(fasta_path):
    """Load sequences from a FASTA file (used for cached re-runs)."""
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_path, "fasta"))
    return records


def build_kmer_set(records, k):
    """Build set of all k-mers from all allergen sequences."""
    kmer_set = set()
    for rec in records:
        seq = str(rec.seq).upper().replace("-", "").replace("*", "")
        for i in range(len(seq) - k + 1):
            kmer_set.add(seq[i:i + k])
    return kmer_set


def screen_peptides(df, peptide_col, kmer_set, k):
    flags = []
    matched_kmers = []
    for pep in df[peptide_col]:
        pep = str(pep).upper().strip()
        hits = []
        for i in range(len(pep) - k + 1):
            window = pep[i:i + k]
            if window in kmer_set:
                hits.append(window)
        flags.append(len(hits) > 0)
        matched_kmers.append(";".join(hits) if hits else "")
    return flags, matched_kmers


def main():
    script_dir   = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, ".."))

    pdf_path      = os.path.join(project_root, ALLERGEN_PDF)
    fasta_cache   = os.path.join(project_root, ALLERGEN_FASTA_OUT)
    input_path    = os.path.join(project_root, INPUT_CSV)
    safe_path     = os.path.join(project_root, OUTPUT_SAFE_CSV)
    flagged_path  = os.path.join(project_root, OUTPUT_FLAGGED_CSV)

    # Validate inputs
    if not os.path.isfile(pdf_path):
        print(f"ERROR: AllergenOnline PDF not found: {pdf_path}")
        print("Place AllergenOnlineV24.pdf in data/external/")
        sys.exit(1)

    if not os.path.isfile(input_path):
        print(f"ERROR: Peptide input not found: {input_path}")
        print("Run 15_prepare_mhcii_safety.py first.")
        sys.exit(1)

    if NCBI_EMAIL == "bidhanji@gmail.com":
        print("ERROR: Set NCBI_EMAIL in the config section to your real email.")
        print("NCBI requires this for Entrez API access.")
        sys.exit(1)

    # Step 1: Get allergen sequences (use cache if available)
    if os.path.isfile(fasta_cache):
        print(f"Cached allergen sequences found: {fasta_cache}")
        print("Loading from cache (delete this file to force re-fetch from NCBI)...")
        records = load_fasta_sequences(fasta_cache)
        print(f"  {len(records)} sequences loaded")
    else:
        print("Step 1: Extracting accession numbers from AllergenOnline PDF...")
        accessions = extract_accessions_from_pdf(pdf_path)
        print(f"  {len(accessions)} unique accessions found")

        print("Step 2: Fetching sequences from NCBI Entrez...")
        records = fetch_sequences_from_ncbi(
            accessions, NCBI_EMAIL, NCBI_BATCH_SIZE, NCBI_SLEEP, fasta_cache
        )

    if not records:
        print("ERROR: No sequences loaded. Check internet connection and NCBI_EMAIL.")
        sys.exit(1)

    # Step 3: Build 6-mer table
    print(f"\nBuilding {KMER_SIZE}-mer lookup table...")
    kmer_set = build_kmer_set(records, KMER_SIZE)
    print(f"  {len(kmer_set):,} unique {KMER_SIZE}-mers indexed from {len(records)} allergen sequences")

    # Step 4: Load and screen peptides
    df = pd.read_csv(input_path)
    print(f"\nScreening {len(df)} MHC-II candidate peptides...")

    if PEPTIDE_COL not in df.columns:
        print(f"ERROR: Column '{PEPTIDE_COL}' not found. Columns: {list(df.columns)}")
        sys.exit(1)

    flags, matched = screen_peptides(df, PEPTIDE_COL, kmer_set, KMER_SIZE)
    df["_allergenic_flag"]      = flags
    df["_matched_allergen_kmers"] = matched

    df_safe    = df[~df["_allergenic_flag"]].copy()
    df_flagged = df[ df["_allergenic_flag"]].copy()

    os.makedirs(os.path.dirname(safe_path), exist_ok=True)
    df_safe.to_csv(safe_path,    index=False)
    df_flagged.to_csv(flagged_path, index=False)

    total  = len(df)
    n_safe = len(df_safe)
    n_flag = len(df_flagged)

    print(f"""
Allergenicity Screening Complete (FAO/WHO 6-mer method)
-------------------------------------------------------
Allergen sequences used : {len(records)}
6-mers in lookup table  : {len(kmer_set):,}
Peptides screened       : {total}
Non-allergenic (safe)   : {n_safe}  ({100*n_safe/total:.1f}%)
Flagged (allergenic)    : {n_flag}  ({100*n_flag/total:.1f}%)

Outputs:
  Safe    --> {OUTPUT_SAFE_CSV}
  Flagged --> {OUTPUT_FLAGGED_CSV}
  (Flagged file shows exact matching 6-mer windows in _matched_allergen_kmers)

Next step: toxicity screening on {OUTPUT_SAFE_CSV}
""")


if __name__ == "__main__":
    main()