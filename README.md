# Capripox T-Cell Epitope Discovery Pipeline
Location: C:\Users\ACER\Desktop\Capripox_TCell_Epitope_Discovery

## Execution Order
1. Run `scripts/01_fetch_genomes.py` to populate `data/raw/`.
2. Run `scripts/02_filter_and_align.py` to generate `data/processed/` files.
3. Update `tables/` with final candidate epitopes before drafting the `paper/`.
