import os
import platform
from pathlib import Path

def get_true_desktop_path():
    home = Path.home()
    
    # Check for Windows OneDrive redirection, a common trap
    if platform.system() == "Windows":
        onedrive_desktop = home / "OneDrive" / "Desktop"
        local_desktop = home / "Desktop"
        
        if onedrive_desktop.exists():
            print("Warning: OneDrive Desktop detected. Routing to local Desktop to prevent sync locks.")
            # Force local desktop to avoid cloud sync interference during heavy I/O operations
            return local_desktop if local_desktop.exists() else onedrive_desktop
            
    # Standard fallback for macOS and Linux
    return home / "Desktop"

def build_research_environment():
    desktop_path = get_true_desktop_path()
    project_root = desktop_path / "Capripox_TCell_Epitope_Discovery"
    
    # Define the exact hierarchy required for reproducible computational biology
    directories = [
        project_root / "config",
        project_root / "data" / "raw",          # Unmodified GenBank downloads go here
        project_root / "data" / "processed",    # Cleaned FASTA and CSV outputs go here
        project_root / "scripts",               # Your Python automation code
        project_root / "tables",                # Final manuscript tables
        project_root / "figures",               # High-resolution plots (PDF/SVG)
        project_root / "paper"                  # Manuscript drafts and references
    ]
    
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
        print(f"Created: {directory}")

    # Create a strict .gitignore to prevent uploading massive genomic data to GitHub
    gitignore_content = """# Ignore all raw genomic data
data/raw/*
!data/raw/.gitkeep

# Ignore processed large files
*.gb
*.fasta
*.csv

# Python artifacts
__pycache__/
*.py[cod]
*.so
.env
venv/
.DS_Store
"""
    with open(project_root / ".gitignore", "w", encoding="utf-8") as f:
        f.write(gitignore_content)

    # Create a centralized configuration file for your pipeline thresholds
    config_content = """# Centralized pipeline parameters
# Modify these values here, not in your Python scripts
accessions:
  lsdv_ref: "MN072619.1"
  sppv_ref: "AY077833"
  gtpv_ref: "AY077835"

thresholds:
  conservation_identity: 90.0
  mhc1_percent_rank: 0.5
  mhc2_percent_rank: 2.0
  proteasome_cleavage: 0.50
  tap_transport: 0.05
  vaxijen_score: 0.40

target_alleles:
  bola_class1:
    - "BoLA-1:02301"
    - "BoLA-2:01201"
    - "BoLA-A11"
    - "BoLA-T2c"
  bola_class2:
    - "BoLA-DRB3*00201"
    - "BoLA-DRB3*01101"
"""
    with open(project_root / "config" / "config.yaml", "w", encoding="utf-8") as f:
        f.write(config_content)

    # Create placeholder files to enforce disciplined documentation
    readme_content = f"""# Capripox T-Cell Epitope Discovery Pipeline
Location: {project_root}

## Execution Order
1. Run `scripts/01_fetch_genomes.py` to populate `data/raw/`.
2. Run `scripts/02_filter_and_align.py` to generate `data/processed/` files.
3. Update `tables/` with final candidate epitopes before drafting the `paper/`.
"""
    with open(project_root / "README.md", "w", encoding="utf-8") as f:
        f.write(readme_content)

    # Create empty marker files so Git tracks the directories
    (project_root / "data" / "raw" / ".gitkeep").touch()
    (project_root / "data" / "processed" / ".gitkeep").touch()
    (project_root / "figures" / ".gitkeep").touch()

    print(f"\nSuccess. Environment built at: {project_root}")
    print("नेपालको इन्टरनेट र पावर कटको अवस्थालाई ध्यानमा राख्दै, सधैं data/raw फोल्डरको बैकअप राख्नुहोस्।") 
    # [Translation: Keeping the reality of Nepal's internet and power cuts in mind, always keep a backup of the data/raw folder.]

if __name__ == "__main__":
    build_research_environment()