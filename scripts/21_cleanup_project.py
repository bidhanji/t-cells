import os
import shutil
from pathlib import Path

def organize_and_cleanup():
    project_root = Path(".")
    processed_dir = project_root / "data" / "processed"
    
    # 1. Create organized folder for safety results
    safety_results_dir = processed_dir / "safety_results"
    safety_results_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("ORGANIZING SAFETY RESULTS & CLEANING REDUNDANT FILES")
    print("="*70)
    
    # 2. Move safety result folders to organized location
    folders_to_move = {
        'allertop': 'MHC-I AllerCatPro results',
        'toxinpred': 'MHC-I ToxinPred results',
        'toxinpredmhcii': 'MHC-II ToxinPred results'
    }
    
    print("\n📁 Organizing safety screening results...")
    for folder_name, description in folders_to_move.items():
        src = processed_dir / folder_name
        dst = safety_results_dir / folder_name
        
        if src.exists() and src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.move(str(src), str(dst))
            print(f"  ✅ Moved {folder_name} -> safety_results/{folder_name}")
        else:
            print(f"  ⚠️  {folder_name} not found (already moved or doesn't exist)")
    
    # 3. Delete intermediate batch folders
    batch_folders_to_delete = [
        'ultimate_web_safe_batches',
        'mhcii_safety_batches',
        'mhcii_allercatpro_batches',
        'mhcii_toxinpred_batches',
        'safety_batches',
        'safety_batches_final',
        'web_ready_batches',
        'mhcii_micro_batches'
    ]
    
    print("\n🗑️  Deleting intermediate batch folders...")
    for folder_name in batch_folders_to_delete:
        folder_path = processed_dir / folder_name
        if folder_path.exists() and folder_path.is_dir():
            try:
                shutil.rmtree(folder_path)
                print(f"  ✅ Deleted folder: {folder_name}")
            except PermissionError:
                print(f"  ⚠️  Could not delete {folder_name} (Close any open files in it)")
    
    # 4. Delete intermediate files
    files_to_delete = [
        'unique_strong_binders.fasta',
        'promiscuous_epitopes.fasta',
        'mhcii_unique_peptides.fasta',
        'mhcii_unique_peptides_clean.fasta',
        'mev_construct_sequence.txt',
        'mev_construct_optimized.txt',
        'top_selected_epitopes.csv',
        'top_selected_epitopes_optimized.csv',
        'mhcii_local_safety_results.csv',
        'mhcii_safe_by_motif.txt'
    ]
    
    print("\n🗑️  Deleting intermediate files...")
    for file_name in files_to_delete:
        file_path = processed_dir / file_name
        if file_path.exists() and file_path.is_file():
            try:
                file_path.unlink()
                print(f"  ✅ Deleted file: {file_name}")
            except PermissionError:
                print(f"  ⚠️  SKIPPED: {file_name} (CLOSE THIS FILE IN EXCEL/NOTEPAD FIRST!)")
    
    # 5. Summary
    print("\n" + "="*70)
    print("CLEANUP COMPLETE")
    print("="*70)
    print("✅ SAFETY RESULTS ORGANIZED IN: data/processed/safety_results/")
    print("✅ REDUNDANT BATCHES REMOVED")
    print("\n⚠️  If any file says 'SKIPPED', close it in Excel/Notepad and run this script again.")
    print("="*70)

if __name__ == "__main__":
    organize_and_cleanup()