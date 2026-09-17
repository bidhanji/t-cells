import os

def create_chimerax_script():
    complex_pdb = "data/processed/docked_complex.pdb"
    
    if not os.path.exists(complex_pdb):
        print("ERROR: docked_complex.pdb not found in data/processed/")
        print("Please make sure the file exists before running this script.")
        return
    
    abs_pdb = os.path.abspath(complex_pdb)
    abs_pdb = abs_pdb.replace("\\", "/")
    
    output_dir = os.path.abspath("data/processed").replace("\\", "/")
    
    cxc_script = f"""# ChimeraX Visualization Script for MEV-TLR4 Docking
# Generated automatically

# Open the docked complex
open "{abs_pdb}"

# Wait for model to load
wait 2

# Color receptor (Chain A) blue and show as surface
color /A dodgerblue
surface /A
transparency /A 30

# Color ligand (Chain B) red and show as sticks
color /B red
style /B stick

# Set background to white for publication
set bgColor white

# Improve lighting
lighting soft

# Focus on the complex
view

# Save high-resolution image
save "{output_dir}/Figure3_MEV_TLR4_Docking.png" width 3000 height 3000 supersample 3

# Save a second angle (rotated 90 degrees)
turn y 90
save "{output_dir}/Figure3_MEV_TLR4_Docking_Side.png" width 3000 height 3000 supersample 3

# Save a third angle (top view)
turn y 90
turn x 90
save "{output_dir}/Figure3_MEV_TLR4_Docking_Top.png" width 3000 height 3000 supersample 3

# Print confirmation
echo "All images saved successfully to {output_dir}/"
"""
    
    cxc_path = "data/processed/visualize_docking.cxc"
    with open(cxc_path, "w", encoding="utf-8") as f:
        f.write(cxc_script)
    
    print("=" * 70)
    print("CHIMERAX SCRIPT CREATED SUCCESSFULLY")
    print("=" * 70)
    print(f"Script saved to: {os.path.abspath(cxc_path)}")
    print()
    print("NOW DO THIS:")
    print("1. Open ChimeraX")
    print("2. Click File -> Open")
    print(f"3. Navigate to: {os.path.abspath('data/processed')}")
    print("4. Select: visualize_docking.cxc")
    print("5. Click Open")
    print()
    print("ChimeraX will automatically:")
    print("  - Load the docked complex")
    print("  - Color receptor blue (surface)")
    print("  - Color ligand red (sticks)")
    print("  - Set white background")
    print("  - Save 3 high-resolution images (3000x3000)")
    print()
    print("Images will be saved as:")
    print(f"  - {output_dir}/Figure3_MEV_TLR4_Docking.png")
    print(f"  - {output_dir}/Figure3_MEV_TLR4_Docking_Side.png")
    print(f"  - {output_dir}/Figure3_MEV_TLR4_Docking_Top.png")
    print("=" * 70)

if __name__ == "__main__":
    create_chimerax_script()