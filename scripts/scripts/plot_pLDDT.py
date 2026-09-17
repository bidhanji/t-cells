import numpy as np
import matplotlib.pyplot as plt

pdb = r"C:\Users\ACER\Desktop\Capripox_TCell_Epitope_Discovery\data\processed\mev_construct_3D.pdb"

resids, bvals, seen = [], [], set()
for line in open(pdb):
    if line.startswith("ATOM"):
        rid = int(line[22:26])
        if rid not in seen:
            seen.add(rid)
            resids.append(rid)
            bvals.append(float(line[60:66]))

plddt = np.array(bvals)
print(f"Residues: {len(plddt)}, mean pLDDT: {plddt.mean():.1f}, min: {plddt.min():.1f}, max: {plddt.max():.1f}")

fig, ax = plt.subplots(figsize=(10, 3), dpi=300)
colors = np.where(plddt > 90, "#0566e8",
         np.where(plddt > 80, "#4dc6ff",
         np.where(plddt > 70, "#f7e04a", "#fb9b06")))
ax.bar(resids, plddt, color=colors, width=1.0)
ax.axhline(70, ls="--", c="gray", lw=0.8)
ax.axhline(50, ls="--", c="gray", lw=0.8)
ax.set_xlabel("Residue")
ax.set_ylabel("pLDDT")
ax.set_ylim(0, 100)
plt.tight_layout()
plt.savefig(r"C:\Users\ACER\Desktop\Capripox_TCell_Epitope_Discovery\figures\Fig2_C_pLDDT_profile.png")
plt.close()
print("Saved Fig2_C_pLDDT_profile.png")