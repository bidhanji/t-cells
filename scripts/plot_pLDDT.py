import json
import numpy as np
import matplotlib.pyplot as plt

jf = r"C:\Users\ACER\Downloads\test_16bca_0.result\test_16bca_0\test_16bca_0_scores_rank_001_alphafold2_ptm_model_2_seed_000.json"

d = json.load(open(jf))
plddt = np.array(d["plddt"])
print(f"Residues: {len(plddt)}, mean pLDDT: {plddt.mean():.1f}, min: {plddt.min():.1f}, max: {plddt.max():.1f}")

fig, ax = plt.subplots(figsize=(10, 3), dpi=300)
colors = np.where(plddt > 90, "#0566e8",
         np.where(plddt > 80, "#4dc6ff",
         np.where(plddt > 70, "#f7e04a", "#fb9b06")))
ax.bar(range(1, len(plddt)+1), plddt, color=colors, width=1.0)
ax.axhline(70, ls="--", c="gray", lw=0.8)
ax.axhline(50, ls="--", c="gray", lw=0.8)
ax.set_xlabel("Residue"); ax.set_ylabel("pLDDT"); ax.set_ylim(0, 100)
plt.tight_layout()
plt.savefig(r"C:\Users\ACER\Desktop\Capripox_TCell_Epitope_Discovery\figures\Fig2_C_pLDDT_profile.png")
plt.close()
print("Saved Fig2_C_pLDDT_profile.png")
