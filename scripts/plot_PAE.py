import json
import numpy as np
import matplotlib.pyplot as plt

jf = r"C:\Users\ACER\Downloads\test_16bca_0.result\test_16bca_0\test_16bca_0_scores_rank_001_alphafold2_ptm_model_2_seed_000.json"

d = json.load(open(jf))
pae = np.array(d["pae"])
print("PAE shape:", pae.shape)

plt.figure(figsize=(6, 5), dpi=300)
plt.imshow(pae, cmap="bwr", vmin=0, vmax=30)
plt.colorbar(label="PAE (?)")
plt.xlabel("Scored residue"); plt.ylabel("Scored residue")
plt.tight_layout()
plt.savefig(r"C:\Users\ACER\Desktop\Capripox_TCell_Epitope_Discovery\figures\Fig2_D_PAE.png")
plt.close()
print("Saved Fig2_D_PAE.png")
