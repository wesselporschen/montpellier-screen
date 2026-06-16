from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
import numpy as np

# load smiles
mols = []
names = []

with open("./montpellier_ligands.smi") as f:
    for line in f:
        parts = line.strip().split()

        if not parts:
            continue

        smi = parts[0]
        name = parts[1] if len(parts) > 1 else smi

        mol = Chem.MolFromSmiles(smi)

        if mol:
            mols.append(mol)
            names.append(name)

print(f"Loaded {len(mols)} molecules")

# Morgan fingerprints (ECFP4 equivalent)
fpgen = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048
)

fps = [fpgen.GetFingerprint(m) for m in mols]

# pairwise tanimoto matrix
n = len(fps)
sim_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i, n):
        sim = TanimotoSimilarity(fps[i], fps[j])
        sim_matrix[i, j] = sim
        sim_matrix[j, i] = sim

# summary stats excluding diagonal
vals = sim_matrix[np.triu_indices(n, k=1)]

print(f"Mean similarity:   {vals.mean():.3f}")
print(f"Median similarity: {np.median(vals):.3f}")
print(f"Max similarity:    {vals.max():.3f}")
