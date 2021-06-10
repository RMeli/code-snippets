from rdkit import Chem
from rdkit.Chem import Draw

smiles = ["O=C=O", "C=C", "C1NC1", "C1CNC1", "c1cNcc1"]
mols = [Chem.MolFromSmiles(smi) for smi in smiles]
img = Draw.MolsToGridImage(mols, molsPerRow=3)
img.save("images/smiles_to_image.png")
