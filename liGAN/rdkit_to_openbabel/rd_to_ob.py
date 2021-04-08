from openbabel import openbabel as ob
from openbabel import pybel

from rdkit import Chem


def rd_mol_to_ob_mol(rd_mol):
    """
    liGAN conversion between RDKit and Open Babel

    https://github.com/mattragoza/liGAN
    """
    ob_mol = ob.OBMol()
    rd_conf = rd_mol.GetConformer(0)

    for i, rd_atom in enumerate(rd_mol.GetAtoms()):
        atomic_num = rd_atom.GetAtomicNum()
        rd_coords = rd_conf.GetAtomPosition(i)
        x = rd_coords.x
        y = rd_coords.y
        z = rd_coords.z
        ob_atom = ob_mol.NewAtom()
        ob_atom.SetAtomicNum(atomic_num)
        ob_atom.SetAromatic(rd_atom.GetIsAromatic())
        ob_atom.SetVector(x, y, z)

    for k, rd_bond in enumerate(rd_mol.GetBonds()):
        i = rd_bond.GetBeginAtomIdx() + 1
        j = rd_bond.GetEndAtomIdx() + 1
        bond_type = rd_bond.GetBondType()
        bond_flags = 0
        if bond_type == Chem.BondType.AROMATIC:
            bond_order = 1
            bond_flags |= ob.OB_AROMATIC_BOND
        elif bond_type == Chem.BondType.SINGLE:
            bond_order = 1
        elif bond_type == Chem.BondType.DOUBLE:
            bond_order = 2
        elif bond_type == Chem.BondType.TRIPLE:
            bond_order = 3
        ob_mol.AddBond(i, j, bond_order, bond_flags)

    return ob_mol


molfile = "molecules/molecule.sdf"

# I/O with Open Babel
ob_mol = next(pybel.readfile("sdf", molfile))
ob_mol.write("sdf", "molecules/ob_mol.sdf", overwrite=True)

# I/O with RDKit
rd_mol = next(Chem.SDMolSupplier(molfile, removeHs=False))
rd_writer = Chem.SDWriter("molecules/rd_mol.sdf")
rd_writer.write(rd_mol)
rd_writer.close()

# I with RDKit / O with Open Babel
rd_mol = next(Chem.SDMolSupplier(molfile, removeHs=False))
ob_mol = rd_mol_to_ob_mol(rd_mol)
ob_mol = pybel.Molecule(ob_mol)
ob_mol.write("sdf", "molecules/converted_nokek.sdf", overwrite=True)

# I with RDKit / O with Open Babel
# Kekulize RDKit molecule and remove aromatic flags
rd_mol = next(Chem.SDMolSupplier(molfile, removeHs=False))
Chem.Kekulize(rd_mol, clearAromaticFlags=True)
ob_mol = rd_mol_to_ob_mol(rd_mol)
ob_mol = pybel.Molecule(ob_mol)
ob_mol.write("sdf", "molecules/converted_kek.sdf", overwrite=True)
