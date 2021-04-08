# RDKit to Open Babel Conversion

## Problem

It is possible to obtain an atomic density grid from an Open Babel molecule using [libmolgrid](https://github.com/gnina/libmolgrid) as follows:

```python
# Ligand GNINA typer from file
lig_map = molgrid.FileMappedGninaTyper(ligmapfile)

# Create openbabel molecule
m = pybel.readstring('smi', 'c1c(Cl)cccc1CO')
m.addh()

# Create coodinate set
c = molgrid.CoordinateSet(m, lig_map)

# Create grid
ex = molgrid.Example()
ex.coord_sets.append(c)
for i in range(batch_size):
    grid_maker.forward(ex, grid_true[i])
```

However, it appears that the grid obtained loading `molecules/molecule.sdf` with Open Babel is different from the grid obtained by loading `molecules/molecule.sdf` using RDKit and convert it to an Open Babel molecule using [liGAN](https://github.com/mattragoza/liGAN)'s `molecule.rd_mol_to_ob_mol()`. This causes the generated molecules to be dependent on the library parsing SDF files.

## Solution

The problem is that a bond with a RDKit type `Chem.BondType.AROMATIC` is converted to a bond order of `1` in Open Babel (with an additional flag for aromaticity), while the original molecule was Kekulized. Kekulizing the RDKit molecule while removing aromaticity flags before the conversion seems to solve the problem:

```python
Chem.Kekulize(rd_mol, clearAromaticFlags=True)
```
