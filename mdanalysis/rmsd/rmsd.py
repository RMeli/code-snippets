import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

import numpy as np

from spyrmsd import rmsd

from MDAnalysis.tests.datafiles import PSF, DCD
from matplotlib import pyplot as plt

verbose = True
symmetry = True
step = 2

top_file = "5JWT.tpr"
traj_file = "5JWT.xtc"

iframe = 0
select = "resname LIG"

traj = mda.Universe(top_file, traj_file)
ref = mda.Universe(top_file, traj_file)

traj.trajectory[-1]

# Set reference frame
ref.trajectory[iframe]

# Align trajectory to reference
# _ = align.AlignTraj(traj, ref, select='name CA').run()

# Select reference atoms
# TODO: Make sure this does not change?
ref_selection = ref.select_atoms(select)
ref_positions = ref_selection.atoms.positions

# Select atoms for RMSD calculation
selection = traj.select_atoms(select)
n_atoms = len(selection)

if verbose:
    print(f"Number of atoms in '{select}': {n_atoms}")
    print(ref_positions.shape)

# Build adjacency matrix for selection
# FIXME: Ensure no bonds go outside selection
try:
    A = np.zeros((n_atoms, n_atoms))

    # FIXME: More efficient way of building the adjacency matrix
    for bond in selection.bonds:
        for i, ai in enumerate(selection.atoms):
            for j, aj in enumerate(selection.atoms):
                if ai in bond and aj in bond and i != j:
                    A[i, j] = 1
    if verbose:
        print(f"Adjacency matrix for '{select}'\n{A}")

except mda.exceptions.NoDataError:
    print(f"No bond information available for '{select}'")
    exit(-1)

# FIXME: Hard-coded atomic numbers for benzene
# FIXME: Need a consistent way to get them from selection!
atomicnums = np.array([6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1])

time = []
rmsds = []
rmsd_nosymm = []
for ts in traj.trajectory[::step]:
    # Align protein backbones
    old_rmsd, new_rmsd = align.alignto(traj, ref, select="backbone", strict=True)
    assert new_rmsd <= old_rmsd

    if symmetry:
        r = rmsd.symmrmsd(
            selection.atoms.positions, ref_positions, atomicnums, atomicnums, A, A
        )
    else:
        r = rmsd.rmsd(selection.atoms.positions, ref_positions, atomicnums, atomicnums)

    time.append(ts.time)
    rmsds.append(r)

    rmsd_nosymm.append(
        rmsd.rmsd(selection.atoms.positions, ref_positions, atomicnums, atomicnums)
    )

# Check rmsd_nosymm is the same as rms.RMDS
rmsd_mda = rms.RMSD(
    traj, ref, select="backbone", groupselections=["resname LIG"], ref_frame=0
)
R = rmsd_mda.run()
assert np.allclose(rmsd_nosymm, R.results.rmsd[::step, -1])

plt.plot(time, rmsds, label="symmrmsd")
plt.plot(time, rmsd_nosymm, label="rmsd")
plt.plot(time, R.results.rmsd[::step, -1], label="mda")
plt.xlabel("Time")
plt.ylabel("RMSD")
plt.legend()
plt.savefig("rmsd.png")
