import molgrid
import torch
import argparse
import os

from openbabel import pybel

p = argparse.ArgumentParser(description="Visualize point cloud (PCD file) with Open3D")
p.add_argument("sdf", type=str, help="SDF file")
p.add_argument("-r", "--resolution", type=float, default=0.5, help="Grid resolution")
p.add_argument(
    "-d", "--dimension", type=float, default=23.5, help="Grid dimension (1D)"
)
p.add_argument(
    "-m", "--ligmap", type=str, default="files/ligmap", help="Ligand types file"
)
p.add_argument("-o", "--output", type=str, default=None, help="Output file")
p.add_argument("--dx", action="store_true", help="Output grids as DX files")

args = p.parse_args()

system = os.path.splitext(os.path.basename(args.sdf))[0]

if args.output is None:
    args.output = f"{system}.pcd"

resolution = args.resolution
dimension = args.dimension

gm = molgrid.GridMaker(resolution=resolution, dimension=dimension)

t = molgrid.FileMappedGninaTyper(args.ligmap)

# Grid dimensions (including types)
gdims = gm.grid_dimensions(t.num_types())

# Pre-allocate grid
# Only one example (batch size is 1)
grid = torch.zeros(1, *gdims, dtype=torch.float32, device="cuda:0")

obmol = next(pybel.readfile("sdf", args.sdf))
obmol.addh()
print(obmol, end="")

# Use OpenBabel molecule object (obmol.OBmol) instead of PyBel molecule (obmol)
cs = molgrid.CoordinateSet(obmol.OBMol, t)

ex = molgrid.Example()
ex.coord_sets.append(cs)

c = ex.coord_sets[0].center()  # Only one coordinate set
print("center:", tuple(c))

# https://gnina.github.io/libmolgrid/python/index.html#the-transform-class
transform = molgrid.Transform(
    c, random_translate=0.0, random_rotation=False,  # float  # bool
)
transform.forward(ex, ex)

# Compute grid
gm.forward(ex, grid[0])

print("grid.shape:", grid.shape)

if args.dx:
    # https://gnina.github.io/libmolgrid/python/index.html#molgrid.write_dx_grids
    # Grid4f is different from Grid4fCUDA
    # If a function takes Grid4f as input, torch.Tensor need to be moved to the CPU
    molgrid.write_dx_grids(
        f"grids/{system}",
        t.get_type_names(),
        grid[0].cpu(),
        ex.coord_sets[0].center(),
        0.5,
    )

cloud = torch.logical_and(grid[0] >= 0.4, grid[0] <= 0.6)


def grid_lims(c):
    """
    Compute grid limits given coordinate of center.
    """
    return c - dimension / 2.0, c + dimension / 2.0


# Number of voxels in one direction
n_steps = round(dimension / resolution) + 1

x1, x2 = grid_lims(c[0])
x = torch.linspace(x1, x2, steps=n_steps)

y1, y2 = grid_lims(c[1])
y = torch.linspace(y1, y2, steps=n_steps)

z1, z2 = grid_lims(c[2])
z = torch.linspace(z1, z2, steps=n_steps)

X, Y, Z = torch.meshgrid(x, y, z)

# Total number of points in the cloud
n_points = torch.sum(cloud)

with open(args.output, "w") as fout:

    fout.write(
        f"""VERSION .7
    FIELDS x y z rgb
    SIZE 8 8 8 4
    TYPE F F F F
    COUNT 1 1 1 1
    WIDTH {n_points}
    HEIGHT 1
    VIEWPOINT 0 0 0 1 0 0 0
    POINTS {n_points}
    DATA ascii\n"""
    )

    for i, name in enumerate(t.get_type_names()):
        idx = cloud[i]  # Indices of cloud points

        for x, y, z in zip(
            X[idx].cpu().numpy(), Y[idx].cpu().numpy(), Z[idx].cpu().numpy()
        ):
            # Miltiplication if to differentiate colors more
            fout.write(f"{x:.5f} {y:.5f} {z:.5f} {i * 1000}\n")
