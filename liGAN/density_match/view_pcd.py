import open3d as o3d
import numpy as np
import argparse

p = argparse.ArgumentParser(description="Visualize point cloud (PCD file) with Open3D")
p.add_argument("files", type=str, nargs="+", help="PCD files")

args = p.parse_args()

pcds = []

for f in args.files:
    pcd = o3d.io.read_point_cloud(f)

    # Multiply colors column by 10000
    # to get different colors in visualization
    # from RGB values as consecutive indices
    # pcd.colors = o3d.utility.Vector3dVector(np.asarray(pcd.colors) * 1000)

    pcds.append(pcd)

o3d.visualization.draw_geometries(pcds)
