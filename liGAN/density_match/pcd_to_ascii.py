import open3d as o3d
import os
import argparse

p = argparse.ArgumentParser(
    description="Convert point cloud file to ASCII")
p.add_argument("file", type=str, help="file to convert to ASCII")
p.add_argument("-o", "--output", type=str, default=None, help="output file")

args = p.parse_args()

if args.output is None:
    fname, ext = os.path.splitext(args.file)

    args.output = f"{fname}_ascii{ext}"

pcd = o3d.io.read_point_cloud(args.file)
o3d.io.write_point_cloud(args.output, pcd, write_ascii=True)
