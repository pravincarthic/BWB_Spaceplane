import gmsh
import numpy as np

# User input prompts
input_file = input("Enter the path to the Gmsh .msh file: ")
output_file = input_file.replace('.msh', '.su2')

# Initialize and load
gmsh.initialize()
gmsh.open(input_file)

# 1. Pull absolute list of all node tags across all dimensions (1D, 2D, and 3D)
node_tags, coords, _ = gmsh.model.mesh.getNodes()

# 2. Bulk multiply the raw coordinates uniformly (e.g., mm -> m)
scaled_coords = np.array(coords) * 0.001
print(f"Scaling all nodes by a factor of 0.001 (mm to m). Total nodes: {len(node_tags)}")

# 3. Update the positions of all individual nodes in Gmsh's database
for tag, coord in zip(node_tags, scaled_coords.reshape(-1, 3)):
    gmsh.model.mesh.setNode(tag, coord.tolist(), [])

gmsh.model.geo.synchronize()
gmsh.write(output_file)
gmsh.finalize()

print(f"Successfully scaled all 1D, 2D, and 3D nodes. Converted {input_file} to {output_file}")