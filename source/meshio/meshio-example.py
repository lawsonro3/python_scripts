import meshio as mio

# Read mesh
mesh = mio.readVtp('./NonSymmetrical.dat')

# Write meshes
#mesh.writeVtp('./NonSymmetrical-demoOutput.vtp')
#mesh.writeNemohMesh('./NonSymmetrical-demoOutput.dat')
#mesh.writeGdf('./NonSymmetrical-demoOutput.gdf')