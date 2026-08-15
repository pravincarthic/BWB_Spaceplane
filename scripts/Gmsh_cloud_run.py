#!/usr/bin/env python
###
### SALOME v9.16.0 - Validated Gmsh 2.2 Meshing Pipeline
###

import sys
import math
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()

###
### GEOM component
###
import GEOM
from salome.geom import geomBuilder
import SALOMEDS

geompy = geomBuilder.New()

# 1. Coordinate Base & Direction Vectors
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

# 2. Import and Orient Spaceplane CAD
BWB_Spaceplane_V1 = geompy.ImportSTEP("/home/user/Design_without_cavity_v1_20260806.step", False, True)
geompy.Rotate(BWB_Spaceplane_V1, OY, -5.0 * math.pi / 180.0)

# 3. Create Far-Field Domain & Fluid Subtraction
Cone_1 = geompy.MakeCone(O, OX, 13.000055, 39.000165, 53.0)
Cut_1 = geompy.MakeCutList(Cone_1, [BWB_Spaceplane_V1], True)

# 4. Construct Geometric Boundary Face Groups
Outlet = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Outlet, [10])

Inlet = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Inlet, [12])

Atmosphere = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Atmosphere, [3])

Solid_Walls = geompy.CreateGroup(Cut_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Solid_Walls, [15, 25, 32, 51, 56, 73, 80, 88, 91, 94, 97, 100, 103, 105, 108])

# 5. Register CAD Entities in Study
geompy.addToStudy(O, 'O')
geompy.addToStudy(OX, 'OX')
geompy.addToStudy(OY, 'OY')
geompy.addToStudy(OZ, 'OZ')
geompy.addToStudy(BWB_Spaceplane_V1, 'BWB Spaceplane V1')
geompy.addToStudy(Cone_1, 'Cone_1')
geompy.addToStudy(Cut_1, 'Cut_1')
geompy.addToStudyInFather(Cut_1, Outlet, 'Outlet')
geompy.addToStudyInFather(Cut_1, Inlet, 'Inlet')
geompy.addToStudyInFather(Cut_1, Atmosphere, 'Atmosphere')
geompy.addToStudyInFather(Cut_1, Solid_Walls, 'Solid_Walls')

###
### SMESH component
###
import SMESH
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Cut_1, 'Mesh_1')

# 1. Initialize Gmsh Hypothesis
GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
Gmsh_Parameters = GMSH.Parameters()

# 2. Numerical Discretization Algorithms
Gmsh_Parameters.Set2DAlgo(6)              # Frontal-Delaunay
Gmsh_Parameters.Set3DAlgo(10)             # Parallel Delaunay (HXT)
Gmsh_Parameters.SetIs2d(0)                # 3D Domain Flag
Gmsh_Parameters.SetMinSize(0.0008)        # 0.8 mm Surface Resolution
Gmsh_Parameters.SetMaxSize(0.15)           # 500 mm Farfield Limit
Gmsh_Parameters.SetMeshCurvatureSize(20)  # 20 Elements per 2*pi Arc
Gmsh_Parameters.SetSmouthSteps(5)

# 3. Viscous Layer Inflation (S = 0.042m, N = 40, r = 1.18 -> y1 ~ 10 um)
Viscous_Layers_1 = GMSH.ViscousLayers(
    0.04965, 45, 1.16,
    [15, 25, 32, 51, 56, 73, 80, 88, 91, 94, 97, 100, 103, 105, 108],
    0, smeshBuilder.SURF_OFFSET_SMOOTH
)

# 4. Map Boundary Groups to Mesh
Outlet_1 = Mesh_1.GroupOnGeom(Outlet, 'Outlet', SMESH.FACE)
Inlet_1 = Mesh_1.GroupOnGeom(Inlet, 'Inlet', SMESH.FACE)
Atmosphere_1 = Mesh_1.GroupOnGeom(Atmosphere, 'Atmosphere', SMESH.FACE)
Solid_Walls_1 = Mesh_1.GroupOnGeom(Solid_Walls, 'Solid_Walls', SMESH.FACE)

# 5. Compute Volume Grid
print("Computing 3D Mesh via Gmsh HXT...")
isDone = Mesh_1.Compute()

# 6. Export to Gmsh 2.2 Format
if isDone:
    print("Mesh generation succeeded. Exporting to Gmsh 2.2 format...")
    try:
        Mesh_1.ExportGMSHIO(r'/home/user/Gmsh_without_cavities_AoA_neg5_20260811.msh', 'Gmsh 2.2', Mesh_1)
        print("Export completed successfully.")
    except Exception as err:
        print(f"ExportGMSHIO failed: {err}")
else:
    print("Mesh computation failed. Check face IDs and surface curvature limits.")

# 7. Publish to SALOME Object Tree
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Viscous_Layers_1, 'Viscous Layers_1')
smesh.SetName(Outlet_1, 'Outlet')
smesh.SetName(Atmosphere_1, 'Atmosphere')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(GMSH.GetAlgorithm(), 'GMSH')
smesh.SetName(Solid_Walls_1, 'Solid_Walls')

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()