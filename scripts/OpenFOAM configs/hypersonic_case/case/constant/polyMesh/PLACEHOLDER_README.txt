THIS DIRECTORY IS EMPTY ON PURPOSE.

The mesh could not be uploaded to generate this case (upload size limits), so
the polyMesh/ directory (points, faces, owner, neighbour, boundary) has NOT
been created here.

To populate it, edit the MESH_FILE placeholder in scripts/setup.sh to point
to your actual .msh file, then run:

    bash scripts/setup.sh

setup.sh calls gmshToFoam on your mesh and writes points/faces/owner/
neighbour/boundary into this directory automatically. Once that succeeds,
delete this placeholder file - its presence does not affect the solver, it
is just a reminder.

IMPORTANT: after conversion, open constant/polyMesh/boundary and confirm the
actual patch names gmshToFoam produced. The boundary condition files in 0/
and the wall-normal-distance/refinement notes in docs/ assume patch names:

    farfield   (freestream / far-field boundary)
    wall       (vehicle surface, no-slip)
    outlet     (downstream outflow, if present as a separate patch)

If your Gmsh physical group names differ (e.g. "inlet", "body", "symm1"),
either rename them in Gmsh before conversion, or rename/copy the boundary
condition entries in case/0/* to match.
