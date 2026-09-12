# Theoretical Validation & Governing Equations

## 1. Mack Second-Mode Instability Damping
The acoustic impedance $Z_w = p'/v_w'$ applied at the wall boundary conditions models the porous or micro-cavity metasurface without requiring full spatial discretization of billions of internal cavity cells:
$$v_w' = A \cdot p_w'$$
Where $A$ represents the complex acoustic admittance of the surface, tuned to absorb frequencies in the $100\text{ kHz} \text{ to } 1\text{ MHz}$ band.

## 2. Lift-to-Drag ($L/D$) Optimization
By maintaining laminar flow over a larger fraction of the windward centerbody, skin-friction coefficient $C_f$ is significantly reduced compared to fully turbulent assumptions:
$$C_f \propto \frac{1}{\sqrt{Re_x}} \quad \text{(Laminar)} \quad \text{vs.} \quad C_f \propto \frac{1}{Re_x^{0.2}} \quad \text{(Turbulent)}$$
This preserves a higher effective $L/D$ ratio across the reentry corridor, enabling shallow, low-g atmospheric descent paths.

---

# Getting Started

## Prerequisites
* **Environment:** OpenFOAM **v2412** (stock OpenCFD). No `hyStrath` / `hy2Foam` build is needed - `rhoCentralFoam` ships with OpenFOAM.
* **CAD / Meshing:** `SALOME NETGEN 1D-2D-3D` or `SALOME GMSH` with high-refinement surface feature edges.
* **HPC Environment:** SLURM workload manager with Intel MPI or OpenMPI.
## Quick Start: Running a rhoCentralFoam Simulation
1. **Clone the repository:**
   ```bash
   git clone [https://github.com/pravincarthic/BWB_Spaceplane/](https://github.com/pravincarthic/BWB_Spaceplane/)
2. **Convert Mesh to OpenFOAM Format:**
   If using SALOME NETGEN `.unv` mesh:
   ```bash
   ideasUnvToFoam ../../../mesh/salome_netgen_mesh.unv
3. **Convert Mesh to OpenFOAM Format:**
   Or if using Gmsh `.msh` mesh:
   ```bash
   ideasUnvToFoam ../../../mesh/salome_netgen_mesh.unv
4. **Verify Mesh Quality:**
   ```bash
   checkMesh
5. **Execute Solver (Parallel):**
   ```bash
   decomposePar -fileHandler collated
   mpirun -np 384 rhoCentralFoam -parallel -fileHandler collated > log.rhoCentralFoam 2>&1 &
6. **Post-Processing:**
   With collated output there is no `reconstructPar` step - ParaView reads the
   decomposed case directly. Open the `case.foam` file that ships with the case
   and select "Decomposed Case" in the reader panel:
   ```bash
   paraview case/case.foam
