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
* **Environment:** OpenFOAM (strictly v1706) with `hyStrath` / `hy2Foam` extension compiled.
* **CAD / Meshing:** `SALOME NETGEN 1D-2D-3D` or `SALOME GMSH` with high-refinement surface feature edges.
* **HPC Environment:** SLURM workload manager with Intel MPI or OpenMPI.
## Quick Start: Running a hy2Foam Simulation
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
   decomposePar
   mpirun -np $NPROC hy2Foam -parallel > log.hy2Foam 2>&1 &
6. **Post-Processing:**
   Reconstruct case and launch ParaView (in case of ParaView GUI, create an empty `.foam` file and open it.):
   ```bash
   reconstructPar
   paraFoam
