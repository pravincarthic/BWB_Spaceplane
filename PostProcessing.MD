# Post Processing - What to do

OpenFOAM v1706 will produce four main types of output data.
Because we are running an uncollated parallel run across 384 processors, the files will be distributed between separate core processor directories and automated global output directories.
Here is exactly what files and folders will appear in case/ root directory after execution:
------------------------------
## 1. The Raw Decomposed Fields (processor0 to processor383)
OpenFOAM will populate each of 384 processor folders with time-step subdirectories matching core simulation outputs (e.g., folder names like 0.01, 0.02, etc.). Inside every time folder, one will see fractured, uncollated pieces of primary physical fields:

* p / T / U: Localised arrays for pressure, temperature, and velocity vectors.
* Tt / Tv: Requested total temperature and vibrational/electronic temperature fields tracked by hypersonic thermophysical solver.

------------------------------
## 2. 3D Shock Wave Snapshots (postProcessing/shockWave3D)
Even though the run is uncollated, OpenFOAM will automatically stitch 3D view outputs together into a single global folder. Inside postProcessing/shockWave3D/, there will be exactly 1,000 subdirectories (named by the exact time-step number, such as 300, 600, 900, etc.).

* bowShockIso.vtk: Inside each time directory, this is the lightweight file. This is the complete, 3D rotatable geometric surface shell capturing where the pressure hits target 1000 Pa limit. One can load this file sequence directly into ParaView on PC to choose animation perspectives.

------------------------------
## 3. Aerodynamic Performance Logs (postProcessing/forces and forceCoeffs)
Inside these two folders, OpenFOAM will generate plain-text tracking files formatted cleanly as tab-separated data matrices.

* forces.dat: Contains an active timeline tracking X, Y, and Z total forces divided cleanly into pressure forces and viscous (skin friction) forces acting across wall boundary patch.
* forceCoeffs.dat: Tracks aerodynamic dimensionless performance numbers over time. It will print out active $C_d$ (Drag Coefficient aligned to the X-axis) and $C_l$ (Lift Coefficient aligned to the Z-axis) values based on reference inputs. We can import this small file straight into MATLAB, Python, or Excel to plot lift-to-drag metrics.

------------------------------
## 4. Convergence Limits Log (postProcessing/fieldMinMax)

* fieldMinMax.dat: This text file records a single line per printed time-step, tracking the absolute lowest and highest values of pressure (p), temperature (T), total temperature (Tt), and vibrational temperature (Tv) across the entire cell domain.
* Why this is critical for Mach 10 simulations: If simulation is going to explode numerically, we can see the maximum values in this log file spike to impractical levels (e.g., temperatures climbing past 50,000K or negative pressures appearing) right before the code throws a fatal crash error.


Map output files into the specific sections of project report:
------------------------------
## 📊 1. Aerodynamic Performance Metrics (The Core Results Section)
Reviewers will look for global vehicle forces to assess the vehicle's aerodynamic efficiency at a $-5^\circ$ Angle of Attack.

* The Source File: postProcessing/forceCoeffs/0/forceCoeffs.dat
* What to extract: Lift Coefficient ($C_l$), Drag Coefficient ($C_d$), and the Lift-to-Drag ratio ($L/D$).
* How to present it:
* Plot of Time vs. $C_l,$ $C_d$ line chart to prove simulation reached a steady-state converged value.
   * Include a summary data table comparing final $C_l$ and $C_d$ against theoretical or historical hypersonic baseline data (e.g., modified Newtonian impact theory).

------------------------------
## 🔎 2. Skin Friction vs. Wave Drag Analysis (Aerodynamic Drag Breakdown)
At Mach 10, total drag is a fierce battle between pressure drag (shock wave formation) and viscous drag (skin friction from the hypersonic boundary layer).

* The Source File: postProcessing/forces/0/forces.dat
* What to extract: The separation of forces into pressure and viscous columns for $X$-axis (drag direction).
* How to present it:
* A Pie Chart or Stacked Bar Chart showing the percentage breakdown of total drag (e.g., "65% Wave Drag vs. 35% Viscous Drag"). This provides excellent physical insight into vehicle's hypersonic shape optimization.

------------------------------
## 🌟 3. The 3D Flow Topology & Shock Structure (Main Visual Figures)
This provides the visual evidence of $-5^\circ$ AoA flight profile, capturing the asymmetry of the flow field.

* The Source File: postProcessing/shockWave3D/<timeStep>/bowShockIso.vtk
* How to present it:
* Figure A (3D Isometric View): Rendering the full 3D 1000 Pa pressure envelope wrapping around vehicle. This visually demonstrates the expansion of the shock layer over the top of the craft and the compression beneath it.
   * Figure B (Centerline 2D Slice): Cutting right down the center plane ($Z$-normal) inside ParaView. Coloring this slice using your Temperature ($T$) field. This highlights the crisp stand-off distance of bow shock, the hot stagnation points, and the expansion fan trailing off the vehicle's aft region.

------------------------------
## 📊 4. Thermal & Chemical Equilibrium Checks (Hypersonic Physics Validation)
Because we are running a hypersonic solver tracking $Tt$ (Total/Translational Temperature) and $Tv$ (Vibrational/Electronic Temperature), we must present that the gas model correctly captures high-temperature, non-equilibrium effects.

* The Source File: postProcessing/fieldMinMax/0/fieldMinMax.dat
* What to extract: The maximum values of $T$, $Tt$, and $Tv$ at final steady-state time steps.
* How to present it in the paper:
* Thermal Non-Equilibrium Validation: Plot a line profile or a comparison showing where $Tt$ and $Tv$ diverge (behind the strong bow shock) and where they equilibrate. This proves the multi-temperature thermodynamic solver settings are physically accurate for a 130,000 ft altitude atmosphere.

------------------------------
## 🛡️ 5. Mesh Independence Study (The Methodology / Verification Section)
Prove 31 million was the correct size to choose (and not just an arbitrary guess).

* The Source File: log.checkMesh.log (generated in parallel) and final forceCoeffs.dat.
* How to present it:
* Include checkMesh metrics (max skewness, max non-orthogonality) in a table to prove grid quality.
   * Create a standard Mesh Convergence Plot: Plot total grid cell count on the $X$-axis (e.g., 10M, 20M, 31M cells) against Lift Coefficient ($C_l$) on the $Y$-axis. Show that 31M cell grid yields the same $C_l$ value as a finer grid, proving that spatial resolution has fully resolved the Mach 10 shock layer.

