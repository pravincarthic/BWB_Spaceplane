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
