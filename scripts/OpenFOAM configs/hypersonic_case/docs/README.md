# BWB spaceplane - rhoCentralFoam case, OpenFOAM v2412

Mach 9.99 at 130,000 ft over a 48 m blended-wing-body spaceplane at -5 deg
angle of attack. ILES on a 31M cell wall-resolved hybrid mesh, 384 cores,
instrumented for second-Mack-mode PSE analysis and for shock-wave
visualisation in ParaView.

Full numbers and their derivations: `docs/PARAMETERS.md`.
Design rationale: `docs/SETTINGS.md`.
What to check before trusting a result: `docs/VALIDATION.md`.

## What this case is

| | |
|---|---|
| Solver | `rhoCentralFoam`, stock OpenFOAM v2412 |
| Thermo | `janaf`, Thigh 6000 K |
| Transport | `polynomial`, degree 7, valid to 6000 K |
| Energy | `sensibleInternalEnergy` |
| Flux scheme | `Kurganov` (KNP) with `vanLeer` limiters |
| Turbulence | laminar / Stokes - ILES, no explicit SGS model |
| Timestep | fixed 45 ns, 675,000 steps, CFL 0.37 target |
| Cores | 384, `ptscotch`, collated output |

No hyStrath or hy2Foam component is used. Everything is core OpenFOAM.

## Layout

```
case/
  0/                     p, U, T                (that is all - see below)
  constant/
    thermophysicalProperties   janaf + polynomial air, 200-6000 K
    turbulenceProperties       laminar / Stokes (ILES)
    polyMesh/                  written by scripts/setup.sh
  system/
    controlDict          45 ns fixed step, 675k steps, collated, 20 writes
    fvSchemes            Kurganov + vanLeer
    fvSolution           diagonal + smoothSolver diffusion correction
    decomposeParDict     384 ranks, ptscotch
    monitorFunctions     derived fields, Courant, forces, heat flux, averages
    probesPSE            second-mode wall probes and BL rakes
    surfacesShock        shock isosurfaces and cutting planes for ParaView
  case.foam              open this in ParaView
scripts/
  setup.sh               gmshToFoam, checkMesh, dictionary validation
  make_pse_probes.py     regenerate probesPSE from the real wall surface
  run.sh                 decomposePar + mpirun
  submit_384.slurm       batch submission
  postProcess.sh         what to open and what to check
  cleanup.sh             archive a finished run
docs/
```

`0/` holds three fields. The previous hy2Foam configuration carried 24 - two
temperatures plus 11 species plus a vibrational temperature per species - and
all of those are gone with the two-temperature model.

## Workflow

```bash
source /opt/openfoam2412/etc/bashrc        # site-specific

export MESH_FILE=/path/to/your.msh
bash scripts/setup.sh                      # convert, check, validate

python3 scripts/make_pse_probes.py         # exact probe placement (recommended)

bash scripts/run.sh 384                    # or: sbatch scripts/submit_384.slurm
bash scripts/postProcess.sh
```

`scripts/setup.sh` prints the mesh patch names. They must be `Inlet`,
`Outlet`, `Atmosphere` and `Solid_Walls` - those names appear in `0/*` and in
every function object.

## Three things to check before the production run

**1. The minimum cell size.** The timestep is fixed, so nothing adapts if the
Courant number climbs. CFL 0.37 at 45 ns requires the smallest cell in the
domain to be at least 0.424 mm. The mesh generator floor is 0.5 mm, giving
Co = 0.314. `scripts/setup.sh` prints this check against the `checkMesh`
output; if any cell is smaller, reduce `deltaT` in `system/controlDict`.

**2. The mesh units.** The `.msh` from `scripts/final_meshing.py` is in
**millimetres** - it scales the STEP geometry by 1000 so the Gmsh sizing
constants can be written in mm. The OpenFOAM case is in metres, so the
conversion must run as `gmshToFoam -scale 0.001`, which `scripts/setup.sh`
does by default. Get this wrong and the solver runs without complaint on a
48 km vehicle. `setup.sh` step 3a and `make_pse_probes.py` both abort if the
converted mesh is not metre-scale.

**3. The probe locations.** The coordinates shipped in `system/probesPSE`
are derived from the vehicle bounding box and rely on `patchProbes` snapping
to the nearest wall face. Run `scripts/make_pse_probes.py` once the mesh
exists to replace them with exact on-surface points and true wall-normal
rakes.

## Output

| Where | What | Cadence |
|---|---|---|
| time directories | full volume fields | 20 over the run |
| `postProcessing/shockSurfaces/` | shock isosurfaces, cutting planes (`.vtp`) | every 1000 steps |
| `postProcessing/wallSurface/` | surface p, heat flux, shear (`.vtp`) | every 1000 steps |
| `postProcessing/pseWall*/` | wall p, T, U, rho time series | every 10 steps |
| `postProcessing/pseBLRakes/` | wall-normal profiles | every 100 steps |
| `postProcessing/force*/` | loads and coefficients | every 100 steps |
| `postProcessing/courantMonitor/` | max Courant number | every 100 steps |

Open `case/case.foam` in ParaView for the volume fields; open the `.vtp`
series directly for the shock animation. Nothing needs `reconstructPar`.

## What changed from the previous configuration

This case previously targeted `hy2Foam` from the hyStrath suite: an
11-species reacting air model with a two-temperature (translational plus
vibrational) formulation, Blottner-Eucken transport, Park 1993 chemistry, and
three custom shared libraries.

It now targets stock `rhoCentralFoam`. Files removed:
`constant/transportProperties`, `constant/chemistryProperties`,
`constant/thermoDEM`, `constant/thermo2TModel`, `constant/hTCProperties`,
`constant/chemDicts/`, and 22 of the 24 fields in `0/`.

The physical justification for dropping the reacting two-temperature model is
in `docs/SETTINGS.md`: the peak temperature this case actually reaches is
4331 K, below the point where air dissociation materially changes the
aerodynamics, and a frozen-composition calorically imperfect gas captures the
cp(T) and mu(T) variation that does matter.
