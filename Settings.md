# OpenFOAM Directory & Key Settings

Stock OpenCFD **OpenFOAM v2412**, solver **`rhoCentralFoam`**. No hyStrath /
hy2Foam component and no custom shared libraries.

Full working for every number below: `scripts/OpenFOAM configs/hypersonic_case/docs/PARAMETERS.md`.

## system/
*   **controlDict**
    *   Solver: `rhoCentralFoam`
    *   `deltaT`: `4.5e-08` (45 ns), **fixed** - `adjustTimeStep`: `no`
    *   Total steps: `675,000`; `endTime`: `0.030375` s (2.0 flow-throughs)
    *   CFL target `0.37`; actual `0.314` at the 0.5 mm mesh minimum
    *   Volume writes: every `33,750` steps, 20 over the run
    *   `fileHandler`: `collated` (no `reconstructPar` in the workflow)
    *   Function objects split into `monitorFunctions`, `probesPSE`, `surfacesShock`
*   **fvSchemes**
    *   `fluxScheme Kurganov` - Kurganov-Noelle-Petrova central-upwind (KNP)
    *   `vanLeer` / `vanLeerV` TVD limiters on every reconstruction; these are
        the only subgrid dissipation (ILES closure)
*   **fvSolution**
    *   `diagonal` on the conservative variables; `smoothSolver` on the
        implicit diffusion correction. Energy variable is `e`.
*   **decomposeParDict**
    *   `384` subdomains, `ptscotch`, about 81,000 cells per rank
*   **probesPSE**
    *   Second-Mack-mode wall probes at 10-step (450 ns) cadence, Nyquist
        1.111 MHz; boundary-layer rakes at 100 steps
*   **surfacesShock**
    *   Shock isosurfaces and cutting planes written every 1000 steps as
        binary `.vtp` for ParaView

## constant/
*   **thermophysicalProperties**
    *   `hePsiThermo` / `pureMixture` / `polynomial` / `janaf` /
        `perfectGas` / `sensibleInternalEnergy`
    *   Single frozen-composition air pseudo-species, `molWeight` 28.95999
    *   `janaf` valid **200-6000 K**, Tcommon 1000 K
    *   `polynomial` transport (degree 7) valid **200-6000 K**, fitted to
        Sutherland below 500 K and Blottner-Wilke above 1000 K, kappa by
        modified Eucken. Max fit error 1.06 percent on mu, 1.41 on kappa.
*   **turbulenceProperties**
    *   `simulationType`: `laminar`
    *   `laminar { model Stokes; }` - ILES, no explicit SGS model

## 0/
*   **p, U, T** - and nothing else.
    *   `p` 292.12 Pa, `U` (3169.7 0 0) m/s, `T` 250 K
    *   Atmosphere, Inlet and Outlet at **250 K**; `Solid_Walls` isothermal at
        **1500 K**
    *   The 11 species mass fractions and the `Tt`/`Tv_*` two-temperature
        fields are gone with the reacting two-temperature model.
