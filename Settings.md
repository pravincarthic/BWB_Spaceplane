# OpenFOAM Directory & Key Settings
## system/
*   **controlDict**
    *   Solver: `hy2Foam`
    *   `endTime`: `0.03` (30 ms)
    *   `adjustTimeStep`: `yes`
    *   `maxCo`: `0.4`
*   **fvSchemes**
    *   High-order WENO/Kurganov schemes for convective flux discretization
*   **fvSolution**
    *   Linear matrix solver tolerances and preconditioning options

## constant/
*   **thermophysicalProperties**
    *   11-species reacting air model with Gupta transport properties
*   **chemistryProperties**
    *   Park / Dunn finite-rate chemical kinetic reaction mechanisms
*   **turbulenceProperties**
    *   `SimulationType`: `LES`
    *   `LESModel`: `WALE` or implicit LES (ILES)

## 0/
*   **U, p, T, Tvib, Te**
    *   Thermo-chemical initial and boundary conditions
*   **Ydefault...**
    *   Mass fractions for 11 species (N2, O2, NO, N, O, N2+, O2+, NO+, N+, O+, e-)
