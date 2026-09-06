# OpenFOAM Directory & Key Settings
## system/
*   **controlDict**
    *   Solver: `hy2Foam`
    *   `endTime`: `0.03` (30 ms)
    *   `adjustTimeStep`: `yes`
    *   `maxCo`: `0.4`
*   **fvSchemes**
    *   Kurganov-Noelle-Petrova scheme for convective flux discretization using TVD 2nd order limiters
*   **fvSolution**
    *   Linear matrix solver tolerances and preconditioning options

## constant/
*   **thermophysicalProperties**
    *   11-species reacting air model with Gupta transport properties
*   **chemistryProperties**
    *   Park / Dunn finite-rate chemical kinetic reaction mechanisms
*   **turbulenceProperties**
    *   `simulationType`: `laminar`
    *   `laminarModel`: `Stokes`

## 0/
*   **U, p, T, Tv<sub>(per element)</sub>, Te**
    *   Thermo-chemical initial and boundary conditions
*   **Ydefault...**
    *   Mass fractions for 11 species (N2, O2, NO, N, O, N2+, O2+, NO+, N+, O+, e-)
