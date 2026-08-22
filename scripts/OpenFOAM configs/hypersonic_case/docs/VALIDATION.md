# Validation Checklist

## 1. Mesh conversion (before running the solver)

- `gmshToFoam` completed without errors (check `scripts/setup.sh` output)
- `checkMesh -allTopology -allGeometry` reports no topology errors and
  acceptable non-orthogonality/skewness for LES (non-orthogonality
  ideally below ~65 deg, skewness below ~4 for stable LES)
- Boundary patch names in `constant/polyMesh/boundary` match the names
  used in `case/0/*` (`farfield`, `wall`, `outlet`) - see README
  placeholder note if they don't. The repo's own `genericCase` example
  uses different names (`inlet`, `cylinder`) for its own geometry, so
  don't assume ours match yours either.

## 2. Species conservation

- sum(Y_i) = 1.0 at every cell, every timestep (mass fraction
  conservation) - check explicitly in post-processing, nothing in this
  case enforces it as a hard constraint at the solver level
- Freestream species (far from the body) should remain at Y_N2=0.767,
  Y_O2=0.233 with all others near zero - if freestream cells show
  significant dissociation, something is wrong with boundary conditions
  or chemistry activation

## 3. Shock structure (Tt, Tv)

- Tt should rise sharply across the bow shock, consistent with roughly
  the ~5000 K normal-shock estimate in `PARAMETERS.md` (oblique shock
  angle at -5 deg AoA will change the exact value - the normal-shock
  number is a sanity-check order of magnitude, not a target)
- Tv should visibly lag Tt immediately behind the shock, then relax
  toward Tt further downstream (V-T relaxation, `constant/thermo2TModel`,
  MillikanWhitePark model) - if Tv jumps to match Tt instantly,
  relaxation is not being resolved (check grid resolution in the shock-
  normal direction)
- Density should show a clear jump from rho_inf = 0.004071 kg/m^3 to a
  higher post-shock value

## 4. Dissociation front

- Y_N2, Y_O2 should decrease behind the shock; Y_N, Y_O should increase
  correspondingly
- The dissociation front should track with the Tv field (Tv-controlled
  rate via the ParkTTv vibration-chemistry coupling in
  `constant/chemistryProperties`), not instantaneously with Tt

## 5. Ionization (sanity check, not expected to be large)

- Y_N2+, Y_O2+, Y_NO+, Y_N+, Y_O+, Y_e- should remain small through most
  of the domain, consistent with the ~5000 K post-shock estimate being
  below typical strong-ionization onset (~8000-10000 K)
- If ion/electron mass fractions are unexpectedly large somewhere, check
  local Tt in that region before assuming a chemistry bug

## 6. ILES-specific checks (no explicit SGS model)

- `constant/turbulenceProperties` is `simulationType laminar;` - there is
  no eddy-viscosity term added anywhere. Resolved turbulent kinetic energy
  should not artificially collapse to near-zero away from walls/shocks -
  if it does, the KNP+Minmod/vanLeer scheme combination in
  `system/fvSchemes` is
  too dissipative for this mesh resolution and a less diffusive limiter
  should be tried
- Conversely, unphysical oscillations near the bow shock or any shock-
  shock interaction region point the other way - under-dissipative for
  this mesh/flow combination
- Since there is no SGS coefficient to tune (correction #4 in
  `SETTINGS.md`), the flux/limiter choice in `fvSchemes` is the only knob
  available for either direction

## 7. Wall boundary (isothermal, 1500 K)

- Confirm `Tt` and `Tv` both reach ~1500 K at the wall in the converged
  solution (fixedValue BC) - if the near-wall solution asymptotes to a
  different value, check for BC application issues at the actual
  converted-mesh wall patch name

## 8. Solver convergence

- Chemistry solver (`chemistryType.chemistrySolver Euler2Implicit` in
  `constant/chemistryProperties`) should not be repeatedly failing/
  retrying - if it is, `initialChemicalTimeStep` (currently 1e-9 s) may
  need adjustment
- `Yi` solver (PBiCG/DILU, tolerance 1e-12) residuals should trend down,
  not plateau or diverge
- `maxCo = 0.5` in `controlDict` is a conservative starting point given
  the complex/extreme flow assumption - if the run is stable, you may be
  able to increase it for faster wall-clock turnaround

## 9. Package/reproducibility

- Confirm you edited `MESH_FILE` in `scripts/setup.sh` before running -
  it is the only remaining placeholder in this package (chemistry/thermo
  data is real, copied from the repo, not a placeholder - see README)
- Record your actual `endTime` and reference length `L` (used for
  Reynolds number) once determined - both were left as placeholders in
  this package (see README)
