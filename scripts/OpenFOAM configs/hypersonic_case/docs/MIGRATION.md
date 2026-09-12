# Migration: hy2Foam (hyStrath) to rhoCentralFoam (stock v2412)

This package previously targeted `hy2Foam` from the ivanZanardi/hypersonicfoam
port of the hyStrath suite, and required `libstrathFiniteVolume.so` and
`libstrathForces.so` built against v2412 by hand. It now runs on stock
OpenCFD OpenFOAM v2412 with no custom libraries.

The reasoning for the model change is in `docs/SETTINGS.md` section 1. The
short version: the peak temperature this case reaches is 4331 K, which is
below the point where air dissociation materially changes the aerodynamics of
interest, so the two-temperature 11-species reacting model is not earning its
cost or its build dependency. What does matter at 4331 K is the caloric
imperfection (cp rises 32 percent from freestream to shock layer) and the
high-temperature viscosity, and `janaf` plus a `polynomial` transport fit
capture both.

## Removed

| File | Why |
|---|---|
| `constant/transportProperties` | hyStrath rarefaction and Gupta collision data |
| `constant/chemistryProperties` | finite-rate chemistry control |
| `constant/chemDicts/hTCReactionsEarth93` | Park 1993 11-species reactions |
| `constant/thermoDEM` | per-species decoupled-energy-mode thermo |
| `constant/thermo2TModel` | vibrational-translational relaxation model |
| `constant/hTCProperties` | high-temperature chemistry control |
| `0/Tt`, `0/Tv`, `0/Tv_*` (13 files) | two-temperature fields |
| `0/N2 O2 NO N O N2+ O2+ NO+ N+ O+ e-` | 11 species mass fractions |
| `case/README.v2412.md` | described the port that is no longer used |

`0/` went from 24 fields to 3: `p`, `U`, `T`.

## Replaced

| File | Change |
|---|---|
| `constant/thermophysicalProperties` | `heRho2Thermo`/`reacting2Mixture`/`BlottnerEucken`/`decoupledEnergyModes` becomes `hePsiThermo`/`pureMixture`/`polynomial`/`janaf`/`sensibleInternalEnergy` |
| `constant/turbulenceProperties` | unchanged in content (laminar/Stokes, ILES); comments updated |
| `system/controlDict` | `hy2Foam` becomes `rhoCentralFoam`; adaptive timestep becomes fixed 45 ns; `libs` removed; collated added; function objects split out |
| `system/fvSchemes` | species and `Yi_h` div schemes removed; `reconstruct(U)` becomes `vanLeerV`; `fluxScheme Kurganov` retained |
| `system/fvSolution` | `rhoEv.*` and `Yi` solvers removed; energy variable is `e` |
| `system/decomposeParDict` | 384 and ptscotch retained; the `preservePatches` constraint on `Atmosphere` removed as inapplicable |
| `scripts/*.sh` | rewritten for the new solver, 384 ranks, collated, no `reconstructPar` |

## Added

| File | Purpose |
|---|---|
| `case/0/T` | the single static temperature |
| `case/system/monitorFunctions` | derived fields, Courant monitor, forces, wall heat flux, mean flow |
| `case/system/probesPSE` | second-mode wall probes and boundary-layer rakes |
| `case/system/surfacesShock` | shock isosurfaces and cutting planes for ParaView |
| `case/case.foam` | ParaView entry point |
| `scripts/make_pse_probes.py` | regenerates the probe file from the real wall surface |
| `scripts/submit_384.slurm` | batch submission |
| `docs/MIGRATION.md` | this file |

## Behaviour changes worth knowing

- **The timestep no longer adapts.** It is fixed at 45 ns because uniform
  probe sampling is required for the PSE FFT. The safety margin comes from
  the 0.5 mm mesh floor; see `docs/SETTINGS.md` section 8.
- **`reconstructPar` is no longer run.** Collated output is read directly by
  ParaView.
- **The forces objects use stock `libforces.so`**, not `libstrathForces.so`.
- **`writeCompression` is now off.** With collated binary writes, gzip
  serialises through one thread and costs more than it saves.
- **Validation status.** The dictionaries have been checked for structural
  correctness (balanced delimiters, resolvable includes, no stray `FoamFile`
  header in the included function files) but have not been parsed by
  OpenFOAM itself or run, because no OpenFOAM installation was available.
  `scripts/setup.sh` runs `foamDictionary -expand` on each dictionary as its
  step 5, which is the first real parse.
