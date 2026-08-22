# Modeling Choices and Rationale

This case was checked against the actual hy2Foam source
(https://github.com/ivanZanardi/hypersonicfoam, cloned and inspected
directly) after an earlier draft of this package used fabricated/assumed
dictionary structures. This file documents both the physical/numerical
choices you asked for and the corrections made once the real repo was
available - both matter for understanding why the files look the way
they do.

## Corrections made after checking the real repo

An earlier draft of this case guessed at several hy2Foam conventions
before the repo was available. Once cloned, these were found to be wrong
and have been fixed:

1. **Chemistry/thermo data was placeholder-only before; now it's real.**
   `constant/chemDicts/hTCReactionsEarth93`, `constant/thermoDEM`, and
   `constant/thermo2TModel` are copied verbatim from the repo's
   `genericCase` example - the actual Park (1993) 11-species air
   mechanism and matching thermo/relaxation data, not invented numbers.
   `constant/thermophysicalProperties` now points at them via
   `$FOAM_CASE`-relative paths that work immediately.

2. **Species field names were wrong.** The mechanism uses literal
   `N2+`, `O2+`, `NO+`, `N+`, `O+`, `e-` (with `+`/`-` characters) as
   species names - not the `N2p`/`em`-style names used in the earlier
   draft. `case/0/` field files are named to match exactly, since hy2Foam
   associates field names to species by literal string match.

3. **`thermodynamicProperties` was the wrong filename/structure.** The
   real dictionary is `constant/thermophysicalProperties`, with a
   `thermoType` block using this solver's actual type strings
   (`heRho2Thermo`, `reacting2Mixture`, `BlottnerEucken`,
   `decoupledEnergyModes`, `sensible2InternalEnergy`, `perfect2Gas`,
   `advancedSpecie`) - confirmed against the repo's own file, not the
   generic OpenFOAM `heRhoThermo`/`multiSpecies` strings used before.

4. **"Implicit LES" needed a different mechanism than expected.**
   Checking `hyStrath/src/TurbulenceModels/compressible/
   turbulentFluidThermoModels/turbulentFluidThermoModels.C` confirmed this
   build DOES register real explicit LES models (WALE, Smagorinsky, kEqn,
   dynamicKEqn, dynamicLagrangian, several DES variants) - but there is no
   "off"/"none" LES model. Setting `LESModel none;` (the earlier draft's
   approach) would fail at runtime. True ILES in this codebase means
   `simulationType laminar;` with `laminarModel Stokes;` (standard
   continuum Navier-Stokes-Fourier, no eddy-viscosity term at all) - the
   KNP+vanLeer numerical dissipation, plus the real Blottner-Eucken
   molecular transport, is then the only dissipation mechanism present.
   That combination is what "implicit" actually means here. The fabricated
   `constant/LESProperties` file has been deleted - it was never a real
   dictionary this solver reads.

5. **`fluxScheme Kurganov` already was KNP - no change needed, just
   clarified.** There's no literal `"KNP"` string accepted by the solver;
   `numerics/KNP-KT.H` in the repo confirms `fluxScheme Kurganov` directly
   implements the Kurganov-Noelle-Petrova central-upwind scheme (source
   comment: `a_pos = ap/(ap - am); //- Eq.9 KNP`).

6. **`fvSolution`'s solved-variable set was wrong.** The earlier draft
   solved `Tt`/`Tv` directly via a linear solver. The real solver's
   primary conservative variables are `rho`, `rhoU`, `rhoE`, and
   `rhoEv.*` (per-mode vibrational energy density); `Tt` and `Tv` are
   *derived* from those via the thermo package, not solved directly. The
   real solver block instead needs `(h|e).*` and `Yi` entries, confirmed
   against the repo's own `fvSolution`.

7. **`fvConstraints` doesn't apply.** The repo's README states hyStrath
   targets OpenFOAM v1706; `fvConstraints` was introduced around v1812,
   after this solver's target version. It has been removed rather than
   left in as dead weight.

8. **`controlDict` was missing a required `libs` entry.** hy2Foam needs
   `libs ("libstrathFiniteVolume.so");` loaded at runtime - confirmed
   against the repo's own case examples and added.

## KNP flux scheme + minmod (U) / vanLeer (rest) limiters

`fluxScheme Kurganov;` in `system/fvSchemes` selects the KNP central-
upwind scheme (see correction #5 above). Momentum (`U`) uses Minmod/
MinmodV - the most diffusive of the standard TVD limiters, more so than
vanLeer - while every other reconstructed/limited field (rho, T, species)
uses vanLeer. This split (tighter limiting on momentum only) is a
reasonable way to add extra stability specifically where shock-shock
interaction would hit hardest (the momentum equation) without over-
damping the density/temperature/species fields elsewhere.

## laminarModel: Stokes, not "linearViscous"

You asked for `linearViscous` instead of `Stokes`. Checked directly
against `turbulentFluidThermoModels.C` (same file that confirmed the LES
model list) - this v1706-era hyStrath build registers exactly two
laminarModel choices: `Stokes` and `Maxwell`. There is no `linearViscous`
option to select; setting it would fail at run time with an unrecognized-
model error, not silently do something close enough.

This isn't actually a substitution, though - it's the same thing under a
different name. Stokes' hypothesis IS the linear (Newtonian) viscous
stress-strain relation, tau = 2*mu*devSymm(gradU), with zero eddy
viscosity. There's no separate "linearViscous" class in this model
hierarchy because Stokes already means exactly that. `laminarModel
Stokes;` was kept as the correct (and only compiling) choice for a pure
ILES run - real molecular viscosity only, no turbulence closure term of
any kind.

## Implicit LES (ILES) only

See correction #4 above for the mechanism. Practical consequence: if you
see excessive damping of resolved turbulent structures in post-processing,
the fix is a less dissipative limiter on the affected field - there is no
SGS coefficient to tune, because there is no SGS model in this
configuration.

## Isothermal wall, 1500 K

Applied via `fixedValue` to both `Tt` and `Tv` at the wall patch (`0/Tt`,
`0/Tv`), rather than the rarefied `nonEqSmoluchowskiJumpT` slip-temperature
BC the repo's own `genericCase` example uses at its `cylinder` patch. That
BC targets near-continuum-breakdown (rarefied) flow; this case is dense
continuum (130,000 ft, y+<1 wall-resolved mesh), so a plain isothermal
`fixedValue` is the physically appropriate and simpler choice. `Tv` is set
equal to `Tt` at the wall (full vibrational accommodation assumed - a
common simplifying choice absent better data on wall vibrational
relaxation).

## Mesh classification: Hybrid, y+ < 1

You described the mesh as having prisms, tets, and 2D (surface/boundary)
elements, with y+ < 1 in the wall-normal direction - this matches "Hybrid"
(prism layers for the boundary layer, tet core for the bulk domain) rather
than an extreme boundary-layer-only refinement (which would imply
y+ << 0.1). This classification mainly matters for the numerical-scheme
tuning discussion below, since this repo's turbulence factory doesn't use
a mesh-type-dependent coefficient lookup the way a WALE-based approach
would.

## Complexity level: Complex (leaning down from "complex/extreme")

You described this as "smooth version of the spaceplane," uncertain
between complex and extreme. The settings assume **complex** (shock-shock
interaction, entropy-layer effects possible near leading edges/control
surfaces) rather than extreme (massively separated, highly unsteady) -
"smooth" suggests fewer geometric discontinuities that typically drive
extreme separation. This is a judgment call: if you see signs of extreme
unsteadiness during the run (large-scale separation, strong shock
oscillation, poor convergence), consider a more diffusive limiter and a
smaller `maxCo`.

## Why 11-species real gas at this condition

Post-shock translational temperature estimate is roughly 5000 K (normal-
shock calculation, see `PARAMETERS.md`) - high enough that O2/N2
dissociation is physically expected behind the bow shock (a single-species
ideal-gas model would meaningfully mis-predict density and entropy-layer
size here), but likely too low for significant ionization. Expect ion/
electron mass fractions to remain small through most of the domain;
double-check against local Tt if post-processing shows otherwise.

## Reaction rates, thermo curve fits, and collision data: real, not re-typed

`constant/chemDicts/hTCReactionsEarth93`, `constant/thermoDEM`,
`constant/thermo2TModel`, and `constant/transportProperties` all ship with
real, published data (Park 1993 reaction rates; Gupta/Yos/Thompson 1989/
1990 collision integrals) copied directly from the repo rather than
transcribed from memory - an 11-species, two-temperature mechanism is far
too easy to get subtly wrong by hand, and using the actual repo's own
validated data avoids that risk entirely.

## v1706 syntax audit

Beyond checking against the hyStrath repo itself, this case was re-checked
against actual OpenFOAM version history (not just this fork) for anything
that might have drifted to post-v1706 syntax:

- **`laminarModel` vs `model`**: OpenFOAM renamed the turbulence-model-
  selection key at some point after v1706 - newer versions (confirmed
  directly against a v2412 example) use `model Maxwell;` inside the
  `laminar{}` block. v1706-era code uses `laminarModel Stokes;` (matching
  the `RASModel`/`LESModel` pattern this repo's own genericCase uses).
  `constant/turbulenceProperties` already used `laminarModel`, confirmed
  correct for this version.
- **Minmod/MinmodV/vanLeer/vanLeerV**: confirmed identical scheme
  registration from OpenFOAM-4.x (2016, pre-v1706) through OpenFOAM-dev
  today - no renaming, no version risk.
- **freestreamVelocity/freestreamPressure/waveTransmissive**: confirmed
  present in OpenFOAM-4.x tutorials, unchanged since - safe for v1706.
- **`thermophysicalProperties`**: diffed byte-for-byte against the real
  repo's own file - identical.

No file changes resulted from this pass - everything checked was already
correct for v1706, but hadn't previously been verified against actual
OpenFOAM version history rather than just this repo's own conventions.
