# Validation

What to check, in order, and what the answer should be. Anything marked
ESTIMATE is a value derived here rather than measured, and should be replaced
by the run's own number once available.

## Before the first step

**Patch names.** `scripts/setup.sh` prints them. They must be exactly
`Inlet`, `Outlet`, `Atmosphere`, `Solid_Walls`. Every field file and every
function object refers to those names.

**Minimum cell size.** From `checkMesh`. The fixed 45 ns step needs every
cell to be at least 0.424 mm for the 0.37 CFL target. The mesh generator
floor is 0.5 mm. If `checkMesh` reports smaller, reduce `deltaT` - nothing in
the case will catch it at runtime.

**Mesh quality.** `checkMesh -allTopology -allGeometry` must pass. Watch max
non-orthogonality (below 70 is comfortable for the `corrected` schemes) and
max skewness (below 4).

**Dictionaries parse.** `scripts/setup.sh` runs `foamDictionary -expand` on
each one. The `-expand` matters: it is what resolves the three `#include`
directives in `controlDict`.

**Probe placement.** Run `scripts/make_pse_probes.py` and read its output.
It reports how many stations it could place and names any it could not. A
station it could not place is one where the bounding-box assumption was
wrong, which is exactly what it exists to find.

## First 1000 steps

**Courant number.** `postProcessing/courantMonitor/`. Expect a maximum near
0.31. Anything above about 0.6 means the run is heading for divergence and
`deltaT` must be reduced; the step is fixed, so nothing will adapt.

**Temperature bounds.** `postProcessing/fieldMinMax/`. `T` must stay inside
`[200, 6000]` K or janaf aborts. Expect the maximum to climb towards 4331 K
as the bow shock forms. A maximum approaching 6000 K is the warning that the
janaf ceiling is about to be hit.

**Pressure positivity.** From the same file. A negative `p` means the vanLeer
limiter has failed somewhere, almost always at a sharp leading edge on a
skewed cell.

**No solver divergence.** The `e` and `U` diffusion solves should converge in
a handful of sweeps. Rising iteration counts mean the diffusive stability
limit is being approached, which at a fixed step means the same remedy.

## Physics checks against the predicted values

These are the numbers to compare the run against. All are ESTIMATES computed
from the thermally perfect janaf gas at frozen composition
(`docs/PARAMETERS.md` has the working).

| Quantity | Predicted | Where to read it |
|---|---|---|
| Post-shock `T2` | 4256 K | stagnation streamline, symmetry plane |
| Stagnation `T0` | 4331 K | `fieldMinMax` maximum of `T` |
| Post-shock `p2` | 35,445 Pa | `fieldMinMax` maximum of `p` |
| `p2/p1` | 121.3 | - |
| `rho2/rho1` | 7.13 | - |
| Post-shock `u2` | 445 m/s | symmetry plane |

The run should come out **below** the predicted `T0`, not above. Real air at
4331 K has some O2 dissociation absorbing energy that this frozen gas cannot;
the frozen model is the conservative-high case for temperature. A run that
exceeds 4331 K anywhere other than a transient start-up spike indicates a
numerical problem, not physics.

**Shock stand-off.** Compare the `shockIsoP_mid` isosurface at the nose
against a Billig correlation for the nose radius. Agreement within 10-20
percent is the expectation for a blunt body at this Mach number; a large
discrepancy usually means the mesh is too coarse where the shock sits, not
that the physics is wrong.

**Shock thickness.** The KNP scheme should smear the shock over 3 to 5 cells.
Measure it on `planeSymmetryY0`. Substantially more than that means the
limiter is being triggered too broadly.

## Boundary layer

**y+ below 1.** `postProcessing/yPlus/`, over `Solid_Walls`. This is the
condition that makes the case wall-resolved and it is what the PSE profiles
depend on. Regions above 1 are regions where the stability analysis is not
trustworthy.

**Boundary-layer thickness.** From `postProcessing/pseBLRakes/`. Compare
against the Eckert reference-temperature prediction:

| x [m] | delta predicted [mm] |
|---|---|
| 1.0 | 29.4 |
| 3.0 | 50.9 |
| 10.0 | 92.9 |
| 24.0 | 143.9 |
| 47.0 | 201.3 |

These are flat-plate estimates on a three-dimensional body, so agreement to
within a factor of about 1.5 is what to expect. Use the measured values, not
these, to re-target the PSE frequency band.

**Rake coverage.** Each rake is 0.6 m long. Confirm the profile actually
reaches the freestream at every station - if the measured delta at x = 47 m
is much larger than 201 mm the rakes are too short and
`--rake-length` in `scripts/make_pse_probes.py` needs raising.

**Wall heat flux.** `postProcessing/wallHeatFlux/`. The peak is at the nose.
This is the primary TPS output and it is the check that the polynomial
kappa(T) is doing its job - a Sutherland model would report this tens of
percent low.

## Second Mack mode

**Spectra.** FFT the wall pressure from `postProcessing/pseWallWindward/`.
Expect a distinct narrow peak that shifts to lower frequency going aft,
tracking the table in `docs/PARAMETERS.md` (about 54 kHz at x = 1 m falling
to about 8 kHz at x = 47 m). A peak that does not shift with x is not the
second mode.

**Amplitude growth.** The peak amplitude should grow exponentially with x
over the unstable region. That growth rate is what the PSE N-factor is
compared against.

**Mode shape.** From the rakes at the same station: the second mode has a
characteristic wall-pressure maximum with a phase reversal across the
critical layer. If the eigenfunction does not show that structure, the peak
is something else - a first mode, an entropy-layer disturbance, or numerical
noise.

**Two-dimensionality.** Compare phase across `pseWallSpanwise` at x = 20 m
and x = 35 m. If it is coherent across the span, `beta = 0` in the PSE
configuration is justified. If not, `n_beta` must be opened up.

**Noise floor.** Confirm the broadband floor in the spectrum sits well below
the peak. If it does not, the disturbance is not resolved above the numerical
noise and the run needs a finer mesh, not more time.

## Loads

`postProcessing/forceCoeffs/`. Cd and Cl should settle after roughly one
flow-through (0.0152 s, step 337,500). If they are still drifting at the end
of the run, the run is too short for the loads to be meaningful, whatever the
stability results show.

`Aref` is 912.1 m2 and is flagged in `PARAMETERS.md` as needing confirmation
against CAD - every coefficient scales directly with it.

## Known limitations of this case

Stated plainly so they are not discovered later:

- **Frozen composition.** No dissociation. The stagnation region is the one
  place where this is a real approximation, and it makes the predicted
  temperature there too high. Do not use this case as the source for a
  stagnation-point heating number.
- **2.0 flow-throughs.** Enough to establish the flow and to give a long
  probe record, not enough for converged turbulence statistics. The
  `fieldAverage` output is a laminar mean baseflow for PSE.
- **Transition is not modelled.** This is a laminar ILES. It computes the
  baseflow and resolves the disturbance; it does not predict where transition
  occurs. That is what the PSE N-factor is for.
- **Shock isosurface thresholds are estimates.** The `|grad rho|` values of
  5.0 and 1.0 kg/m4 come from the normal-shock density jump, not from a
  measurement. Retune them on the first written frame.
- **`Aref` is unconfirmed.** See above.
